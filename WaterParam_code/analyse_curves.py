#!/usr/bin/env python
"""Usage:
    analyse_curves.py <infiles> <Efile> [--rc <rc> --deg <d> --method <m>]
                      [--plot --T <T>]

Read in all the energy curves for all blob combinations
and produce a temperature-weighted master curve.
Fit the master curve by a parabola.
Two method to weigh:
1. One partition function Z per blob pair
2. One partition function per distance

Options:
    --T <T>        Temperature [default: 300]
    --plot         Plot master curve
    --fit          Fit the curve within a given distance cutoff
    --rc <rc>      Distance cutoff for fitting [default: 6e-10]
    --deg <d>      Degree of fitting polynomial [default: 2]
    --method <m>   Two ways to calculate Z [default: 1]

pv278@cam.ac.uk, 03/04/16
"""
import numpy as np
import matplotlib.pyplot as plt
import glob, sys
from docopt import docopt

kB = 1.38e-23
H2J = 1.602e-19 * 27.211     # Hartrees to Joules
AA2m = 1e-10                 # AA to metres


def get_master_curve(infiles, E, T):
    """Produce master curve from list of curves, energies of separate clusters
    and temperature. Arguments:
    * infiles: list of files storing (N, 2) matrix of [dist (AA), energy (H)]
    * E: (N, 2) array of [blob number, energy]
    * T: temperature
    """
    ME = np.asarray(np.loadtxt(infiles[0]))
    ME[:, 1] = 0.0
    E0 = 2 * np.min(E[:, 1])
    Z = 0.0                  # normalisation, i.e. partition function

    for i in range(len(infiles)):
        b1, b2 = [int(n) for n in infiles[i].rstrip(".out").split("_")[1:]]
        Ebp = np.loadtxt(infiles[i])
        # rescale wrt energy minimum to prevent underflow
        scaledE = (E[b1, 1] + E[b2, 1] - E0) * H2J
        Boltzmann = np.exp(- scaledE / (kB*T))
        print(Ebp[:, 1] - (E[b1, 1] + E[b2, 1]))
        ME[:, 1] += (Ebp[:, 1] - (E[b1, 1] + E[b2, 1])) * Boltzmann
        Z += Boltzmann
    
    ME[:, 1] /= Z
    return ME


def get_master_curve2(infiles, E, T):
    """Compared with 1st func, partition function is different 
    for each blob distance.
    Arguments:
    * infiles: list of files storing (N, 2) matrix of [dist (AA), energy (H)]
    * E: (N, 2) array of [blob number, energy]
    * T: temperature
    """
    Nb = len(E)              # blobs
    Nbp = len(infiles)       # blob pairs
    Nr = len(open(infiles[0]).readlines())
    ME = np.asarray(np.loadtxt(infiles[0]))
    ME[:, 1] = 0.0
    E0 = 2 * np.min(E[:, 1])
    Z = np.zeros(Nr)        # normalisation, i.e. partition function

    for i in range(Nbp):
        b1, b2 = [int(n) for n in infiles[i].rstrip(".out").split("_")[1:]]
        Ebp = np.loadtxt(infiles[i])
        Boltzmann = np.exp(- (Ebp[:, 1] - E0) * H2J / (kB*T))
        print(Ebp[:, 1] - (E[b1, 1] + E[b2, 1]))
        ME[:, 1] += (Ebp[:, 1] - (E[b1, 1] + E[b2, 1])) * Boltzmann
        Z += Boltzmann
    
    ME[:, 1] /= Z
    return ME


def fit_energy(ME, rc, degree=2):
    """
    Fit master curve by a quadratic function
    Arguments:
    * ME: master curve (n, 2) matrix
    * rc: discard all values beyond this cutoff
    """
    ME = ME[ME[:, 0] <= rc]    # take points less than rc
    dr = ME[1, 0] - ME[0, 0]
    N = len(ME)

    # recreate the right side of the parabola for correct fit
    # INCLUDE ALSO STEPS IN BETWEEN
    new_c = np.zeros((N+N-1, 2))
    new_c[:N] = ME
    new_c [N:, 1] = ME[::-1][1:, 1]
    new_c[N:, 0] = np.arange(ME[N-1, 0]+dr, ME[N-1, 0]+dr+(N-1)*dr, dr)

    coeffs, err, _, _, _ = np.polyfit(new_c[:, 0], new_c[:, 1], \
                                      degree, full=True)
    return coeffs, err, new_c


def quadratic_fit(A, rc):
    return A * (1 - r/rc)**2


def quartic_fit(A1, A2, rc):
    return A1 * (1 - r/rc)**4 + A2 * (1 - r/rc)**2


def gen_log(coeffs="fff", deg=1.0, kT=1.0):
    logname = "fit.log"
    s = "Fitting of energy curve\n"
    s += "=======================\n"
    s += "Degree of fit: %i\n" % deg
    s += "Coeffs: %s\n" % coeffs
    s += "kT: %.4e\n" % kT
    open(logname, "w").write(s)
    print("Log saved in %s." % logname)


if __name__ == "__main__":
    args = docopt(__doc__)
    infiles = glob.glob(args["<infiles>"])
    energies = np.loadtxt(args["<Efile>"])
    T = float(args["--T"])
    method = int(args["--method"])
    if method not in [1, 2]:
        sys.exit("Choose method 1 or 2.")
    
    print("====== Fitting water blob energies =====")
    if not infiles:
        sys.exit("No files captured, aborting.")
    print(len(infiles), "energy curves will be analysed.")

    if method == 1:
        master_curve = get_master_curve(infiles, energies, T)
    if method == 2:
        master_curve = get_master_curve2(infiles, energies, T)
    print(master_curve)
    master_curve[:, 0] *= AA2m
    master_curve[:, 1] *= H2J
    
    outname = "master_curve.out"
    np.savetxt(outname, master_curve)
    print("Master curve saved in", outname)
    
    # fitting
    rc = float(args["--rc"])
    if rc <= 4.0e-10:
        sys.exit("Too small cutoff rc.")
    if rc > 20e-10:
        print("WARNING: Radius larger than 20 AA!")
    deg = int(args["--deg"])
    if deg % 2 != 0:
        sys.exit("Only polynomial degreee allowed.")

    print("Fitting curve, rc: %.2e" % rc)
    coeffs, err, new_c = fit_energy(master_curve, rc, degree=deg)
    print("Fit:", coeffs)
    print("Fit error:", err)
    print("a(DPD) = %.2f" % (2 * coeffs[-3] * rc**2 / (kB*T)))
#    gen_log(coeffs=coeffs, deg=deg, kT=kB*T)

    fname = "fit_rc%.2f.out" % (rc / AA2m)
    np.savetxt(fname, new_c)
    print("Master curve and fit saved in %s." % fname)

    fname = "coeffs_d%i_rc%.2f.out" % (deg, rc / AA2m)
    np.savetxt(fname, np.vstack((np.arange(deg, -1, -1), coeffs)).T, fmt="%.4e")
    print("Fit coefficients saved in %s." % fname)

    if args["--plot"]:
        plt.plot(master_curve[:, 0], master_curve[:, 1])
        plt.plot(new_c[:, 0], new_c[:, 1])
        x_fit = np.linspace(4e-10, 8e-10, 101)
        plt.plot(x_fit, np.polyval(coeffs, x_fit), "r-")

        plt.xlim([4e-10, 8e-10])
        plt.ylim([np.min(master_curve[:, 1]), np.max(master_curve[:, 1])])
#        plt.ylim([-4e-15, -3.99e-15])
        plotname = "plot_master_curve.png"
        plt.savefig(plotname)
        print("Plot of master curve saved in", plotname)
 

