#!/usr/bin/env python
"""Usage:
    gen_mw.py <input> [--xyz <xyz>] [--stream <stream>]

Generate methanol/water mixture for DL_MESO software.

Options:
    --xyz <xyz>         Produce xyz file
    --stream <stream>   Choose log [default: sys.stdout]

pv278@cam.ac.uk, 26/07/16
"""
import numpy as np
from math import sqrt
import yaml, sys
from collections import namedtuple
from analyse_mw_data import mol_size, get_chi, aii_from_k
import dlms_lib as dlms
from docopt import docopt

NA = 6.022e23
R = 8.314
kB = 1.381e-23
AMU = 1.66e-27

subst = namedtuple("subst", ["rho", "Mm", "sol", "Hv", "kappaT"])
water = subst(1000.0, 18e-3, 47.8, 40.65e3, 4.65e-10)
meoh = subst(792.0, 32.04e-3, 29.7, 35.3e3, 1.22e-9)


def wout(s, stream=sys.stdout):
    stream.write(s)
    stream.flush()


def calc_num_molecules(V, eta):
    """Calculate number of water and meoh molecules given 
    their ratio and total volume of the box"""
    Nw = int(V * NA / (water.Mm / water.rho + 1/eta * meoh.Mm / meoh.rho))
    Nm = int(Nw / eta)
    return Nw, Nm


def method_1(kT=1.0, rc=1.0, gamma=4.5):
    """Assumptions:
    * aii = 25 kT
    * delta(a) ~chi from Flory-Huggins"""
    a_ii = 25.0
    chi = get_chi(water, meoh)
    a_ij = {}
    a_ij[("M", "M")] = [a_ii * kT, rc, gamma]
    a_ij[("W", "W")] = [a_ii * kT, rc, gamma]
    a_ij[("M", "W")] = [a_ii * kT + 3.27 * chi, rc, gamma]
    return a_ij


def method_2(mols_in_beads, rho_DPD=3.0, kT=1.0, rc=1.0, gamma=4.5):
    """Assumptions:
    * aii based on compressibility
    * aij = sqrt(aii ajj) + 3.27 chi
    """
    a_ij = {}
    aW = mols_in_beads["W"] * aii_from_k(water) / rho_DPD
    aM = mols_in_beads["M"] * aii_from_k(meoh) / rho_DPD
    delta_aMW = 3.27 * get_chi(water, meoh)
    a_ij[("W", "W")] = [aW * kT, rc, gamma]
    a_ij[("M", "M")] = [aM * kT, rc, gamma]
    a_ij[("M", "W")] = [sqrt(aW*aM) * kT + delta_aMW, rc, gamma]
    return a_ij


if __name__ == "__main__":
    args = docopt(__doc__)
    np.random.seed(1234)
    data = yaml.load(open(args["<input>"]).read())

    rho_DPD = data["dpd-density"]
    L = data["box-size"] * 1e-9
    T = data["temperature"]
    eta = data["water-meoh-ratio"]
    # DPD scales
    rc = (mol_size(water)**3 * data["waters-in-bead"] * rho_DPD)**(1./3)
    m0 = water.Mm * 1e3 * AMU
    tau = sqrt(m0 * rc**2 / (kB * T))

    V = L**3
#    Vw, Vm = V * wr, L**3 * (1-wr)
#    Nw = water.rho * Vw * NA / water.Mm
#    Nm = meoh.rho * Vm * NA / meoh.Mm
    Nw, Nm = calc_num_molecules(V, eta)
    w_in_b, m_in_b = data["waters-in-bead"], data["meohs-in-bead"]
    Nwb = Nw // w_in_b
    Nmb = Nm // m_in_b
    mwb = w_in_b * water.Mm * 1e3 * AMU
    mmb = m_in_b * meoh.Mm * 1e3 * AMU
    
    L_DPD = np.round(L / rc, 1)
    dt_DPD = data["dt"]
    Ttot = data["run-time"] * 1e-9    # in s
    Ttot_DPD = Ttot / tau
    Nsteps = int(Ttot_DPD // dt_DPD) // 1000 * 1000   # round to 1000

    wout("===== Generating water-methanol mixture ====\n")
    wout("DPD params | rc: %.2e | tau: %.2e\n" % (rc, tau))
    wout("Box size in SI/DPD units: %.2e / %.2f\n" % (L, L_DPD))
    wout("METHOD: %i\n" % data["method"])
    wout("Number of molecules | water: %i | meoh: %i\n" % (Nw, Nm))
    wout("Beads per molecule | water: %i | meoh: %i\n" % (w_in_b, m_in_b))
    wout("Number of beads | water: %i | meoh: %i\n" % (Nwb, Nmb))
    wout("Bead masses | water: %.2e | meoh: %.2e\n" % (mwb, mmb))

    # producing xyz file
    mols_in_beads = {"W": w_in_b, "M": m_in_b}
    bead_types = ["W", "M"]
    bead_pop = {"W": Nwb, "M": Nmb}

    xyz = np.vstack( (np.random.rand(Nwb, 3), np.random.rand(Nmb, 3)) )
    atom_ids = ["W"] * Nwb + ["M"] * Nmb
    masses = dict( (i, 1.0) for i in range(1, len(bead_types)+1) )
    wout("%i beads created\n" % (Nwb + Nmb))

    # interactions
    if data["method"] == 1:
        a_ij = method_1()
    elif data["method"] == 2:
        a_ij = method_2(mols_in_beads, rho_DPD=rho_DPD)
    else:
        raise NotImplementedError

    # writing into file
    fname = "CONFIG"
    dlms.save_config(fname, atom_ids, xyz)

    field_string = "bla\n\n" +\
                   dlms.species2str(bead_types, bead_pop) +\
                   dlms.inter2str(a_ij) + \
                   "close\n"

    fname = "FIELD"
    open(fname, "w").write(field_string)
    print("FIELD file saved.")
    
    dlms.gen_control(L_DPD, dt_DPD, Nsteps, thermo=500, traj_after=0)



