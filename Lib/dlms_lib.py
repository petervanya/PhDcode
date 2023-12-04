#!usr/bin/env python
"""
Collection of functions to produce strings in the DL_MESO 
formatting style.

pv278@cam.ac.uk, 17/06/16
"""
import numpy as np
import sys
from collections import OrderedDict


# ===== manipulate output
def read_xyzfile(outfile):
    """Read one xyz outfile into a numpy matrix.
    Return one (n, 4) matrix."""
    try:
        A = open(outfile, "r").readlines()[2:]
    except FileNotFoundError:
        sys.exit("File %s not found." % outfile)
    A = [line.split() for line in A]
    A = np.array(A, order="F").astype(float)
    return A


def read_xyzfile2(outfile):
    """Read one xyz outfile into a numpy matrix.
    Return vector of names and (n, 3) xyz matrix."""
    try:
        A = open(outfile, "r").readlines()[2:]
    except FileNotFoundError:
        sys.exit("File %s not found." % outfile)
    A = open(outfile, "r").readlines()[2:]
    A = np.array([line.split() for line in A]).astype(float)
    names, xyz = A[:, 0].astype(int), A[:, 1:4]
    return names, xyz


def save_xyzfile(fname, mat):
    """Take xyz matrix [ids, x, y, z] and save into fname"""
    N = len(mat)
    with open(fname, "w") as f:
        f.write(str(N) + "\nbla\n")
        for i in range(N):
            f.write("%i\t%f\t%f\t%f\n" % \
                    (mat[i, 0], mat[i, 1], mat[i, 2], mat[i, 3]))


# ===== FIELD file
def species2str(beads):
    """
    beads: dict of bead populations, e.g. {"A": 500, "B": 500}
    Ordering on line: bead type, mass, charge, number, freeze (or not)
    """
    s = "SPECIES %i\n" % len(beads)
    beads = OrderedDict(sorted(beads.items()))
    for bt, bp in beads.items():
        s += "%s    1.0 0.0 %i" % (bt, bp)
        if bt in ["E", "P"]:   # freeze platinum or electrode beads
            s += " 1\n"
        else:
            s += " 0\n"
    return s + "\n"


def species2str2(beads):
    """
    Include freezing parameter.
    beads: dict of bead populations, e.g. {"A": [500, 0], "B": [500, 0]}
    Ordering on line: bead type, mass, charge, population, freeze (or not)
    """
    s = "SPECIES %i\n" % len(beads)
    beads = OrderedDict(sorted(beads.items()))
    for bt, v in beads.items():
        s += "%s    1.0 0.0 %6i %i\n" % (bt, v[0], v[1])
    return s + "\n"


def inter2str(a_ij):
    """Dictionary of all pair interactions.
    * key: pair, e.g. "A B"
    * values: array of [coefficient, rc, gamma]
    Method: add mddpd in the future
    """
    s = "INTERACTIONS %i\n" % len(a_ij)
    a_ij = OrderedDict(sorted(a_ij.items()))
    for k, v in a_ij.items():
        if isinstance(k, tuple):
            s += "%s %s    %s   %.3f  %.1f  %.2f\n" \
                % (k[0], k[1], "dpd", v[0], v[1], v[2])
        else:
            s += "%s    %s   %.3f  %.1f  %.2f\n" \
                % (k, "dpd", v[0], v[1], v[2])
    return s + "\n"


def inter2str_mdpd(a_ij):
    """Dictionary of all pair interactions for many-body DPD.
    * key: pair, e.g. "A B"
    * values: array of [coefficient, rc, gamma]
    Method: add mddpd in the future
    """
    s = "INTERACTIONS %i\n" % len(a_ij)
    a_ij = OrderedDict(sorted(a_ij.items()))
    for k, v in a_ij.items():
        if isinstance(k, tuple): # MODIFY!
            s += "%s %s    %s   " % (k[0], k[1], "mdpd")
            s += "%.3f  %.3f  %.1f  %.1f  %.1f  %.1f  %.2f\n" \
                % (v[0], v[1], v[2], v[3], v[4], v[5], v[6])
        else:
            s += "%s    %s   " % (k, "mdpd")
            s += "%.3f  %.3f  %.1f  %.1f  %.1f  %.1f  %.2f\n" \
                % (v[0], v[1], v[2], v[3], v[4], v[5], v[6])
    return s + "\n"


def chi2da_mdpd(A, B):
    """Coefficient converting chi to da for rd = 0.75"""
    rho = 2.901 + 0.68 * (-A) * (B - 4.09)**(-0.716)
    return - 0.259 + 0.196 * rho 


def mol2str(molname, Nmols, bead_list, bond_mat, bond_type="harm", \
    k0=4.0, r0=0.1):
    """Input:
    * molecule name
    * Nmols: number of molecules of this type
    * bead_list: list of bead types in one molecule, each one char
    * bond_mat: (Nbonds, 2) matrix of connected beads
    * bond_type
    * k0: spring constant
    * r0: equilibrium distance
    """
    s = molname + "\n"
    s += "nummols %s \n" % str(Nmols)
    s += "beads %i\n" % len(bead_list)
    for n in bead_list:
        s += n + "\n"
    s += "bonds %i\n" % len(bond_mat)
    for i in range(len(bond_mat)):
        s += "%s  %.2i %.2i %.3f %.3f\n" % \
             (bond_type, bond_mat[i, 0], bond_mat[i, 1], k0, r0)
    s += "finish\n"
    return s


# ===== CONFIG file
def read_config(conffile, nafion=False, shift=False):
    """
    Read CONFIG file. 

    Options
    =======
    * nafion: use bead ordering 'ABCWEP'
    * shift: shift all beads by half-box
    """
    try:
        conf_str = open(conffile, "r").readlines()
    except FileNotFoundError:
        sys.exit("File not found: %s." % conffile)

    levcfg, imcon = (int(s) for s in conf_str[1].split())
    if imcon == 0:     # no box coordinates
        conf_str = np.array(conf_str[2:])
    else:              # skip over box coordinates
        box = np.array([line.split() for line in conf_str[2:5]]).astype(float)
        L = np.diag(box)
        print("System size:", L)
        conf_str = np.array(conf_str[5:])

    print("Levcfg = %i" % levcfg)
    if levcfg == 0:
        N = len(conf_str) // 2
    elif levcfg == 1:
        print("File contains velocities.")
        N = len(conf_str) // 3
    elif levcfg >= 2:
        print("File contains velocities and forces.")
        N = len(conf_str) // 4
    
    mask = np.arange(N) * (levcfg + 2)
    names = [conf_str[i].split()[0] for i in mask]
    names_dict = {}

    if nafion:
        print("Nafion bead order 'ABCWEP' used.")
        for i, bt in enumerate("ABCWEP"):
            names_dict[bt] = i + 1
    else:
        for i in enumerate(sorted(set(names))):   # bead types to numbers
            names_dict[i[1]] = i[0] + 1

    names = [names_dict[i] for i in names]

    xyz = np.array([[float(j) for j in conf_str[i].split()] for i in mask+1])

    if shift:
        if imcon == 0:
            L = np.array([np.round((max(xyz[:, i]) - min(xyz[:, i])), 1) \
                    for i in range(3)])
        xyz += L
        xyz = xyz % L
        print("Shifted box by", L / 2)

    vel = None
    forces = None

    if levcfg == 1:
        vel = np.array([[float(j) for j in conf_str[i].split()] \
                for i in mask+2])
        forces = None

    if levcfg >= 2:
        vel = np.array([[float(j) for j in conf_str[i].split()] \
                for i in mask+2])
        forces = np.array([[float(j) for j in conf_str[i].split()] \
                for i in mask+3])

    return names, xyz, vel, forces


def save_config(fname, names, xyz, box=None, nafion=False):
    """
    Save positions into file.

    Options
    =======
    * box: matrix of form diag(Lx, Ly, Lz)
    * nafion: use bead ordering 'ABCWEP'
    """
    N = len(xyz)
    imcon = 0 if box is None else 1
    conf_str = "bla\n" + "0\t%i\n" % imcon

    if box is not None:
        for i in range(len(box)):
            conf_str += "%f\t%f\t%f\n" % (box[i, 0], box[i, 1], box[i, 2])

    names_dict = {}
    if nafion:
        print("Nafion bead order 'ABCWEP' used.")
        for i, bt in enumerate("ABCWEP"):
            names_dict[i+1] = bt
        names = [names_dict[i] for i in names]

    for i in range(N):
        conf_str += "%s        %i\n" % (names[i], i+1)
        # careful about the spaces!
        conf_str += "    %.10f    %.10f    %.10f\n" \
            % (xyz[i, 0], xyz[i, 1], xyz[i, 2])

    open(fname, "w").write(conf_str)
    print("Initial configuration saved in %s." % fname)


# ===== CONTROL file
def gen_control(L, dt, steps, thermo=100, halo=2.5, traj_after=20000, \
    kT=1.0, method="dpd", rd=0.75):
    s = "bla\n\n"
    
    s += "volume %.2f\n" % L**3
    s += "temperature 1.0\n"
    s += "cutoff 1.0\n"
    if method == "mdpd":
        s += "manybody cutoff %.2f\n" % rd
    s += "boundary halo %.1f\n\n" % halo
    
    s += "timestep %.3f\n" % dt
    s += "steps %i\n" % steps
    s += "equilibration steps 0\n"
    s += "scale temperature every 10\n"
    s += "trajectory %i 100\n" % traj_after
    s += "stats every 100\n"
    s += "stack size 100\n"
    s += "print every %i\n\n" % thermo
    
    s += "job time 1000000.0\n"
    s += "close time 1.0\n\n"
    
    s += "ensemble nvt dpdvv\n\n"
    s += "finish\n"
    
    print("CONTROL: Box size: %.1f | Timestep: %.3f | Num steps: %i" % \
         (L, dt, steps))
    open("CONTROL", "w").write(s)
    print("CONTROL file saved.")


