#!usr/bin/env python
"""
Collection of functions to produce strings in the DL_MESO formatting
style.


pv278@cam.ac.uk, 17/06/16
"""
import numpy as np


# ===== FIELD file
def species2str(bead_types, bead_pop):
    """
    * bead_types: e.g. "ABCWEP"
    * bead_pop: dict of bead population
    """
    s = "SPECIES %i\n" % len(bead_types)
    for b in bead_types:
        s += b + "    1.0 0.0 " + str(bead_pop[b])
        if b == "E" or b == "P":
            s += " 1\n"
        else:
            s += " 0\n"
    return s + "\n"


def inter2str(a_ij, method="dpd"):
    s = "INTERACTIONS %i\n" % len(a_ij)
    for k, v in a_ij.items():
        if isinstance(k, tuple):
            s += "%s %s    %s   %.3f  %.1f  %.2f\n" % \
                 (k[0], k[1], method, v[0], v[1], v[2])
        else:
            s += "%s    %s   %.3f  %.1f  %.2f\n" % (k, method, v[0], v[1], v[2])
    return s + "\n"


def mol2str(molname, Nmols, bead_list, bond_mat, bond_type="harm", k0=4.0, r0=0.1):
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
def save_config(fname, names, xyz, imcon=0):
    """Save positions into file
    * imcom: include box coords (0 or 1)"""
    N = len(xyz)
    conf_str = "pokus\n" + "0\t%i\n" % imcon
    if imcon == 1:
        box = L*np.eye(3)
        for i in range(len(box)):
            conf_str += "%f\t%f\t%f\n" % (box[i, 0], box[i, 1], box[i, 2])

    for i in range(N):
        conf_str += "%s        %i\n" % (names[i], i+1)
        # careful about the spaces
        conf_str += "    %.10f    %.10f    %.10f\n" % (xyz[i, 0], xyz[i, 1], xyz[i, 2])

    open(fname, "w").write(conf_str)
    print("Initial configuration saved in %s" % fname)


