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
    * bead_pop: dict of bead population, e.g. {"A": 500, "B": 500}
    Ordering: bead type, mass, charge, number, freeze or not
    """
    s = "SPECIES %i\n" % len(bead_types)
    for b in bead_types:
        s += b + "    1.0 0.0 " + str(bead_pop[b])
        if b == "E" or b == "P":   # freeze platinum or electrode beads
            s += " 1\n"
        else:
            s += " 0\n"
    return s + "\n"


def species2str2(beads):
    """
    beads: dict of bead populations, e.g. {"A": 500, "B": 500}
    Ordering on line: bead type, mass, charge, number, freeze (or not)
    """
    s = "SPECIES %i\n" % len(beads)
    for bt, bp in beads.items():
        s += "%s    1.0 0.0 %i" % (bt, bp)
        if bt in ["E", "P"]:   # freeze platinum or electrode beads
            s += " 1\n"
        else:
            s += " 0\n"
    return s + "\n"


def inter2str(a_ij):
    """Dictionary of all pair interactions.
    * key: pair, e.g. "A B"
    * values: array of [coefficient, rc, gamma]
    Method: add mddpd in the future
    """
    s = "INTERACTIONS %i\n" % len(a_ij)
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
def save_config(fname, names, xyz, imcon=0):
    """Save positions into file
    * imcom: include box coords (0 or 1)"""
    N = len(xyz)
    conf_str = "bla\n" + "0\t%i\n" % imcon
    if imcon == 1:
        box = L*np.eye(3)
        for i in range(len(box)):
            conf_str += "%f\t%f\t%f\n" % (box[i, 0], box[i, 1], box[i, 2])

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


