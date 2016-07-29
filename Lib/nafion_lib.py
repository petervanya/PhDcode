#!/usr/bin/env python
"""
A collection of functions regarding Nafion DPD simulation

pv278@cam.ac.uk, 13/05/16
"""
import sys, yaml

elem_wts = yaml.load(open(sys.path[0]+"/atomic_weights.yaml").read())
AMU = 1.66e-27


def nafion_bead_wt(bead, arg="wet"):
    """Return weigths in atomic units of beads A, B, C, W, E, P"""
    if bead == "A":          # (CF2)6
        return 6*elem_wts["C"] + 12*elem_wts["F"]
    elif bead == "B":        # O CF2 C(CF3)F O CF2
        return 2*elem_wts["O"] + 4*elem_wts["C"] + 8*elem_wts["F"]
    elif bead == "C":
        if arg == "dry":   # CF3 SO3H
            return elem_wts["C"] + 3*elem_wts["F"] + elem_wts["S"] + \
                   3*elem_wts["O"] + elem_wts["H"]
        else:              # CF3 SO3H 3H2O
            return elem_wts["C"] + 3*elem_wts["F"] + elem_wts["S"] + \
                   6*elem_wts["O"] + 7*elem_wts["H"]
    elif bead == "W":        # (H2O)6
        return 12*elem_wts["H"] + 6*elem_wts["O"]
    elif bead == "E":        # carbon black
        return 54 * elem_wts["C"]
    elif bead == "P":        # Pt
        return 30 * elem_wts["Pt"]
    else:
        print("No such bead.")
        return 0


def num_poly_chains(rho, V, Nmc=15, arg="wet"):
    """Calculate number of chains in the given volume
    and using density from papers, using SI units
    TODO:
    * customise the numbers of beads, not just for Nafion"""
    tot_mass = rho * V  
    one_chain_mass = (3*nafion_bead_wt("A") + nafion_bead_wt("B") + \
                     nafion_bead_wt("C")) * Nmc * AMU
    return tot_mass/one_chain_mass


