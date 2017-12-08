#!/usr/bin/env python
"""Usage:
    gen_mdpd_ptfe_water.py <input> [--xyz]

Generate Nafion input files to feed into a DPD software suite
with Nafion in left half and pure water in right half
in the framework of MDPD.
Read input from input.yaml file.

Options:
    --xyz     Save coordinates into xyz file

pv278@cam.ac.uk, 07/12/17
"""
import numpy as np
from math import pi, cos, sin, acos
import sys
import yaml
from docopt import docopt
import lmp_lib as ll
import dlms_lib as dlms

NA = 6.022e23
AMU = 1.66e-27
r_DPD = 8.14e-10               # DPD distance unit
k0 = 4.0
m0 = 6 * 18 * AMU
kT = 300.0 * 1.381e-23
tau = np.sqrt(m0 * r_DPD**2 / kT)


def calc_nc_nw(N, Nmc, Nbm, lmbda):
    """Return number of water beads and chains based on 
    universal DPD density, number of beads and water uptake"""
    const = (lmbda - 3) / 6.0
    Nc = N / (Nmc * (const + Nbm))
    Nw = const * Nmc * Nc
    return int(Nc), int(Nw)


def grow_polymer(Nbm, Nc, Nmc, L, Ls, mu):
    """Generate xyz matrix from a given number of chains,
    by stacking one bead next to another in distance mu"""
    Nbc = Nbm * Nmc       # num beads per chain
    xyz = np.zeros((Nc * Nbc, 3))
    for i in range(Nc):
        xyz[i*Nbc : (i+1)*Nbc] = grow_one_chain(Nbm, Nmc, L, mu)
    xyz[:, 0] *= Ls / L        # fit into proper volume
    return xyz


def grow_one_chain(Nbm, Nmc, L, mu):
    """Return xyz matrix of one polymer chain.
    Return (Nbc*len(beads_in_m), 3) xyz matrix
    """
    N = Nbm * Nmc   # beads per chain
    xyz = np.zeros((N, 3))
    xyz[0] = np.random.rand(3) * L
    for i in range(1, N):
        th = acos(1 - 2 * np.random.rand())
        phi = 2 * pi * np.random.rand()
        r = mu
        new_pos = np.array([sin(th) * cos(phi), sin(th) * sin(phi), cos(th)])
        xyz[i] = xyz[i-1] + r * new_pos
        xyz[i] = np.where(xyz[i] > L, L, xyz[i])   # on boundary coord = L/0
        xyz[i] = np.where(xyz[i] < 0.0, 0.0, xyz[i])
    return xyz


def gen_water_beads(Nw, L, Ls):
    """Generate xyz matrix from a given number of water beads.
    Return (Nw, 3) xyz matrix
    Nw: number of water beads"""
    if Nw == 0:
        return np.empty((0, 3))
    xyz = np.zeros((Nw, 3))
    xyz = np.random.rand(Nw, 3)
    xyz[:, 0] *= Ls
    xyz[:, 1:3] *= L
    return xyz


def gen_water_slab(Nws, L, Lws, Ls):
    """Return (Nsup, 3) xyz matrix"""
    xyz = np.random.rand(Nws, 3)
    xyz[:, 0] *= Lws
    xyz[:, 0] += Ls
    xyz[:, 1:3] *= L
    return xyz


def gen_bonds(Nmc, Nc, mono_beads, start=0):
    """
    Generate bonds for one chain.
    * mono_beads: e.g. "AAABC", length Nbm (Number of beads per monomer)
    * Nmc: num monomers in a chain
    * Nc: number of chains
    Return (Nmc*Nbm*Nc-1, 2) matrix: [bead1, bead2]
    """
    Nbc = Nmc * len(mono_beads)   # beads per chain
    bond_mat = np.zeros(((Nbm*Nmc-1)*Nc, 2))
    Nboc = Nbc - 1                # bonds per chain
    for i in range(Nc):
        bond_mat[Nboc*i : Nboc*(i+1)] =\
                 gen_bonds_one_chain(Nmc, mono_beads, start=start+i*Nbc)
    return bond_mat


def gen_bonds_one_chain(Nmc, mono_beads, start=0):
    """
    Generate bonds for one chain.
    * mono_beads: e.g. "AAABC", length Nbm (Number of beads per monomer)
    * Nmc: num monomers in a chain
    Return (Nmc * Nbm, 2) matrix: [bead1, bead2]
    Build bond block purely from string of beads, where 1st bead is always
    connected to -2nd bead from previous monomer
    """
    Nbm = len(mono_beads)  # number of beads per monomer
    mono_bond_block = np.c_[np.arange(1, Nbm+1), np.arange(0, Nbm)]
    mono_bond_block[0, 1] = -2
    bond_mat = np.zeros((Nbm*Nmc, 2), dtype=int)
    for i in range(Nmc):
        bond_mat[Nbm*i : Nbm*(i+1)] = mono_bond_block + Nbm*i
    bond_mat += start * np.ones((len(bond_mat), 2), dtype=int)
    return bond_mat[1:]    # reject 1st dummy bond


def gen_inter_coeffs(atoms_yaml, bead_types, a_DPD, b_DPD, chi2da, gamma=4.5, rc=1.0):
    """
    Generate atomic params a_ij for all possible combinations 
    given number of atom types. Read custom bonds from input.yaml file
    * bead_types, e.g. "ABCW"
    """
    Nbt = len(bead_types)
    a_ij = {}
    for i in range(Nbt):
        for j in range(i+1):
            pair = "%s %s" % (bead_types[j], bead_types[i])
            if pair in atoms_yaml.keys():
                a_ij[pair] = [a_DPD + chi2da * atoms_yaml[pair], b_DPD, \
                        0, 0, 0, rc, gamma]
            elif pair[::-1] in atoms_yaml.keys():
                a_ij[pair] = [a_DPD + chi2da * atoms_yaml[pair[::-1]], b_DPD, \
                        0, 0, 0, rc, gamma]
            else:
                a_ij[pair] = [a_DPD, b_DPD, 0, 0, 0, rc, gamma]
    return a_ij


if __name__ == "__main__":
    args = docopt(__doc__)
    try:
        data = yaml.load(open(args["<input>"]))
    except IOError:
        sys.exit("Input file not found:", args["<input>"])
    np.random.seed(data["seed"])
    a_DPD, b_DPD = data["A"], data["B"]
    rho_DPD = data["density"]
    chi2da = data["chi2da"]
    gamma = data["gamma"] 
    r0 = data["equilibrium-dist"]
    k0 = data["bond-coeff"]
    L = data["box-size"]
    lmbda = data["water-uptake"]
    if lmbda < 3:
        sys.exit("Water uptake should be more than 3, aborting.")

    print("=== Creating input file for Nafion for DL_MESO ===")
    print("tau: %.2e | Box size: %.1f | lmbda: %.1f" % \
         (tau, L, lmbda))
    print("A: %.2f | B: %.2f | density: %.2f" % (a_DPD, b_DPD, rho_DPD))

    # ===== setting numbers
    mono_beads = data["mono-beads"]
    Nbm = len(mono_beads)
    Nmc = data["polymerisation"]
    Ls, Lws = L/2, L/2       # slab and water slab width
    N = int(rho_DPD * L**2 * Ls)   # num. poly beads
    Nc, Nw = calc_nc_nw(N, Nmc, Nbm, lmbda)
    Nws = int(rho_DPD * L**2 * Lws)  # num. beads in water slab
    Nbc = Nbm * Nmc
    print("Nafion slab size: %.2f | Water slab size: %.2f" % (Ls, Lws))
    print("Monomers per chain: %i, Beads per monomer: %i" % (Nmc, Nbm))
    print("%i polymer chains created | Water beads: %i" % (Nc, Nw))
    bead_types = "ABCWEP"
    Nbt = len(bead_types)
    bead_pop = {"A": 0, "B": 0, "C": 0, "W": Nw+Nws, "E": 0, "P": 0}

    # ===== beads inter params and bond params
    poly_xyz = grow_polymer(Nbm, Nc, Nmc, L, Ls, mu=0.5)
    wb_xyz = gen_water_beads(Nw, L, Ls)
    ws_xyz = gen_water_slab(Nws, L, Lws, Ls)
    xyz = np.r_[wb_xyz, ws_xyz, poly_xyz]
    # keep order, chains at the end
    atom_ids_l = ["W"] * Nw + ["W"] * Nws + list(mono_beads)*Nmc*Nc
    print("%i beads created, density: %.2f" % (len(xyz), len(xyz)/L**3))

    bt2num = {}
    for i, bt in enumerate(bead_types): bt2num[bt] = i+1
    atom_ids_n = [bt2num[bt] for bt in atom_ids_l]

    a_ij = gen_inter_coeffs(data["chi-params"], bead_types, \
            a_DPD, b_DPD, chi2da, gamma)
    masses = dict( (i, 1.0) for i in range(1, Nbt+1) )

    # ==== printing
    dlms.save_config("CONFIG", atom_ids_l, xyz, L * np.eye(3))

    bond_mat = gen_bonds_one_chain(Nmc, mono_beads, start=0)
    print("FIELD: %i bonds in a chain" % len(bond_mat))
    bead_list = list(mono_beads * Nmc)
    nafion_mol_str = dlms.mol2str("nafion", Nc, bead_list, bond_mat, \
                                  bond_type="harm", k0=k0, r0=r0)
    field_string = "bla\n\n" +\
                   dlms.species2str(bead_pop) +\
                   dlms.inter2str_mdpd(a_ij) + \
                   "MOLECULES 1\n" + \
                   nafion_mol_str + "\n" + \
                   "close\n"

    open("FIELD", "w").write(field_string)
    print("FIELD file saved.")

    if args["--xyz"]:
        fname = "nafion.xyz"
        ll.save_xyzfile(fname, np.c_[atom_ids_n, xyz])
        print("xyz file saved in", fname)
 
 
