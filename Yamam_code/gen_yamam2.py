#!/usr/bin/env python
"""Usage:
    dlms_yamam.py <input> (dlms | lammps) [--xyz <xyz>]

Generate Nafion input files to feed into a DPD software suite.
Options:
* LAMMPS data file, or 
* DL_MESO input CONFIG and FIELDS files.
Read input from input.yaml file.

Options:
    --xyz <xyz>    Save coordinates into xyz file

pv278@cam.ac.uk, 17/07/16
"""
import numpy as np
from numpy import pi, cos, sin
import sys
import yaml
from docopt import docopt
import lmp_lib as ll
import dlms_lib as dlms

NA = 6.022e23
AMU = 1.66e-27
rho_DPD = 3.0
a_DPD = 104.0
rc = 7.1e-10               # DPD distance unit
k0 = 100.0                  # from the Yamam paper
rho_el = {"carbon": 2000.0, "quartz": 2648.0, "silica": 2196.0}
rho_Pt = 21450.0
elem_wts = {"C": 12.01, "O": 15.99, "Si": 28.08}


def set_electrodes(data, L):
    """From yaml file extract:
    * Lcl: electrode layer width
    * Lpt: Pt layer width on electrodes
    * Nsup: number of support beads
    * Npt: number of platinum beads
    * Scheme:
        *-----------*
        |C |     |C | L-Lpt |
        |--|PTFE |--|        } L
        |Pt|     |Pt| Lpt   |
        *-----------*
         Lcl
    """
    Lcl = data["electrodes"]["width"]  # electrode layer width
    if Lcl > L/2:
        sys.exit("Electrode layer thicker than membrane, aborting.")

    if Lcl > 1e-5:
        Pt_ratio = data["electrodes"]["Pt-ratio"]
        elmat = data["electrodes"]["material"]
        Lpt = L * Pt_ratio
        Vcl = 2 * Lcl * L**2
        print("Electrodes on | Width: %.2f | Material: %s" % (Lcl, elmat))

        Nelb = int(rho_DPD * int(Vcl))     # not (4*pi/3*rc**3), must be cubes
        Nsup = int((1-Pt_ratio) * Nelb)
        Npt = int(Pt_ratio * Nelb)
        if Nsup%2 != 0: Nsup += 1
        if Npt%2 != 0: Npt += 1
        print("Electrode beads: %i | Support: %i | Pt: %i" % (Nelb, Nsup, Npt))
    else:
        Lcl, Lpt = 0.0, 0.0
        Nsup, Npt = 0, 0
        print("Electrodes off.")
    return Lcl, Lpt, Nsup, Npt


def calc_nc_nw(N, Nmc, Nbm, lmbda):
    """Return number of water beads and chains based on 
    universal DPD density, number of beads and water uptake"""
    const = (lmbda-3)/6.0
    Nc = N/(Nmc*(const + Nbm))
    Nw = const*Nmc*Nc
    return int(Nc), int(Nw)


def grow_polymer(Nbm, Nc, Nmc, L, Lcl, mu):
    """Generate xyz matrix from a given number of chains,
    by stacking one bead next to another in distance mu"""
    Nbc = Nbm*Nmc       # num beads per chain
    xyz = np.zeros((Nc*Nbc, 3))
    for i in range(Nc):
        xyz[i*Nbc : (i+1)*Nbc] = grow_one_chain(Nbm, Nmc, L, Lcl, mu)
    return xyz


def grow_one_chain(Nbm, Nmc, L, Lcl, mu=1.0):
    """Return xyz matrix of one polymer chain.
    Return (Nbc*len(beads_in_m), 3) xyz matrix
    """
    N = Nbm*Nmc   # beads per chain
    xyz = np.zeros((N, 3))
    xyz[0] = np.random.rand(3)*L
    for i in range(1, N):
        th = np.random.rand()*pi
        phi = np.random.rand()*2*pi
        r = mu
        new_bead_pos = [r*cos(th), r*sin(th)*cos(phi), r*sin(th)*sin(phi)]
        xyz[i] = xyz[i-1] + new_bead_pos
        xyz[i] = np.where(xyz[i] > L, L, xyz[i])   # on boundary coord = L/0
        xyz[i] = np.where(xyz[i] < 0.0, 0.0, xyz[i])
    xyz[:, 0] = xyz[:, 0]*(L-2*Lcl)/L + Lcl        # fit into proper volume
    return xyz


def gen_water_beads(Nw, L, Lcl):
    """Generate xyz matrix from a given number of water beads.
    Return (Nw, 3) xyz matrix
    Nw: number of water beads"""
    if Nw == 0:
        return np.empty((0, 3))
    xyz = np.zeros((Nw, 3))
    xyz = np.random.rand(Nw, 3)
    xyz[:, 0] *= L - 2*Lcl
    xyz[:, 0] += Lcl
    xyz[:, 1:3] *= L
    return xyz


def gen_el_support(Nsup, L, Lcl, Lpt):
    """Generate electrode support on both sides of the membrane,
    squeeze into Lsup in z-dir and shift by Lpt in z-dir.
    Return (Nsup, 3) xyz matrix
    Nsup: number of electrode beads"""
    Lsup = L - Lpt
    if Nsup == 0:
        return np.empty((0, 3), dtype=float)
    xyz1 = np.random.rand(Nsup/2, 3)
    xyz1[:, 0] *= Lcl
    xyz1[:, 1] *= Lsup
    xyz1[:, 1] += Lpt
    xyz1[:, 2] *= L
    xyz2 = np.random.rand(Nsup/2, 3)
    xyz2[:, 0] *= Lcl
    xyz2[:, 0] += L - Lcl
    xyz2[:, 1] *= Lsup
    xyz2[:, 1] += Lpt
    xyz2[:, 2] *= L
    return np.vstack((xyz1, xyz2))


def gen_platinum(Npt, L, Lcl, Lpt):
    """Randomly generate Pt beads in catalyst layer
    Return (Npt, 3) xyz matrix
    Npt: number of platinum beads"""
    if Npt == 0:
        return np.empty((0, 3), dtype=float)
    xyz1 = np.zeros((Npt/2, 3))
    xyz1 = np.random.rand(Npt/2, 3)
    xyz1[:, 0] *= Lcl
    xyz1[:, 1] *= Lpt
    xyz1[:, 2] *= L
    xyz2 = np.zeros((Npt - Npt/2, 3))
    xyz2 = np.random.rand(Npt - Npt/2, 3)
    xyz2[:, 0] *= Lcl
    xyz2[:, 0] += L - Lcl
    xyz2[:, 1] *= Lpt
    xyz2[:, 2] *= L
    xyz = np.vstack((xyz1, xyz2))
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
    bond_mat = np.zeros(( (Nbm*Nmc-1)*Nc, 2) )
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
    """
    Nbm = len(mono_beads)  # number of beads per monomer
    mono_bond_block = np.array([[1, -2],\
                                [2, 1],\
                                [3, 2],\
                                [4, 3],\
                                [5, 4],\
                                [6, 5]], dtype=int)
    bond_mat = np.zeros((Nbm*Nmc, 2), dtype=int)
    for i in range(Nmc):
        bond_mat[Nbm*i : Nbm*(i+1)] = mono_bond_block + Nbm*i
    bond_mat += start*np.ones((len(bond_mat), 2), dtype=int)
    return bond_mat[1:]    # reject 1st dummy bond


def gen_inter_coeffs(atoms_yaml, bead_types, a_DPD, gamma=4.5, rc=1.0):
    """
    Generate atomic params a_ij for all possible combinations 
    given number of atom types. Read custom bonds from input.yaml file
    * bead_types, e.g. "ABCW"
    * 3.27 = coeff from Groot-Warren, JCP, 1997 for rho=3
    """
    Nbt = len(bead_types)
    a_ij = {}
    for i in range(Nbt):
        for j in range(i+1):
            pair = "%s %s" % (bead_types[j], bead_types[i])
            if pair in atoms_yaml.keys():
                a_ij[pair] = [a_DPD + 3.27*atoms_yaml[pair], rc, gamma]
            elif pair[::-1] in atoms_yaml.keys():
                a_ij[pair] = [a_DPD + 3.27*atoms_yaml[pair[::-1]], rc, gamma]
            else:
                a_ij[pair] = [a_DPD, rc, gamma]
    return a_ij


if __name__ == "__main__":
    args = docopt(__doc__)
    try:
        data = yaml.load(open(args["<input>"]))
    except IOError:
        sys.exit("Input file not found:", args["<input>"])
    np.random.seed(data["seed"])
    gamma = data["gamma"] 
    r0 = data["equilibrium-dist"]
    lmbda = data["water-uptake"]
    L = data["box-size"]
    box = L * np.eye(3)

    print("=== Creating input file for Nafion for %s ===" % \
          ("DL_MESO" if args["dlms"] else "LAMMPS"))
    print("Box size: %.1f | Water uptake: %i" % (L, lmbda))

    # ===== setting numbers
    mono_beads = "AAAABC"
    Nbm = len(mono_beads)
    Nmc = data["mono-per-chain"]
    Lcl, Lpt, Nsup, Npt = set_electrodes(data, L)
    Nelb = Nsup + Npt
    N = int(rho_DPD * L**2*(L-2*Lcl))   # num. polymer beads
#    Nw = int(rho_DPD * L**2*(L-2*Lcl) * wt)
#    Nc, Nw = calc_nc_nw(N, Nmc, Nbm, lmbda)
    Nw = int(rho_DPD * L**2 * (L - 2 * Lcl) / (1 + 4 * Nbm / lmbda))
    Nc = int((N - Nw) / (Nbm * Nmc))
    Nbc = Nbm * Nmc              
    print("Monomers per chain: %i, Beads per monomer: %i" % (Nmc, Nbm))
    print("%i polymer chains created | Water beads: %i" % (Nc, Nw))
    bead_types = "ABCWEP"
    Nbt = len(bead_types)
    bead_pop = {"A": 0, "B": 0, "C": 0, "W": Nw, "E": Nsup, "P": Npt}

    # ===== beads inter params and bond params
    poly_xyz = grow_polymer(Nbm, Nc, Nmc, L, Lcl, mu=1.0)
    wb_xyz = gen_water_beads(Nw, L, Lcl)
    el_xyz = gen_el_support(Nsup, L, Lcl, Lpt)
    pt_xyz = gen_platinum(Npt, L, Lcl, Lpt)
    xyz = np.vstack((wb_xyz, el_xyz, pt_xyz, poly_xyz))  # careful about order!
    print("%i beads created, density: %.2f" % (len(xyz), len(xyz)/L**3))

    atom_ids_l = ["W"]*Nw + ["E"]*Nsup + ["P"]*Npt + list(mono_beads)*Nmc*Nc
    bt2num = {}
    for i, bt in enumerate(bead_types): bt2num[bt] = i+1
    atom_ids_n = [bt2num[bt] for bt in atom_ids_l]

    a_ij = gen_inter_coeffs(data["chi-params"], bead_types, a_DPD, gamma, rc=1.0)
    k_ij = {1: [k0, r0]}
    masses = dict( (i, 1.0) for i in range(1, Nbt+1) )

    # ==== printing
    if args["dlms"]:
        fname = "CONFIG"
        dlms.save_config(fname, atom_ids_l, xyz, box)
        print("Coordinates file saved in %s" % fname)
 
        bond_mat = gen_bonds_one_chain(Nmc, mono_beads, start=0)
        print("FIELD: %i bonds in a chain" % len(bond_mat))
        bead_list = list(mono_beads*Nmc)
        nafion_mol_str = dlms.mol2str("nafion", Nc, bead_list, bond_mat, \
                                      bond_type="harm", k0=k0, r0=r0)
        field_string = "bla\n\n" +\
                       dlms.species2str(bead_pop) +\
                       dlms.inter2str(a_ij) + \
                       "MOLECULES 1\n" + \
                       nafion_mol_str + "\n" + \
                       "close\n"
 
        fname = "FIELD"
        open(fname, "w").write(field_string)
        print("FIELD file saved in %s" % fname)

    elif args["lammps"]:
        # ===== pair coeffs
        a_ij_lmp = {}
        for k in a_ij.keys():
            k_new = " ".join([str(bt2num[i]) for i in k.split()])
            a_ij_lmp[k_new] = a_ij[k]
        for v in a_ij.values():
            v[1], v[2] = v[2], v[1]  # swap gamma/rc

        # ===== molecular ids
        Nmol = Nw + Nelb + Nc  # total num of molecules
        mol_ids = list(range(1, Nw+Nelb+1))
        for i in range(Nw+Nelb+1, Nmol+1): #(1, Nc+1):   # chains
            mol_ids += [i]*Nbc

        xyz_str = ll.atoms2str(np.hstack((np.matrix(mol_ids).T,\
                               np.matrix(atom_ids_n).T, xyz)))
        # ===== bonds
        bond_mat = gen_bonds(Nmc, Nc, mono_beads, start=Nw+Nelb)
        bonds_str = ll.bonds2str2(bond_mat)
        print("%i bonds created." % len(bond_mat))

        data_string = ll.header2str(len(xyz), len(bond_mat), Nbt, len(k_ij), L) + \
                      ll.mass2str(masses) + \
                      ll.pair_dpd_coeffs2str(a_ij_lmp) + \
                      ll.bond_coeffs2str(k_ij) + \
                      "Atoms\n\n" + xyz_str + \
                      "Bonds\n\n" + bonds_str

        fname = "nafion.data"
        open(fname, "w").write(data_string)
        print("Data file saved in", fname)

    if args["--xyz"]:
        fname = args["--xyz"]
        xyz = np.hstack((np.matrix(atom_ids_n).T, xyz))
        ll.save_xyzfile(fname, xyz)
        print("xyz file saved in", fname)
 
 
