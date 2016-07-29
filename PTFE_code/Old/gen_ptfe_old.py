#!/usr/bin/env python
"""Usage:
    gen_ptfe.py <input> [--el <el>] [--parsetopo]
                        [--save <fname>] [--xyz <xyz>] [--units <u>]

Generate LAMMPS data file from input.yaml file.
* Nbm: beads in a monomer
* Nmc: monomers in a chain
* Number of chains (Nc) to fit the density from the papers (dry or wet) or DPD rho=3
* Add water beads
* Option: add electrodes of carbon black and catalyst layer.

Arguments:
    <input>                  Input yaml file

Options:
    --parsetopo              Get bead topology from the string in input.yaml, else use default (PTFE)
    --el <el>                Add electrodes: "carbon", "quartz" or (amorphous) "silica"
    --save <fname>           Save the data file [default: nafion.data]
    --xyz <xyz>              Print as xyz file for VMD view
    --units <u>              SI or DPD [default: DPD]

pv278@cam.ac.uk, 09/11/15
"""
import numpy as np
from math import *
import sys
import yaml
from docopt import docopt
import Lib.parse_topo as pt
from Lib.nafion_lib import nafion_bead_wt, num_poly_chains
import Lib.lmp_lib as ll


kB = 1.38064e-23
NA = 6.022e23
AMU = 1.66e-27
m0 = 6*18*AMU
rc = 8.14e-10               # DPD distance unit
rho_DPD = 3.0
a_DPD = 25.0
k0 = 4.0
rho_dry = 1950.0            # kg/m^3
rho_wet = 1680.0
rho_Pt = 21450.0
rho_el = {"carbon": 2000.0, "quartz": 2648.0, "silica": 2196.0}
elem_wts = {"C": 12.01, "O": 15.99, "Si": 28.08}


def calc_nc_nw(N, Nmc, Nbm, lmbda):
    """Return number of water beads and chains based on 
    universal DPD density, number of beads and water uptake"""
    const = (lmbda-3)/6.0
    Nc = N/(Nmc*(const + Nbm))
    Nw = const*Nmc*Nc
    return int(Nc), int(Nw)


def grow_polymer(beads_in_m, Nc, Nmc, L, Lcl, mu, sigma):
    """Generate xyz matrix from a given number of chains,
    by stacking one bead next to another with N(mu, sigma)"""
    Nbm = len(beads_in_m)
    Nbc = Nbm*Nmc
    mol_ids = range(1, Nc+1)
    xyz = np.zeros((Nc*Nbc, 5))
    for i in range(Nc):
        xyz[i*Nbc : (i+1)*Nbc] = grow_one_chain(beads_in_m, Nmc, L, Lcl, mol_ids[i], mu, sigma)
    return xyz


def grow_one_chain(beads_in_m, Nmc, L, Lcl, mol_id=1, mu=1.0, sigma=0.1):
    """Return xyz matrix of one polymer chain.
    Return (Nbc*len(beads_in_m), 5) xyz matrix with columns:
    * molecule/chain ID
    * atomic ID
    * xyz coords"""
    Nbm = len(beads_in_m)
    types = np.matrix(beads_in_m*Nmc).T
    mol_ids = np.matrix([mol_id]*Nbm*Nmc).T
    xyz = np.zeros((Nbm*Nmc, 3))
    xyz[0] = np.random.rand(3)*L
    for i in range(1, len(types)):
        theta = np.random.rand()*pi
        phi = np.random.rand()*2*pi
        r = mu #+ np.random.randn()*L*sigma
        new_bead_pos = [r*cos(theta), r*sin(theta)*cos(phi), r*sin(theta)*sin(phi)]
        xyz[i] = xyz[i-1] + new_bead_pos
        xyz[i] = np.where(xyz[i] > L, L, xyz[i])     # on the boundary set coordinate to L or 0
        xyz[i] = np.where(xyz[i] < 0.0, 0.0, xyz[i])
    xyz[:, 0] = xyz[:, 0]*(L - 2*Lcl)/L + Lcl
    return np.hstack((mol_ids, types, xyz))


def bonds_mat(Na, Nmc, Nbm=5):
    """
    NOT USED NOW
    Nc -- num chains
    Nmc -- num monomers in a chain
    Nbm -- num beads in monomer
    Types:
    * 1: AA bond
    * 2: AB bond
    * 3: AC bond
    Order:
    count, id, part1, part2
    """
    Nbonds = Nc * (Nmc-1 + Nmc*4)        # between m'mers + 4 in each m'mer
    bonds = np.zeros((Nbonds, 3), dtype=int)
    cnt = 0
    for i in range(Nc):
        for j in range(Nmc):           # within a monomer
            first = Nc*i+j+1           # first bead in a monomer
            bonds[cnt]   = [1, first,   first+1]
            bonds[cnt+1] = [1, first+1, first+2]
            bonds[cnt+2] = [2, first,   first+3]
            bonds[cnt+3] = [3, first+3, first+4]
            cnt += 4
    for i in range(Nc):
        firstm = Nc*i + 1              # first monomer in str
        for j in range(Nmc-1):         # between monomers
            bonds[cnt] = [1, firstm+j*Nbm+2, firstm+j*Nbm+5]
            cnt += 1
    return bonds


def bonds_mat2(raw_topo, Nc):
    """Create bond matrix from topology string and number of chains"""
    bead_list, Nmc = pt.parse_beads(raw_topo)
    Nbm = len(bead_list)
    bonds = pt.construct_bonds(bead_list, Nmc, 0)
    for i in range(1, Nc):
        bonds = np.vstack( (bonds, pt.construct_bonds(bead_list, Nmc, i*Nmc*Nbm)) ) 
    return bonds


def gen_pair_coeffs(bead_types, atoms_yaml, gamma, units="DPD"):
    """
    Generate atomic params a_ij for all possible combinations 
    given number of atom types. Read custom bonds from input.yaml file
    Nafion bead_types = "ABCW"
    """
    a_ij = {}
    Nbt = len(bead_types)
    num2coeff = dict((num, coeff) for (num, coeff) in zip(range(1, Nbt+1), bead_types))
    for i in range(1, Nbt+1):
        for j in range(1, i+1):
            key = "%i %i" % (j, i)
            lkey = "%s %s" % (num2coeff[j], num2coeff[i])
            if lkey in atoms_yaml.keys() or lkey[::-1] in atoms_yaml.keys():
                try:
                    a_ij[key] = [(a_DPD + 3.27*atoms_yaml[lkey]), gamma, 1.0]
                except KeyError: # "B A" -> "A B"
                    a_ij[key] = [(a_DPD + 3.27*atoms_yaml[lkey[::-1]]), gamma, 1.0]
            else:
                a_ij[key] = [a_DPD, gamma, 1.0]

    if units == "SI":
        for k in a_ij.keys():
            a_ij[k][0] *= kBT/rc
            a_ij[k][1] *= m0 / (sqrt(m0 * rc**2/(kBT)))
            a_ij[k][2] *= rc
    return a_ij


def gen_bond_coeffs(bead_types, bonds_yaml, r0, units="DPD"):
    """
    Generate bond coeffs k_ij and r0 for all possible combinations
    given number of atom types. Read custom bonds from input.yaml file
    Nafion bead types = "ABCW"
    """
    k_ij = {}
    bmap = pt.bond_map(bead_types)   # "ABCW" -> 1234
    Nbt = len(bead_types)
    for i in range(Nbt):
        for j in range(i+1): 
            key = bead_types[i] + " " + bead_types[j]
            if key in bonds_yaml.keys():
                k_ij[bmap[key]] = [bonds_yaml[key], r0]
            else:
                k_ij[bmap[key]] = [k0, r0]
    
    if units == "SI":
        for k in k_ij.keys():
            k_ij[k][0] *= kBT/(rc**2)
            k_ij[k][1] *= rc
    return k_ij


def gen_water_beads(Nw, L, Lcl, count=1):
    """Generate xyz matrix from a given number of water beads.
    count -- where to start molecular id counting
    Nw: number of water beads"""
    if Nw == 0:
        return np.empty((0, 5))
    xyz = np.zeros((Nw, 5))
    xyz[:, 2:5] = np.random.rand(Nw, 3)
    xyz[:, 2] *= L - 2*Lcl
    xyz[:, 2] += Lcl
    xyz[:, 3] *= L
    xyz[:, 4] *= L
    xyz[:, 1] = 4      # atom id, MAKE THIS GENERAL
    xyz[:, 0] = range(count, count+Nw)
    return xyz


def gen_electrodes(Nelb, L, Lcl, count=1):
    """Generate electrodes on both sides of the membrane
    Nelb: number of electrode beads"""
    xyz1 = np.zeros((Nelb/2, 5))
    xyz1[:, 2:5] = np.random.rand(Nelb/2, 3)
    xyz1[:, 2] *= Lcl
    xyz1[:, 3] *= L
    xyz1[:, 4] *= L
    xyz2 = np.zeros((Nelb - Nelb/2, 5))
    xyz2[:, 2:5] = np.random.rand(Nelb - Nelb/2, 3)
    xyz2[:, 2] *= Lcl
    xyz2[:, 2] += L - Lcl
    xyz2[:, 3] *= L
    xyz2[:, 4] *= L
    xyz = np.vstack((xyz1, xyz2))
    xyz[:, 1] = 5      # atom id, MAKE THIS GENERAL
    xyz[:, 0] = range(count, count+Nelb)
    return xyz


def gen_platinum(NPt, L, Lcl, count=1):
    """Randomly generate Pt beads in catalyst layer
    NPt: number of platinum beads"""
    if NPt == 0:
        return np.empty((0, 5), dtype=float)
    xyz1 = np.zeros((NPt/2, 5))
    xyz1[:, 2:5] = np.random.rand(NPt/2, 3)
    xyz1[:, 2] = Lcl
    xyz1[:, 3] *= L
    xyz1[:, 4] *= L
    xyz2 = np.zeros((NPt - NPt/2, 5))
    xyz2[:, 2:5] = np.random.rand(NPt - NPt/2, 3)
    xyz2[:, 2] = L - Lcl
    xyz2[:, 3] *= L
    xyz2[:, 4] *= L
    xyz = np.vstack((xyz1, xyz2))
    xyz[:, 1] = 6    # atom id, MAKE THIS GENERAL
    xyz[:, 0] = range(count, count+NPt)
    return xyz


# =====
# ===== Main
# =====
if __name__ == "__main__":
    args = docopt(__doc__)
    try:
        data = yaml.load(open(args["<input>"]))
    except IOError:
        print("Input file not found:", args["<input>"])
        sys.exit()
    np.random.seed(1234)
    
    units = args["--units"]
    kBT = data["temperature"]
    L = data["box-size"]
    tau = 1.0/sqrt(kBT)
    gamma = data["gamma"] 
    r0 = data["equilibrium-dist"]
    lmbda = data["water-uptake"]

    if units == "SI":
        kBT = kB*300.0
        tau = sqrt(m0 * rc**2/(kBT))
#        gamma *= m0/tau

    if lmbda < 3:
        print("Water uptake should be more than 3, aborting.")
        sys.exit()

    print("=== Creating LAMMPS input file for Nafion ===")
    print("Box size:", L,"(DPD) | Temperature:", kBT,"| Tau:", tau)

    # ===== set electrode parameters
    if args["--el"]:
        elmat = args["--el"]
        if not elmat in ["carbon", "quartz", "silica"]:
            print("Choose electrodes from 'carbon', 'quatrz', 'silica' (amorphous). Aborting.")
            sys.exit()
        Lcl = data["electrodes"]["width"]  # catalyst layer width
        if Lcl > L/2:
            print("Catalyst layer thicker than membrane, aborting.")
            sys.exit()

        Pt_amount = data["electrodes"]["Pt-ratio"]
        Vcl = L**2 * 2*Lcl

        print("Electrodes on | CL width on both sides:", Lcl, "(DPD) | Material:", elmat)

        Nelb = int(rho_DPD * int(Vcl))         # not (4*pi/3*rc**3) ) must be cubes
        NPt = int(Pt_amount * Nelb)

        if elmat == "carbon":
            Natoms = NA * rho_el[elmat] * Vcl/(elem_wts["C"]*1e-3)
        elif elmat == "quartz":
            Natoms = NA * rho_el[elmat] * Vcl/((elem_wts["Si"] + 2*elem_wts["O"])*1e-3)
        elif elmat == "silica":
            Natoms = NA * rho_el[elmat] * Vcl/((elem_wts["Si"] + 2*elem_wts["O"])*1e-3)
        print("CL:", int(Natoms/Nelb*rc**3), "electrode atoms per bead at density", rho_el[elmat])
        print("Electrode beads: %i | Platinum beads: %i" % (Nelb, NPt))
    else:
        Lcl, Vcl = 0.0, 0.0
        print("Electrodes off.")

    # ===== set polymer parameters
    if args["--parsetopo"]:
        raw_topo = data["topology"]
        print("Topology:", raw_topo)
        bead_list, Nmc = pt.parse_beads(raw_topo)   # information about connectivity
        Nbm = len(bead_list)
        bead_dict = pt.gen_bead_dict(raw_topo)      # dict of beads in one monomer
        bead_types = sorted("".join(bead_dict.keys()) + "W")
        Nbt = len(bead_types)
        coeff2num = dict((coeff, num) for num, coeff in zip(range(1, Nbt+1), bead_types))
        beads = []
        [[beads.append(coeff2num[k]) for i in range(v)] for k, v in sorted(bead_dict.items())]
    else:
        Nbm, Nmc = 5, data["mono-per-chain"]
        bead_types = "ABCW"
        beads = [1, 1, 1, 2, 3]
        Nbt = len(bead_types)
    print("Nmc: %i, Nbm: %i" % (Nmc, Nbm))

    # ===== setting numbers
    N = int(rho_DPD * (L**3 - Vcl))            # must be cubes
    Nc, Nw = calc_nc_nw(N, Nmc, Nbm, lmbda)
    print(Nc, "polymer chains created")

    Nbc = Nbm*Nmc              
    Nb = Nbc*Nc                
    mu, sigma = 1.0, 0.01   # arbitrary sigma for generating chains

    # ===== beads
    poly_xyz = grow_polymer(beads, Nc, Nmc, L, Lcl, mu, sigma)
    wb_xyz = gen_water_beads(Nw, L, Lcl, count=Nc+1)
    if args["--el"]:
        el_xyz = gen_electrodes(Nelb, L, Lcl, count=Nc+Nw+1)
        pt_xyz = gen_platinum(NPt, L, Lcl, count=Nc+Nw+Nelb+1)
        xyz = np.vstack((poly_xyz, wb_xyz, el_xyz, pt_xyz))
        bead_types += "EP"
        Nbt = len(bead_types)
    else:
        xyz = np.vstack((poly_xyz, wb_xyz))
    print(len(xyz), "beads created, density:", float(len(xyz)) / L**3)

    masses = dict( (i, 1.0) for i in range(1, Nbt+1) )  # all beads weigh the same

    if units == "SI":
        masses = dict( (i, m0) for i in range(1, Nbt+1) )
        xyz *= rc
        L *= rc
    xyz_str = ll.atoms2str(xyz)

    # ===== bonds
    bonds = bonds_mat2(data["topology"], Nc)
    bonds_str = ll.bonds2str(bonds)
    print(len(bonds), "bonds created")

    # ===== pair and bond parameters
    a_ij = gen_pair_coeffs(bead_types, data["chi-params"], gamma, units)
    k_ij = gen_bond_coeffs(bead_types, data["bond-coeffs"], r0, units)

    # ===== putting it together
    final_string = ll.header2str(len(xyz), len(bonds), Nbt, len(k_ij), L) + \
                   ll.mass2str(masses) + \
                   ll.pair_dpd_coeffs2str(a_ij) + \
                   ll.bond_coeffs2str(k_ij) + \
                   "Atoms\n\n" + xyz_str + \
                   "Bonds\n\n" + bonds_str

    fname = args["--save"]
    open(fname, "w").write(final_string)
    print("Data file saved in", fname)

    if args["--xyz"]:
        fname = args["--xyz"]
        ll.save_xyzfile(fname, xyz[:, 1:])
        print("xyz file saved in", fname)


