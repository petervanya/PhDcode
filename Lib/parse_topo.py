#!/usr/bin/env python
"""Usage:
    read_topo.py [--topo <topo>] [<Nm>]

Playing with reading topologies to create a polymer chain from 
the given topology.
Two types of beads: A1, and [A1]
* Nm -- number of monomers in chain (polymerisation)

Options:
    --topo <top>     The topology string [default: (A 3, [B 1], [C 1]) 15]
    <Nm>             Polymerisation
"""
import numpy as np
import re
from docopt import docopt


def bond_map(atom_types="ABCW"):
    """Map bond pairs to single number"""
    N = len(atom_types)
    map = {}
    cnt = 1
    for i in range(N):
        for j in range(i+1):
            map[atom_types[i] + " " + atom_types[j]] = cnt
            cnt += 1
    return map


def prelim_parsing(raw_topo):
    """Create a number of monomers and topology out of the submitted string"""
    topo_str, Nm = raw_topo.split(")")[0], int(raw_topo.split(")")[-1])
    topo_str = topo_str.lstrip("(")
    return topo_str, Nm


def parse_beads(raw_topo):
    """From a string parse the beads and numbers in one monomer"""
    topo_str, Nm = prelim_parsing(raw_topo)
    topo_entries = [s.strip(" ") for s in topo_str.split(",")]
#    print("Topology:", topo_entries)
    Nbm = sum([int(i) for i in re.findall(r"\d+", topo_str)])
    bead_list = []
    cnt = 1
    for topo in topo_entries:
         if topo[0] == "[" and topo[-1] == "]":   # sides
             typ, num = topo.strip("[]")[0], int(topo.strip("[]")[1:])
             for i in range(num):   # MAKE GENERAL -- WHAT IF NUM > 1?
                 bead_list.append([cnt, typ, "side", cnt-1])
                 cnt += 1
         else:                                    # backbone
             typ, num = topo[0], int(topo[1:])
             for i in range(num):
                 if i == 0:
                     bond = cnt - (Nbm - num + 1)
                 else:
                     bond = cnt - 1
                 bead_list.append([cnt, typ, "backbone", bond])
                 cnt += 1
    return bead_list, Nm


def gen_bead_dict(raw_topo):
    """Return dict with bead types and their numbers in one monomer"""
    topo_str, Nm = prelim_parsing(raw_topo)
    topo_entries = [s.strip(" ") for s in topo_str.split(",")]
    bead_dict = {}
    for topo in topo_entries:
         if topo[0] == "[" and topo[-1] == "]":   # sides
             bead_dict[topo.strip("[]")[0]] = int(topo.strip("[]")[1:])
         else:                                    # backbone
             bead_dict[topo[0]] = int(topo[1:])
    return bead_dict


def construct_bonds(bead_list, Nm, start_num=0):
    """From a given bead list construct bonds
    return (Nm*Nbm, 3) matrix, rows: [bond_type, atom1, atom2]"""
    Nbm = len(bead_list)
    bond_mat = np.zeros((Nm*Nbm, 3), dtype=int)
    cnt = 0
    for i in range(Nm):
        for j in range(Nbm):
            if i == 0:
                bond_type = "A A"
            else:
                bond_type = str(bead_list[j][1]) + " " + \
                            str(bead_list[ bead_list[j][3] ][1])
            try:
                bond_num = bond_map("ABCW")[bond_type]        # MAKE "ABCW" GENERAL
            except KeyError:
                bond_num = bond_map("ABCW")[bond_type[::-1]]  # "C A" -> "A C"

            bond_mat[cnt] = [bond_num, \
                             bead_list[j][0] + Nbm*i + start_num, \
                             bead_list[j][3] + Nbm*i + start_num]
            cnt += 1
    return bond_mat[1:]       # discard the first nonexistent "bond"


if __name__ == "__main__":
    args = docopt(__doc__)
    raw_topo = args["--topo"]
    bead_list, Nm = parse_beads(raw_topo)
    bead_dict = gen_bead_dict(raw_topo)
    print("Beads in monomer:\n", bead_dict)
    print("Num beads in monomer:", sum([v for v in bead_dict.itervalues()]))
    print("Testing correct parsing of topology string:\n", np.array(bead_list))
    
#    print(bond_map("ABCW"))
    bond_mat = construct_bonds(bead_list, 3, 0)
    print("Testing connectivity for 3 monomers:\n", bond_mat)
    
