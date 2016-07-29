#!/usr/bin/env python3
"""Usage:
    parse_ms.py <infile> [--boxsize <L>]

[AD HOC] Parse a Materials Studio xcf file.
Create two files, one with atoms and position, the other with bonds.

Options:
    --boxsize <L>       Rescale w.r.t. box size [default: 40]

pv278@cam.ac.uk, 10/05/16
"""
import numpy as np
from docopt import docopt
import lmp_lib as ll


def parse_beads(xmlstring, bead_dict):
    """XML string has a form of lines,
    User ID: human-readable, numbering of beads starting from 1 rising by 1
    Bead ID: XML-readable, general number of a XML element, mixed with Bond ID
    return:
    * xyz matrix (N, 3)
    * bead types vector (N, 1)
    * bead IDs vector (N, 1)
    * molecule IDs, vector (N, 1)
    """
    # Examples of XML bead lines:
    # <Bead ID="60669" Mapping="60486" Parent="6" UserID="93" XYZ="0.335008132054088,0.56406164179294,0.35715898635349"
    #    Connections="60670,60672,60746" ForcefieldType="A" BeadType="388551"/>
    # <Molecule ID="6" Mapping="60486" Parent="2" Children="ci(149):149+60636" Name="Nafion_EW1244_MW11040"/>
    N = len([line for line in xmlstring if "Bead ID=" in line])
    xyz_mat = np.zeros((N, 3), dtype=float)
    bead_IDs = np.zeros(N, dtype=int)
    user_IDs = np.zeros(N, dtype=int)
    mol_IDs = np.zeros(N, dtype=int)
    Nb = {}

    cnt = 0
    for b in bead_dict.keys():
        E = [line for line in xmlstring if "<Bead ID" in line and "ForcefieldType=\"%s\"" % b in line]
        Nb[b] = len(E)
        print("Number of %s beads: %i" % (b, Nb[b]))
        for line in E:
            key = [field for field in line.split() if "XYZ=" in field][0]
            bead_ID = [field for field in line.split() if "ID=" in field][0]
            user_ID = [field for field in line.split() if "UserID=" in field][0]
            mol_ID = [field for field in line.split() if "Parent=" in field][0]
            key = np.array(key[5:-1].split(",")).astype(float)
 
            xyz_mat[cnt] = key
            bead_IDs[cnt] = int(bead_ID[4:-1])
            user_IDs[cnt] = int(user_ID[8:-1])
            mol_IDs[cnt] = int(mol_ID[8:-1])
            cnt += 1
        del E

    bead_types = []
    for b in bead_dict.keys():
        bead_types += [bead_dict[b]]*Nb[b]
    bead_types = np.array(bead_types)
    return xyz_mat, bead_types, bead_IDs, user_IDs, mol_IDs


def parse_beads2(xmlstring, bead_dict):
    """Improved function to parse beads linearly, no by bead type"""
    N = len([line for line in xmlstring if "Bead ID=" in line])
    xyz_mat = np.zeros((N, 3), dtype=float)
    bead_IDs = np.zeros(N, dtype=int)
    bead_types = np.zeros(N, dtype=int)
    user_IDs = np.zeros(N, dtype=int)
    mol_IDs = np.zeros(N, dtype=int)
    Nb = {"A": 0, "B11": 0, "C": 0, "W": 0, "G": 0}

    i = 0
    good = [line for line in xmlstring if "<Bead ID" in line]
    for line in good:
        key = [field for field in line.split() if "XYZ=" in field][0]
        bead_ID = [field for field in line.split() if "ID=" in field][0]
        bead_type = [field for field in line.split() if "ForcefieldType=" in field][0]
        user_ID = [field for field in line.split() if "UserID=" in field][0]
        mol_ID = [field for field in line.split() if "Parent=" in field][0]
        key = np.array(key[5:-1].split(",")).astype(float)

        xyz_mat[i] = key
        bead_IDs[i] = int(bead_ID[4:-1])
        bead_type = bead_type.rstrip(">")
        bead_num = bead_dict[bead_type[16:-1]]
        bead_types[i] = bead_num
        user_IDs[i] = int(user_ID[8:-1])
        mol_IDs[i] = int(mol_ID[8:-1])
        Nb[bead_type[16:-1]] += 1
        i += 1

    for k, v in Nb.items():
        print("Number of %s beads: %i" % (k, v))
    return xyz_mat, bead_types, bead_IDs, user_IDs, mol_IDs



def parse_bonds(xmlstring):
    """Extract bonds from a XML string in the form of lines,
    return (N, 4) matrix with columns:
    * bond ID, starting from 1 and rising by 1
    * bond type = 1
    * bead 1 ID
    * bead 2 ID
    """
    # Example of a XML bond line:
    # <BeadConnector ID="60670" Mapping="60486" Parent="6" Connects="60667,60669"/>
    N = len([line for line in xmlstring if "<BeadConnector" in line])
    parsed_mat = np.zeros((N, 4), dtype=int)
    cnt = 0

    good = [line for line in xmlstring if "<BeadConnector" in line]
    for line in good:
        bond_id = [field for field in line.split() if "ID=" in field][0]
        beads = [field for field in line.split() if "Connects=" in field][0]
        beads = np.array(beads[10:-3].split(",")).astype(int)

        parsed_mat[cnt, 0] = cnt + 1
        parsed_mat[cnt, 1] = 1
        parsed_mat[cnt, 2:] = beads
        cnt += 1

    return parsed_mat


def save_atoms(fname, user_IDs, mol_IDs, bead_types, xyz_mat):
    """Ad hoc function to all information on atoms into file"""
    N = len(user_IDs)
    with open(fname, "w") as f:
         for i in range(N):
             f.write("%i\t%i\t%i\t%.6e\t%.6e\t%.6e\n" % \
                    (user_IDs[i], mol_IDs[i], bead_types[i],\
                     xyz_mat[i, 0], xyz_mat[i, 1], xyz_mat[i, 2]))
    print("xyz matrix was saved into", fname)


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<infile>"]
    bead_dict = {"A": 1, "B11": 2, "C": 3, "W": 4, "G": 5}

    xmlstring = open(infile, "r").readlines()
    print("File: %s | Length: %i" % (infile, len(xmlstring)))
    
    print("Parsing atoms...")
    xyz_mat, bead_types, bead_IDs, user_IDs, mol_IDs = parse_beads2(xmlstring, bead_dict)
    xyz_mat *= float(args["--boxsize"]) * 8.14e-10       # rescale to DPD units
    print("Parsing bonds...")
    bond_mat = parse_bonds(xmlstring)

    N = xyz_mat.shape[0]
    NB = bond_mat.shape[0]
    print("Number of atoms: %i | Number of bonds: %i" % (N, NB))

    print("Cleaning data...")
    mol_IDs += -np.min(mol_IDs) + 1            # Molecule IDs to start with 1
    for ui in user_IDs[(bead_types != 4) & (bead_types != 5)]:  # replace weird user IDs by order 1..N
        np.place(bond_mat[:, 2:], bond_mat[:, 2:]==bead_IDs[ui-1], ui)

    fname_beads = "nafion_ms.xyz"
    save_atoms(fname_beads, user_IDs, mol_IDs, bead_types, xyz_mat)

    fname_bonds = "nafion_ms.bonds"
    np.savetxt(fname_bonds, bond_mat, fmt="%i")
    print("Bond matrix saved into", fname_bonds)


