#!/usr/bin/env python
"""Usage:
    avg_mixture.py <sol1> <nmol1> <sol2> <nmol2> <frame>

Average xyz frame of binary mixtures to their centre of mass.

pv278@cam.ac.uk, 21/06/17
"""
import numpy as np
import sys
from docopt import docopt
from xyz_lib import Atoms

masses = {}
masses["SPC"] = {1: 1.008, 2: 15.999}
masses["TIP3P_EW"] = {1: 1.008, 2: 15.999}
masses["TIP4P_2005"] = {1: 1.008, 2: 15.999}
masses["Methanol"] = {1: 1.008, 2: 1.008, 3: 12.011, 4: 15.999}
masses["Ethanol"] = {1: 1.008, 2: 1.008, 3: 12.011, 4: 12.011, 5: 15.999}
masses["Acetone"] = {1: 12.011, 2: 12.011, 3: 1.008, 4: 15.999}
masses["SPC_T320"] = {1: 1.008, 2: 15.999}
masses["SPC_T340"] = {1: 1.008, 2: 15.999}
masses["SPC_T360"] = {1: 1.008, 2: 15.999}

atoms = {}
atoms["SPC"] = [2, 1, 1]
atoms["TIP3P_EW"] = [2, 1, 1]
atoms["Methanol"] = [3, 4, 1, 2, 2, 2]
atoms["Ethanol"] = [3, 5, 1, 2, 2, 4, 2, 2, 2]
atoms["Acetone"] = [4, 2, 1, 1, 3, 3, 3, 3, 3, 3]
atoms["SPC_T320"] = [2, 1, 1]
atoms["SPC_T340"] = [2, 1, 1]
atoms["SPC_T360"] = [2, 1, 1]

Na = {k: len(v) for (k, v) in atoms.items()}  # num atoms in molecule


args = docopt(__doc__)
sol1 = args["<sol1>"]
sol2 = args["<sol2>"]
Nmol1 = int(args["<nmol1>"])
Nmol2 = int(args["<nmol2>"])
if sol1 not in masses.keys():
    sys.exit("Solvent %s not present." % sol1)
if sol2 not in masses.keys():
    sys.exit("Solvent %s not present." % sol2)

fname = args["<frame>"]
A = Atoms().read(fname)
A.names = list(map(int, A.names))
exp_Nat = Nmol1 * Na[sol1] + Nmol2 * Na[sol2] # expected number of atoms
if len(A.names) != exp_Nat:
    sys.exit("Number of atoms in xyz file %i does not correspond to input %i."\
            % (len(A.names), exp_Nat))

names = np.array([1] * Nmol1 + [2] * Nmol2)
at_masses1 = [masses[sol1][i] for i in atoms[sol1]]
at_masses2 = [masses[sol2][i] for i in atoms[sol2]]
mol_mass1 = sum(at_masses1)
mol_mass2 = sum(at_masses2)
print("Molecular masses: %.3f, %.3f" % (mol_mass1, mol_mass2))
print("Atomic masses: %s, %s" % (at_masses1, at_masses2))

xyz_com = np.zeros((Nmol1 + Nmol2, 3))
for i in range(Nmol1):
    for j in range(Na[sol1]):
        xyz_com[i] += A.xyz[Na[sol1]*i + j] * at_masses1[j] / mol_mass1
for i in range(Nmol2):
    for j in range(Na[sol2]):
        xyz_com[Nmol1+i] += A.xyz[Na[sol1]*Nmol1 + Na[sol2]*i + j] \
                * at_masses2[j] / mol_mass2

B = Atoms(names, xyz_com)
outname = fname.rstrip(".xyz") + "_com.xyz"
B.save(outname, vmd=True)


