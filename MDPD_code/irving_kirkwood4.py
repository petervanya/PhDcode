#!/usr/bin/env python
"""Usage:
    irving_kirkwood.py <frames> [--nu <nu>]
    
Compute surface tension in x-direction on a grid of <nu> points.

Options:
    --nu <nu>    Number of points in u-dir [default: 100]

Last modified: 05/06/17
"""
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt
from dlms_lib import read_xyzfile2
import os, glob, time
from docopt import docopt


def parse_field():
    """Extract coefficients from the field file.
    Assume rc = 1.0"""
    try:
        f = open("FIELD", "r").readlines()
    except FileNotFoundError:
        sys.exit("FIELD file not found in the current dir.")
    for line in range(len(f)):
        if "SPECIES" in f[line].upper():
            Nb = int(f[line].split()[-1])
            sp_lines = f[line+1 : line+Nb+1]
        if "INTERACTIONS" in f[line].upper():
            Ni = int(f[line].split()[-1])
            int_lines = f[line+1 : line+Ni+1]

    sp = {}       # species
    sp_num = {}   # species dict with numbers
    for i in range(Nb):
        tmp = sp_lines[i].split()
        sp[tmp[0]] = i
        sp_num[tmp[0]] = int(tmp[3])

    int_d = {}    # interaction dict
    for i in range(Ni):
        tmp = int_lines[i].split()
        if tmp[2] != "mdpd":
            sys.exit("Interaction type: %s. Should be: mdpd." % tmp[1])
        int_d[tmp[0] + " " + tmp[1]] = tmp[3:5]

    int_A = np.zeros((Nb, Nb))  # interaction matrix
    int_B = np.zeros((Nb, Nb))
    for k, v in int_d.items():
        i, j = [sp[l] for l in k.split()]
        int_A[i, j] = v[0]
        int_A[j, i] = v[0]
        int_B[i, j] = v[1]
        int_B[j, i] = v[1]
    return int_A, int_B


def parse_control():
    """Return volume"""
    rd = 0.0
    try:
        f = open("CONTROL", "r").readlines()
    except FileNotFoundError:
        sys.exit("CONTROL file not found in the current dir.")
    for line in f:
        if "manybody cutoff" in line.lower():
            rd = float(line.split()[-1])
        if "volume" in line.lower():
            tmp = line.split()
            if len(tmp) == 2:
                L = float(tmp[-1])**(1/3)
                L = np.array([L, L, L])
            elif len(tmp) == 4:
                L = np.array(list(map(float, tmp[1:4])))
    return L, rd


#def boxsplit_in_z(xyz, z):
#    """2 for z, 0 for x"""
#    xyz1 = xyz[(xyz[:, 0] < z)]
#    xyz2 = xyz[(xyz[:, 0] > z)]
#    return xyz1, xyz2

def boxsplit_in_z(xyz, p):
    ## p is the distance of the plane
    # need rc ==1 , if rc changes scale "local"
    direction = 0 #2 for z, 0 for x
    L_direction = L[direction]
    local = 2
    far_side = L_direction-local
    if (p > local) and (p < (far_side)):
        xyz1 = xyz[np.logical_and((xyz[:,direction]>p-local),(xyz[:, direction]<p))]
        xyz2 = xyz[np.logical_and((xyz[:,direction]<p+local),(xyz[:, direction]>p))]
    elif p < local:
        xyz1 =  xyz[np.logical_or((xyz[:,direction]>far_side),(xyz[:, direction]<p))]
        xyz2 =  xyz[np.logical_and((xyz[:,direction]<p+local),(xyz[:, direction]>p))]
    elif p > far_side:
        xyz1 =  xyz[np.logical_and((xyz[:,direction]>p-local),(xyz[:, direction]<p))]
        xyz2 =  xyz[np.logical_or((xyz[:,direction]<local),(xyz[:, direction]>p))]

    return xyz1, xyz2


def perodic_1D_difference(xyz1, xyz2, Li):
    outer_product = np.subtract.outer(xyz1, xyz2)
    ## periodic
    outer_product = np.where((np.abs(outer_product) > Li/2), \
            (Li - np.abs(outer_product)) * -1 * \
            np.sign(outer_product), outer_product)
    return outer_product


def difference_matrix(xyz1, xyz2):
    """need to get in the periodic boundary"""
    dx = perodic_1D_difference(xyz1[:, 0], xyz2[:, 0], Lx)
    dy = perodic_1D_difference(xyz1[:, 1], xyz2[:, 1], Ly)
    dz = perodic_1D_difference(xyz1[:, 2], xyz2[:, 2], Lz)
    dr = np.sqrt(dx**2 + dy**2 + dz**2)
    return dr, dx, dy, dz


def local_density(r):
    weight = 15 / (2 * np.pi * rd**3) 
    mask = r < rd         # eliminate too far 
    mask_z = r > 0        # eliminates self interference
    rho = weight * mask * mask_z * (1 - np.divide(np.abs(r), rd))**2
    rho = np.sum(rho, 1)
    return rho


args = docopt(__doc__)
files = sorted(glob.glob(args["<frames>"]))
Nf = len(files)
Nu = int(args["--nu"])

A, B = parse_field()
A, B = A[0, 0], B[0, 0]
L, rd = parse_control()
Lx, Ly, Lz = L
print("Surface tension calculation along x-coord")
print("rd: %.2f | L: %s" % (rd, L))
print("A: %.2f | B: %.2f" % (A, B))

pT_avg = np.zeros(Nu)
du = Lx / Nu
u = np.arange(0, Lx, du)
gammas = np.zeros(Nf)
print("Nu: %i | du: %.3f" % (Nu, du))
   
ti = time.time()
for j in range(Nf):
    _, xyz = read_xyzfile2(files[j])
    dm = difference_matrix(xyz, xyz)
    local_den = local_density(dm[0])
    xyz = np.c_[xyz, local_den]

    pT = np.zeros(Nu)
    P = np.zeros((Nu, 3))
   
    for i in range(Nu):
        left, right = boxsplit_in_z(xyz, u[i])
        r, x, y, z = difference_matrix(left, right)
        
        mask_z = (r > 0)
        mask_a = (r < 1)
        mask_b = (r < rd)
        
        F_A = A * (1 - np.abs(r)) * mask_a * mask_z
        F_B = B * (1 - np.abs(r) / rd) * np.add.outer(left[:, 3], right[:, 3]) * \
                mask_b * mask_z
        F = F_A + F_B

        F_x = np.multiply(F, np.nan_to_num(np.divide(x, np.abs(r))))
        F_y = np.multiply(F, np.nan_to_num(np.divide(y, np.abs(r))))
        F_z = np.multiply(F, np.nan_to_num(np.divide(z, np.abs(r))))
	    
        p_xx = np.sum(np.multiply(x / np.abs(x), F_x))
        p_yy = np.sum(np.multiply(y / np.abs(x), F_y))
        p_zz = np.sum(np.multiply(z / np.abs(x), F_z))
        P[i] = [p_xx, p_yy, p_zz]
        pT[i] = (p_xx - (p_yy + p_zz) / 2) / (Ly * Lz)

    gammas[j] = simps(pT, u) / 2.0
    print(files[j], "| Gamma =", gammas[j])
    pT_avg += pT / Nf
   
tf = time.time()
print("Time: %.2f s." % (tf - ti))

np.savetxt("gamma.out", np.c_[u, pT_avg])
gamma_avg = simps(pT_avg, u) / 2.0
print("Final gamma: %f | Std: %f" % (gamma_avg, np.std(gammas)))


