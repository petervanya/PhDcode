#!/usr/bin/env python
"""Usage:
    gen_Pt_cluster.py water [-s <s>] [-p <p>] [-r <r>] [--save <fname>]
    gen_Pt_cluster.py Pt <cluster> [-s <s>] [--save <fname>]
    gen_Pt_cluster.py carbon <Nh> <Nv> [--addH [--Ha <Ha>]] [-s <s>] [--save <fname>]
    gen_Pt_cluster.py --atom <a> -s <s>

Generate xyz files of some molecules or clusters: 
Pt custer, water molecule, carbon sheet

Arguments:
    <cluster>              Pt cluster type, e.g. 9_10_9
    <Nh>                   Number of carbon rings horizontally
    <Nv>                   Number of carbon rings vertically

Options:
    -s <s>,--shift <s>     Shift 0th atom by "x y z" [default: 0 0 0]
    -p <p>,--posPt <p>     Position of 0th atom in Pt lattice vectors "x y z" [default: 0 0 0]
    -r <r>,--rotate <r>    Rotate atoms by theta, phi [default: 0 0]
    --save <fname>         Save xyz coords
    --addH                 Add hydrogens around carbon sheet perimeter
    --Ha <Ha>              Hydrogen distance from first carbon atom [default: 1.0]

pv278@cam.ac.uk, 17/06/15
"""
from docopt import docopt
import numpy as np
from numpy.linalg import norm
from math import radians, sqrt
from xyzlib import Atoms


def Pt_basis(a=2.775):
    """Basis vectors for Pt lattice"""
    v = a*sqrt(3.0)/2
    h = a*sqrt(2.0/3)
    basis = np.array([[a,   a/2, a/2],
                      [0.0, v,   v/3],
                      [0.0, 0.0, h  ]])
    return basis


def gen_water(shift=np.zeros(3), shiftPt=np.zeros(3), theta=0.0, phi=0.0):
    """Init water molecule, already optimised using 6-311G**"""
    coords = np.zeros((3,3))
    coords[1, 1] += 0.757009
    coords[2, 1] += -0.757009
    coords[1:3, 2] += 0.593565
    names = ["O", "H", "H"]
    mol = Atoms(names, coords)
    mol.rotate(theta, phi)
    mol.shift(np.dot(Pt_basis(), shiftPt))
    mol.shift(shift)
    return mol


def gen_Pt(a=2.775, cluster="3", shift=np.zeros(3), shiftPt=np.zeros(3),\
           theta=0.0, phi=0.0):
    """
    Generate small Pt clusters. Available clusters:
    * one-layer: 3, 4, 68, 12
    * two-layer: 6_3, 8_4, 12_7
    * three-layer: 5_10_5, 9_10_9
    """
    v = sqrt(3.0)/2*a
    h = sqrt(2.0/3.0)*a

    def get_cluster(cluster):
        """get the cluster type from string input"""
        arr = [int(i) for i in cluster.split("_")]
        Natoms = sum(arr)
        return arr, Natoms
        
    def get_coords3():
        coords = np.zeros((3,3))
        coords[1,0] = a
        coords[2,0] = a/2
        coords[2,1] = v
        coords += shift
        return coords
    
    def get_coords4():
        coords = np.zeros((4,3))
        for i in range(2):
          coords[1+i,0] = -a/2 + i*a
          coords[1+i,1] = v
        coords[3,1] = 2*v
        coords += shift
        return coords
        
    def get_coords6():
        coords = np.zeros((6,3))
        for i in range(3):
          coords[i,0] = i*a
        for i in range(2):
          coords[i+3,0] = a/2 + i*a
          coords[i+3,1] = v
        coords[5,0] = a
        coords[5,1] = 2*v
        coords += shift
        return coords
        
    def get_coords7():
        coords = np.zeros((7,3))
        for i in range(2):
          coords[i,0] = i*a
        for i in range(3):
          coords[2+i,0] = -a/2 + i*a
          coords[2+i,1] = v
        for i in range(2):
          coords[5+i,0] = i*a
          coords[5+i,1] = 2*v
        coords += shift
        return coords
        
    def get_coords8():
        coords = np.zeros((8,3))
        coords[1,0] = a
        for i in range(3):
          coords[2+i,0] = -a/2 + i*a
          coords[2+i,1] = v
        for i in range(2):
          coords[5+i,0] = i*a
          coords[5+i,1] = 2*v
        coords[7,0] = a/2
        coords[7,1] = 3*v
    
        coords += shift
        return coords
        
    def get_coords12():
        coords = np.zeros((12,3))
        for i in range(3):
          coords[i,0] = i*a
        for i in range(4):
          coords[3+i,0] = -a/2 + i*a
          coords[3+i,1] = v
        for i in range(3):
          coords[7+i,0] = i*a
          coords[7+i,1] = 2*v
        for i in range(2):
          coords[10+i,0] = a/2 + i*a
          coords[10+i,1] = 3*v
        coords += shift
        return coords
        
    def get_coords10():
        coords = np.zeros((10,3))
        for i in range(3):
          coords[i,0] = a*i
        for i in range(4):
          coords[i+3,0] = -a/2 + a*i
          coords[i+3,1] = v
        for i in range(3):
          coords[i+7,0] = a*i
          coords[i+7,1] = 2*v
        coords += shift
        return coords
    
    def get_coords6_3():
        arr,Natoms = get_cluster(cluster)
        coords = np.zeros((Natoms,3))
        coords[:6,:] = get_coords6(a,[0,0,0])
        coords[6:,:] = get_coords3(a,[a/2,v/3,h])
        coords += repmat(shift,Natoms,1)
        return coords
    
    def get_coords8_4():
        arr,Natoms = get_cluster(cluster)
        coords = np.zeros((Natoms,3))
        coords[:8,:] = get_coords8(a)
        coords[8:,:] = get_coords4(a,[a/2,v/3,h])
        coords += repmat(shift,Natoms,1)
        return coords
        
    def get_coords12_7():
        arr,Natoms = get_cluster(cluster)
        coords = np.zeros((Natoms,3))
        coords[:arr[0],:] = get_coords12(a)
        coords[arr[0]:,:] = get_coords7(a,[a/2,v/3,h])
        coords += repmat(shift,Natoms,1)
        return coords
    
    def get_coords5_10_5():
        arr,Natoms = get_cluster(cluster)
        L1,L2,L3 = arr
        coords = np.zeros((Natoms,3))
        # 1st layer
        for i in range(3):
            coords[i,0] = a*i
            coords[i,1] = 2*v/3
        for i in range(2):
            coords[i+3,0] = a/2 + a*i
            coords[i+3,1] = 5*v/3
        for i in range(5):
            coords[i,2] = -h
        # 2nd layer
        coords[L1:L1+L2,:] = get_coords10(a)
        # 3rd layer
        for i in range(2):
            coords[L1+L2+i,0] = a/2 + a*i
            coords[L1+L2+i,1] = v/3
        for i in range(3):
            coords[L1+L2+i+2,0] = a*i
            coords[L1+L2+i+2,1] = 4*v/3
        for i in range(5):
            coords[L1+L2+i,2] = h
        coords += repmat(shift,Natoms,1)
        return coords
    
    def get_coords9_10_9():
        arr,Natoms = get_cluster(cluster)
        L1,L2,L3 = arr
        coords = np.zeros((Natoms,3))
        # 1st layer
        for i in range(4):
            coords[i,0] = -a/2 + a*i
            coords[i,1] = -v/3
        for i in range(3):
            coords[i+4,0] = a*i
            coords[i+4,1] = -v/3 + v
        for i in range(2):
            coords[i+7,0] = a/2 + a*i
            coords[i+7,1] = -v/3 + 2*v
        for i in range(L1):
            coords[i,2] = -h
        # 2nd layer
        coords[L1:L1+L2,:] = get_coords10(a)
        # 3rd layer
        for i in range(2):
            coords[L1+L2+i,0] = a/2 + a*i
            coords[L1+L2+i,1] = v/3
        for i in range(3):
            coords[L1+L2+i+2,0] = a*i
            coords[L1+L2+i+2,1] = v/3 + v
        for i in range(4):
            coords[L1+L2+i+5,0] = -a/2 + a*i
            coords[L1+L2+i+5,1] = v/3 + 2*v
        for i in range(L3):
            coords[L1+L2+i,2] = h
        coords += repmat(shift,Natoms,1)
        return coords

    # one layer
    if cluster == "3":    mol = Atoms(["Pt"]*3, get_coords3())
    elif cluster == "4":  mol = Atoms(["Pt"]*4, get_coords4())
    elif cluster == "6":  mol = Atoms(["Pt"]*6, get_coords6())
    elif cluster == "7":  mol = Atoms(["Pt"]*7, get_coords7())
    elif cluster == "8":  mol = Atoms(["Pt"]*8, get_coords8())
    elif cluster == "10": mol = Atoms(["Pt"]*10, get_coords10())
    elif cluster == "12": mol = Atoms(["Pt"]*12, get_coords12())
    # two layers
    elif cluster == "6_3":  mol = Atoms(["Pt"]*9, get_coords6_3())
    elif cluster == "8_4":  mol = Atoms(["Pt"]*12, get_coords8_4())
    elif cluster == "12_7": mol = Atoms(["Pt"]*19, get_coords12_7())
    # three layers
    elif cluster == "5_10_5": mol = Atoms(["Pt"]*20, get_coords5_10_5())
    elif cluster == "9_10_9": mol = Atoms(["Pt"]*28, get_coords9_10_9())
    else:
        print "Cluster not implemented, please choose a different one."
        raise SystemExit
    return mol

# ===== carbon
def create_one_ring(a, shift):
    """Create one carbon ring"""
    coords = np.array([[0,      0,           0],
                       [-a/2,   a*sqrt(3)/2, 0],
                       [0.0,    a*sqrt(3),   0],
                       [a,      a*sqrt(3),   0],
                       [a*3/2,  a*sqrt(3)/2, 0],
                       [a,      0,           0]])
    coords += shift
    return coords


def gen_carbon(Nh, Nv, a=1.42, shift=np.zeros(3)):
    """Generate carbon sheet given number of rings"""
    N = Nh*Nv

    def remove_overlapped(coords):
        """Remove overlapped carbon atoms"""
        eps = 1.0e-5
        i = 0
        j = 1
        while(i < len(coords)):
            while(j < len(coords)):
                if norm(coords[i] - coords[j]) < eps:
                    coords = np.delete(coords, j, axis=0)
                j += 1
            i += 1
            j = i + 1
        return coords

    dx = 3*a/2
    dy = a*sqrt(3.0)
    coords = create_one_ring(a, shift)
    for n in range(1, N):
      	j = n%Nv
      	i = n/Nv #(n-j)/Nv
      	base = [i*dx, j*dy + i*dy/2, 0.0]
      	coords = np.vstack((coords, create_one_ring(a, base)))
    coords = remove_overlapped(coords)
    return Atoms(["C"]*len(coords), coords)


def add_hydrogens(rings, a, b):
    """Add hydrogen atoms around carbon sheet
    to make coronene-like molecule"""
    N = len(rings)
    H_list = []
    tol = a + 0.1                  # for each atom, check atoms in this radius
    for i in range(N):
        bonds = []                 # bonds for ith carbon atom
        for j in range(N):
            if norm(rings[j] - rings[i]) < tol and i != j:
                bonds.append(j)
        if len(bonds) == 2:        # edge atoms have two bonds
            v = -(rings[bonds[0]] - rings[i] + rings[bonds[1]] - rings[i])
            H_rings = list(rings[i] + v/norm(v)*b)
            H_list.append(H_rings)
    return Atoms(["H"]*len(H_list), H_list)


if __name__ == "__main__":
    args = docopt(__doc__, version=1.0)
#    print args
    theta, phi = (radians(float(i)) for i in args["--rotate"].split())
    shift = np.array(args["--shift"].split()).astype(float)
    posPt = np.array(args["--posPt"].split()).astype(float)

    if args["water"]:
        mol = gen_water(shift=shift, shiftPt=posPt, theta=theta, phi=phi)
        print mol
        if args["--save"]:
            mol.save(args["--save"])    # "water.xyz"
    elif args["Pt"]:
        cluster = args["<cluster>"]
        mol = gen_Pt(cluster=cluster, shift=shift)
        print mol
        if args["--save"]:
            mol.save(args["--save"])    # "Pt.xyz"
    elif args["carbon"]:
        a = 1.42
        Nh, Nv = int(args["<Nh>"]), int(args["<Nv>"])
        mol = gen_carbon(Nh, Nv, shift=shift)
        if args["--addH"]: 
            b = float(args["--Ha"])
            mol += add_hydrogens(mol.coords, a, b)
        print mol
        if args["--save"]:
            mol.save(args["--save"])
    else:                   # single atom
        name = args["--atom"]
        coords = [0, 0, 0]
        xyz = Atoms(name, coords)
        print xyz
        if args["--save"]:
            xyz.save(args["--save"])


