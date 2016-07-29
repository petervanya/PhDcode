#!/usr/bin/env python
"""Usage: 
    gen_Pt.py <cluster> [-s <s>] [-d <d>]

Generate coordinates of atoms for certain Pt clusters
according to Jacob et al.

Arguments:
    <cluster>              Pt cluster type, e.g. 9_10_9

Options:
    -h --help              Show this message and exit
    -s <s>,--shift <s>     Shift the initial atom by a vector [x,y,z]
    -d <d>,--dir <d>       Directory to save the produced xyz file

pv278@cam.ac.uk, 12/2014
"""
from docopt import docopt
import numpy as np
from numpy.matlib import repmat
from math import sqrt
from iolib import save_xyz


cluster_list = open("cluster_list.out").readlines()


def get_cluster(cluster):
    """get the cluster type from string input"""
    arr = [int(i) for i in cluster.split("_")]
    Natoms = sum(arr)
    return arr, Natoms


def get_coords3(a,shift=[0,0,0]):
    v = sqrt(3.0)/2*a
    coords = np.zeros((3,3))
    coords[1,0] = a
    coords[2,0] = a/2
    coords[2,1] = v
    
    coords += shift
    return coords


def get_coords4(a,shift=[0,0,0]):
    v = sqrt(3.0)/2*a
    coords = np.zeros((4,3))
    for i in range(2):
      coords[1+i,0] = -a/2 + i*a
      coords[1+i,1] = v
    coords[3,1] = 2*v
    
    coords += shift
    return coords


def get_coords6(a,shift=[0,0,0]):
    v = sqrt(3.0)/2*a
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
    

def get_coords7(a,shift=[0,0,0]):
    v = sqrt(3.0)/2*a
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
    

def get_coords8(a,shift=[0,0,0]):
    v = sqrt(3.0)/2*a
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
    

def get_coords12(a,shift=[0,0,0]):
    v = sqrt(3.0)/2*a
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
    

def get_coords10(a,shift=[0,0,0]):
    v = sqrt(3.0)/2*a
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


def get_coords6_3(a,shift=[0,0,0]):
    v = sqrt(3.0)/2*a
    arr,Natoms = get_cluster(cluster)
    coords = np.zeros((Natoms,3))
    coords[:6,:] = get_coords6(a,[0,0,0])
    coords[6:,:] = get_coords3(a,[a/2,v/3,h])
    coords += repmat(shift,Natoms,1)
    return coords


def get_coords8_4(a,shift=[0,0,0]):
    v = sqrt(3.0)/2*a
    arr,Natoms = get_cluster(cluster)
    coords = np.zeros((Natoms,3))
    coords[:8,:] = get_coords8(a)
    coords[8:,:] = get_coords4(a,[a/2,v/3,h])
    coords += repmat(shift,Natoms,1)
    return coords


def get_coords12_7(a,shift=[0,0,0]):
    v = sqrt(3.0)/2*a
    arr,Natoms = get_cluster(cluster)
    coords = np.zeros((Natoms,3))
    coords[:arr[0],:] = get_coords12(a)
    coords[arr[0]:,:] = get_coords7(a,[a/2,v/3,h])
    coords += repmat(shift,Natoms,1)
    return coords


def get_coords5_10_5(a,shift=[0,0,0]):
    v = sqrt(3.0)/2*a
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


def get_coords9_10_9(a,shift=[0,0,0]):
    v = sqrt(3.0)/2*a
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


if __name__  ==  "__main__":
    a = 2.775                    # atom distance in Angstroms
    v = sqrt(3.0)/2*a
    h = sqrt(2.0/3.0)*a

    args = docopt(__doc__,version=2.0)
    
    cluster = args["<cluster>"]
    if args["--shift"]:
        shift = np.array(args["--shift"].split()).astype(float)
        print("Shift =", shift)
    else:
        shift = np.zeros(3)
    
    if args["--dir"]:
        dirr = args["--dir"] + "/"
        print(dirr)
    else:
        dirr = ""
    filename = str(dirr + "Pt.xyz")
    
    # one layer
    if cluster == "3":    coords = get_coords3(a,shift)
    
    elif cluster == "4":  coords = get_coords4(a,shift)
    
    elif cluster == "6":  coords = get_coords6(a,shift)
    
    elif cluster == "7":  coords = get_coords7(a,shift)
    
    elif cluster == "8":  coords = get_coords8(a,shift)
    
    elif cluster == "10": coords = get_coords10(a,shift)
       
    elif cluster == "12": coords = get_coords12(a,shift)
    # two layers
    elif cluster == "6_3":  coords = get_coords6_3(a,shift)
       
    elif cluster == "8_4":  coords = get_coords8_4(a,shift)
       
    elif cluster == "12_7": coords = get_coords12_7(a,shift)
    # three layers
    elif cluster == "5_10_5": coords = get_coords5_10_5(a,shift)
         
    elif cluster == "9_10_9": coords = get_coords9_10_9(a,shift)
       
    else:
        print("Cluster not implemented, please choose from the following:")
        print(", ".join([s.rstrip("\n") for s in cluster_list]))
        raise SystemExit
    
    # save to file
    names = ["Pt"] * len(coords)
    save_xyz(coords, names, filename)
   

