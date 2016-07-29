#!/usr/bin/env python
"""
Analyse the surface tensions of various
spin configurations of Pt clusters w or w/o water

pv278@cam.ac.uk, 07/04/15
"""
from docopt import docopt
import numpy as np
from math import sqrt

def isfloat(value):
    """check if character is float"""
    try:
        float(value)
        return True
    except ValueError:
        return False

def print_table(res):
    header = "Spin \t Energy (eV) \t Surf. tension (J/m**2)"
    print header
    print "-"*len(header)
    N = len(res)
    for i in range(N):
        print res[i][0],'\t {0:.5f} \t {1:.5f}'.format(res[i][1],res[i][2])
    
def parse_data(file):
    data = []
    f = open(file,"r")
    for line in f:
        data.append( line.rstrip("\n").split(" \t ") )
    f.close()
    return data

if __name__ == "__main__":
    #Ewater = -76.0107464919
    Ewater = -76.393551           # ZP corrected B3LYP/LANL2DZ
    a = 2.775
    S = a**2*sqrt(3.0)/2.0*1e-20  # in m**2
    
    data_p = parse_data("Pt4_plain.out")
    data_w = parse_data("Pt4_water.out")
    res = []                # output table
    N = len(data_p)

    for i in range(N):
        Ew = data_w[i][2]   # Pt with water
        Ep = data_p[i][2]   # plain Pt
        if isfloat(Ew) and isfloat(Ep):
            dE = float(Ew) - Ewater - float(Ep)
            res.append([int(data_w[i][0]), dE*27.21138505])
    
    for i in range(len(res)):
        res[i].append(res[i][1]*1.602e-19/S)
    
    print_table(res)
    
    
