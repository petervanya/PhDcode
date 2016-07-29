#!/usr/bin/env python
"""Usage:
    default_input.py 


Source for Hildebrand solubility parameters:
* http://cool.conservation-us.org/byauth/burke/solpar/solpar2.html
* https://en.wikipedia.org/wiki/Hildebrand_solubility_parameter
Compressibility:
* Taravillo et al. http://pubs.acs.org/doi/pdf/10.1021/je060415l
* wiki
* Bulk modulus of water: 2.15e9 | methanol: 8.23e8
Methanol expt:
* Hruby et al. doi:10.1006/jcht.1993.1121
* Kubota et al. doi:10.1007/BF00503224 READ THIS!

pv278@cam.ac.uk, 21/07/16
"""
from docopt import docopt
from collections import namedtuple
from math import sqrt

NA = 6.022e23
AMU = 1.66e-27
R = 8.314
kB = 1.381e-23

# rho (kg/m**3), molar mass (kg/mol), sol (MPa**0.5),
# enthalpy of vaporisation (J/mol), isothermal compressibility
subst = namedtuple("subst", ["rho", "Mm", "sol", "Hv", "kappaT"])
water = subst(1000.0, 18e-3, 47.8, 40.65e3, 4.65e-10)
meoh = subst(792.0, 32.04e-3, 29.7, 35.3e3, 1.22e-9)
Pt = subst(21450.0, 195e-3, 0.0, 510e3, 4.34e-12)
PE = subst(1178.0, 28e-3, 0.0, 0.0, 0.0) # CHECK


def mol_size(subst):
    return (subst.Mm*1e3*AMU/subst.rho)**(1/3)


def number_density(subst):
    return subst.rho * NA / subst.Mm


def Vm(subst):
    """In cm**3"""
    return subst.Mm / subst.rho * 1e6


def solubility(subst, T=300.0):
    """In MPa**0.5"""
    return sqrt((subst.Hv - R*T)/(subst.Mm/subst.rho)/1e6)


def get_chi(s1, s2, T=300.0):
    """Default T = 300 K"""
    return (Vm(s1) + Vm(s2)) / (2*R*T) * (s1.sol - s2.sol)**2


def inv_kappa(subst, T=300.0):
    """Dimensionless compressibility"""
    n = number_density(subst)
    return 1.0 / (n*kB*T*subst.kappaT)


def aii(subst):
    """Like bead interaction parameter
    in the units of kT/rho
    For water, a ~= 75 kT/rho"""
    k = inv_kappa(subst)
    alpha = 0.101            # GW bulgarian constant
    return (k - 1)/(2*alpha)


def travis_aii(subst, mols_per_bead, rho=3.0, T=300.0):
    """Eq. (30) from Travis_JCP_2007
    SI form: 
    return solubility(subst)**2/(alpha * subst.rho**2 * rc**4)
    """
    alpha = 0.101
    rc = (mol_size(subst)**3 * mols_per_bead * rho)**(1./3)
    E0 = kB*T
    s2 = solubility(subst)**2 * 1e6
    s2_dimless = s2 * rc**3 / E0
    return s2_dimless**2 / (alpha * rho**2)


def match_mol_volumes(s1, s2, N=10, tol=0.1):
    """
    * N: upper number of molecules in a bead
    * tol: tolerance in percentage
    """
    if tol < 0. or tol > 1.0: raise ValueError("Enter tol between 0 and 1")
    print("\nMatching volumes:")
    for i in range(1, N+1):
        for j in range(1, N+1):
            V1 = i*mol_size(s1)**3
            V2 = j*mol_size(s2)**3
            r = abs(V1 - V2)/V1
            if r < tol:
                print("%i, %i (%.2e, %.2e) %.2f" %(i, j, V1, V2, r))


if __name__ == "__main__":
    T = 300.0
    print("Molar volumes | water: %.2f | meoh: %.2f | Pt: %.2f" % \
         (Vm(water), Vm(meoh), Vm(Pt)))
    print("Molecular sizes | water: %.2e | meoh: %.2e | Pt: %.2e" % \
         (mol_size(water), mol_size(meoh), mol_size(Pt)))
    print("Solubility from Hidebrand | water: %.2f | meoh: %.2f | Pt: %.2f" % \
          (solubility(water), solubility(meoh), solubility(Pt)))
    print("Dimensionless inverse kappa | water: %.2f | meoh: %.2f | Pt: %.2f" %\
         (inv_kappa(water), inv_kappa(meoh), inv_kappa(Pt)))

    print("=====")
    print("Like bead inter params a | water: %.2f | meoh: %.2f | Pt: %.2f" %\
         (aii(water), aii(meoh), aii(Pt)))

    mpb = 3  # molecules per bead
    print("Travis like bead inter params a | water: %.2e | meoh: %.2e | Pt: %.2e" %\
         (travis_aii(water, mpb), travis_aii(meoh, mpb), travis_aii(Pt, mpb)))
    print("Chi param water/meoh: %.2f" % get_chi(water, meoh))
    print("Chi param water/Pt: %.2f" % get_chi(water, Pt))

    match_mol_volumes(water, meoh, tol=0.15)


