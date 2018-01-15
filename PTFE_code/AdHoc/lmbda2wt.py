#!/usr/bin/env python
"""
Convert water uptake lmbda to water content as a percentage.

Formula for Wu   : wt = 1 / (1 + 6 Nbm / (lmbda - 3))
Formula for Yamam: wt = 1 / (1 + 4 Nbm / lmbda)

03/05/17
"""


def lmbda2wt(l, N, Nmc, Nbm):
    Nc = N / (Nmc * ((l - 3) / 6 + Nbm))
    Nwb = (l - 3) / 6 * Nmc * Nc
    return Nwb / N * (1 + 3 / (l - 3))
    

Nbm = 5
Nmc = 15
rho = 3
L = 40.0
N = rho * L**3
lmbdas = list(range(4, 25, 2))
wts = [lmbda2wt(l, N, Nmc, Nbm) for l in lmbdas]

#import matplotlib.pyplot as plt
#plt.plot(lmbdas, wts, "*-")
#plt.show()

print("Wu parametrisation:")
print("Beads per monomer: %i | Monomers per chain: %i" % (Nbm, Nmc))
for i in range(len(lmbdas)):
    print("%2.i   %.3f" % (lmbdas[i], wts[i]))

print("Yamamoto parametrisation:")
Nbm = 6
lmbdas = list(range(2, 25, 2))
wts_y = [l / (l + 4 * Nbm) for l in lmbdas]
print("Beads per monomer: %i | Monomers per chain: %i" % (Nbm, Nmc))
for i in range(len(lmbdas)):
    print("%2.i   %.3f" % (lmbdas[i], wts_y[i]))
