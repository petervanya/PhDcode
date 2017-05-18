#!/usr/bin/env python
"""
Convert water uptake lmbda to water content as a percentage.

03/05/17
"""
import matplotlib.pyplot as plt


def lmbda2wt(l, N, Nmc, Nbm):
    Nc = N / (Nmc * ((l - 3) / 6 + Nbm))
    Nwb = (l - 3) / 6 * Nmc * Nc
    return Nwb / N * (1 + 3 / (l - 3))
    

Nbm = 5
Nmc = 15
N = 192000
lmbdas = list(range(4, 25))
wts = [lmbda2wt(l, N, Nmc, Nbm) for l in lmbdas]

#plt.plot(lmbdas, wts, "*-")
#plt.show()

for i in range(len(lmbdas)):
    print("%i   %.1f" % (lmbdas[i], 100*wts[i]))
