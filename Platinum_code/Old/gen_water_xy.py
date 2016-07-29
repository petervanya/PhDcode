#!/usr/bin/python
# =====
# script to generate water coordinates shifted by h in xy plane
# 20/02/15
# =====
import sys
import numpy as np
import numpy.matlib
from math import pi, sqrt, sin, cos, radians

def printH2O(coordsH2O,filename):
		f=open(filename,"w")
		s = str("O\t")
		for j in range(3):
			s += "%.3f" % coordsH2O[0,j] + "\t"
		f.write(s+"\n")
		for i in range(1,3):
			s = str("H\t")
			for j in range(3):
				s += "%.3f" % coordsH2O[i,j] + "\t"
			f.write(s+"\n")
		f.close()

# ===== get input
a=2.775       # platinum constant
dist=np.zeros(3)

if sys.argv[1] == "top":
  dist[0] = 0.0
  dist[2] = 0.0
elif sys.argv[1] == "bridge":
  dist[0] = 0.0
  dist[2] = a/2
elif sys.argv[1] == "fcc":
  dist[0] = a*cos(pi/6)/3
  dist[2] = a/2
else:
  raise NotImplementedError

if len(sys.argv)>2:
  dist[1] = sys.argv[2]                   # shift in y-dir, out of Pt atoms plane
else:
  print "No vertical distance specified, use 1.0 A."
  dist[2] = 1.0
print "Shift distance: ",dist

# ===== water coords

theta = radians(104.45)                   # angle between H atoms
l_OH = 0.9584                             # bond length between O and H
coords=np.zeros((3,3))
coords[1,0] += l_OH*sin(theta/2)          # shift x-coord of 1st H
coords[2,0] += -l_OH*sin(theta/2)         # shift x-coord of 2nd H
coords[1:3,1] += l_OH*cos(theta/2)        # shift z-coords of Hs

Dist = np.matlib.repmat(dist,3,1)         # shift all atoms by dist
coords += Dist

# ===== print into file
filename="water.xyz"
printH2O(coords,filename)
print "Water coords printed into file water.xyz"


