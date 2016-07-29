#!/usr/bin/python
# =====
# 20/02/15
# Generate coords of water in xz plane
# 3 cmd ln args -- x,y,z shift from 0,0,0
# =====
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
dist=np.zeros(3)

if len(sys.argv)<4:
	print "You did not specify enough shift components, use default vector: [0,0,1] A."
	dist = [0,0,1]
else:
	for i in range(3):
		dist[i] = float(sys.argv[i+1])
	print "Shift distance: ",dist

# ===== water coords

theta = radians(104.45)                   # angle between H atoms
l_OH = 0.9584                             # bond length between O and H
coords=np.zeros((3,3))
coords[1,1] += l_OH*sin(theta/2)          # shift x-coord of 1st H
coords[2,1] += -l_OH*sin(theta/2)         # shift x-coord of 2nd H
coords[1:3,2] += l_OH*cos(theta/2)        # shift z-coords of Hs

coords += np.matlib.repmat(dist,3,1)      # shift all atoms by dist

# ===== print into file
filename="water.xyz"
printH2O(coords,filename)
print "Water coords printed into file water.xyz"


