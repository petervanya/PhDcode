#!/usr/bin/python
# Script to generate Gaussian input file for 
# Pt clusters 9.10.9 with H2O on top
# 29/01/14
import sys
import numpy as np
from math import sqrt, sin, cos, radians


# ===== Pt
def printPt(coordsPt,filename):
    f=open(filename,"a")
    N_Pt=coordsPt.shape[0]
    for i in range(N_Pt):
      s=""
      for j in range(3):
        s += "%.3f" % coordsPt[i,j] +"\t"
      f.write( str("Pt\t")+s+"\n" )
    f.write("\n")
    f.close()
        
def get_config(config):
    arr = [int(i) for i in config.split("_")]
    N_Pt = sum(arr)
    return arr, N_Pt
    
def get_coords3(a,x0,y0,z0):
		v=sqrt(3.0)/2*a
		coords=np.zeros((3,3))
		for i in range(2):
			coords[i,0] = i*a
		coords[2,0] = a/2
		coords[2,1] = v
		
		coords[:,0] += x0
		coords[:,1] += y0
		coords[:,2] += z0
		return coords

def get_coords4(a,x0,y0,z0):
		v=sqrt(3.0)/2*a
		coords=np.zeros((4,3))
		for i in range(2):
			coords[i,0] = i*a
		for i in range(2):
			coords[2+i,0] = a/2 + i*a
			coords[2+i,1] = v
		
		coords[:,0] += x0
		coords[:,1] += y0
		coords[:,2] += z0
		return coords
    
def get_coords6(a,x0,y0,z0):
		v=sqrt(3.0)/2*a
		coords=np.zeros((6,3))
		for i in range(3):
		  coords[i,0] = i*a
		for i in range(2):
		  coords[i+3,0] = a/2 + i*a
		  coords[i+3,1] = v
		coords[5,0] = a
		coords[5,1] = 2*v
		
		coords[:,0] += x0
		coords[:,1] += y0
		coords[:,2] += z0
		return coords
		
def get_coords7(a,x0,y0,z0):
		v=sqrt(3.0)/2*a
		coords=np.zeros((7,3))
		for i in range(2):
			coords[i,0] = a/2 + i*a
		for i in range(3):
			coords[2+i,0] = i*a
			coords[2+i,1] = v
		for i in range(2):
			coords[5+i,0] = a/2 + i*a
			coords[5+i,1] = 2*v
		
		coords[:,0] += x0
		coords[:,1] += y0
		coords[:,2] += z0
		return coords
		
def get_coords8(a,x0,y0,z0):
		v=sqrt(3.0)/2*a
		coords=np.zeros((8,3))
		for i in range(2):
			coords[i,0] = a/2 + i*a
		for i in range(3):
			coords[2+i,0] = i*a
			coords[2+i,1] = v
		for i in range(3):
			coords[5+i,0] = a/2 + i*a
			coords[5+i,1] = 2*v
		
		coords[:,0] += x0
		coords[:,1] += y0
		coords[:,2] += z0
		return coords
		
def get_coords12(a,x0,y0,z0):
		v=sqrt(3.0)/2*a
		coords=np.zeros((12,3))
		for i in range(3):
			coords[i,0] = a/2 + i*a
		for i in range(4):
			coords[3+i,0] = i*a
			coords[3+i,1] = v
		for i in range(3):
			coords[7+i,0] = a/2 + i*a
			coords[7+i,1] = 2*v
		for i in range(2):
			coords[10+i,0] = a + i*a
			coords[10+i,1] = 3*v
		
		coords[:,0] += x0
		coords[:,1] += y0
		coords[:,2] += z0
		return coords
		
def get_coords10(a,x0,y0,z0):
		v=sqrt(3.0)/2*a
		coords=np.zeros((10,3))
		for i in range(3):
			coords[i,0] = a/2 + i*a
		for i in range(4):
			coords[i+3,0] = a*i
			coords[i+3,1] = v
		for i in range(3):
			coords[i+7,0] = a/2 + a*i
			coords[i+7,1] = 2*v

		coords[:,0] += x0
		coords[:,1] += y0
		coords[:,2] += z0
		return coords
		
# ==== H20
def printH20(coordsH2O,filename):
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
		
def setH2Ocoords(Pt_atom_coords,dist):
		theta = radians(104.45)                   # angle between H atoms
		l_OH = 0.9584                             # bond length between O and H
		coords=np.zeros((3,3))
		for i in range(3):
			coords[i,:] = Pt_atom_coords + dist			# set O atom
		coords[1,0] += l_OH*sin(theta/2)          # shift x-coord of 1st H
		coords[2,0] += -l_OH*sin(theta/2)         # shift x-coord of 2nd H
		coords[1:3,2] += l_OH*cos(theta/2)        # shift z-coords of Hs
		return coords
    
# ===============================================
# main program
# ===============================================
a=2.775           # atom distance in A
v=sqrt(3.0)/2*a
h=sqrt(2.0/3.0)*a

filename=str("test.xyz")
config=sys.argv[1]
arr,N_Pt = get_config(config)
print "Total number of Pt atoms: ",N_Pt


# one layer
if config=="3":
	coords=np.zeros((3,3))
	coords[1,0] = a
	coords[2,0] = a/2
	coords[2,1] = v

elif config=="6":
	coords=np.zeros((6,3))
	for i in range(3):
	 coords[i,0] = i*a
	for i in range(2):
	 coords[i+3,0] = a/2 + i*a
	 coords[i+3,1] = v
	coords[5,0] = a
	coords[5,1] = 2*v

elif config=="8":
	coords=np.zeros((8,3))
	for i in range(2):
	  coords[i,0] = a/2 + i*a
	for i in range(3):
	  coords[2+i,0] = i*a
	  coords[2+i,1] = v
	for i in range(2):
	  coords[5+i,0] = a/2 + i*a
	  coords[5+i,1] = 2*v
	coords[7,0] = a
	coords[7,1] = 3*v
   
elif config=="12":
	coords=np.zeros((12,3))
	for i in range(3):
	  coords[i,0] = a/2 + i*a
	for i in range(4):
	  coords[3+i,0] = i*a
	  coords[3+i,1] = v
	for i in range(3):
	  coords[7+i,0] = a/2 + i*a
	  coords[7+i,1] = 2*v
	for i in range(2):
	  coords[10+i,0] = a + i*a
	  coords[10+i,1] = 3*v

# two layers
elif config=="6_3":
	N_Pt=int(config[0])+int(config[2])
	print "Total number of atoms: ",N_Pt
	coords=np.zeros((N_Pt,3))
	coords[:6,:] = get_coords6(a,0,0,0)
	coords[6:,:] = get_coords3(a,a/2,v/3,h)
   
elif config=="8_4":
   N_Pt=int(config[0])+int(config[2])
   print "Total number of atoms: ",N_Pt
   coords=np.zeros((N_Pt,3))
   coords[:8,:] = get_coords8(a,0,0,0)
   coords[8:,:] = get_coords4(a,a/2,2*v/3,h)
   
elif config=="12_7":
   arr,N_Pt = get_config(config)
   print "Total number of atoms: ",N_Pt
   coords=np.zeros((N_Pt,3))
   coords[:arr[0],:] = get_coords12(a,0,0,0)
   coords[arr[0]:,:] = get_coords7(a,a/2,v/3,h)

# three layers
elif config=="5_10_5":
   arr,N_Pt = get_config(config)
   print "Total number of atoms: ",N_Pt
   L1=arr[0]
   L2=arr[1]
   L3=arr[2]
   coords=np.zeros((N_Pt,3))
   # 1st layer
   for i in range(3):
     coords[i,0] = a/2 + a*i
     coords[i,1] = 2*v/3
   for i in range(2):
     coords[i+3,0] = a + a*i
     coords[i+3,1] = 2*v/3 + v
   # 2nd layer
   for i in range(3):
     coords[L1+i,0] = a/2 + a*i
   for i in range(4):
     coords[L1+i+3,0] = a*i
     coords[L1+i+3,1] = v
   for i in range(3):
     coords[L1+i+7,0] = a/2 + a*i
     coords[L1+i+7,1] = 2*v
   for i in range(10):
     coords[L1+i,2] = h
   # 3rd layer
   for i in range(2):
     coords[L1+L2+i,0] = a + a*i
     coords[L1+L2+i,1] = v/3
   for i in range(3):
     coords[L1+L2+i+2,0] = a/2 + a*i
     coords[L1+L2+i+2,1] = 4*v/3
   for i in range(5):
     coords[L1+L2+i,2] = 2*h
     
elif config=="9_10_9":
   arr,N_Pt = get_config(config)
   print "Total number of atoms: ",N_Pt
   L1=arr[0]
   L2=arr[1]
   L3=arr[2]
   coords=np.zeros((N_Pt,3))
   # 1st layer
   for i in range(4):
     coords[i,0] = a*i
     coords[i,1] = 5*v/3
   for i in range(3):
     coords[i+4,0] = a/2 + a*i
     coords[i+4,1] = 2*v/3
   for i in range(2):
     coords[i+7,0] = a + a*i
     coords[i+7,1] = -v/3
   # 2nd layer
   coords[L1:L1+L2,:] = get_coords10(a,0,0,h)
   # 3rd layer
   for i in range(4):
     coords[L1+L2+i,0] = a*i
     coords[L1+L2+i,1] = v/3
   for i in range(3):
     coords[L1+L2+i+4,0] = a/2 + a*i
     coords[L1+L2+i+4,1] = 4*v/3
   for i in range(2):
     coords[L1+L2+i+7,0] = a + a*i
     coords[L1+L2+i+7,1] = 7*v/3
   for i in range(9):        # shift vertically
     coords[L1+L2+i,2] = 2*h
   
else: 
  raise NotImplementedError

# ===== get H20 coords
dist=2.00         # in Angstroms
ref_Pt_atom=25
dist_vect = np.array([0.0,0.0,dist])
xyzH2O = np.array( setH2Ocoords(coords[ref_Pt_atom-1,:],dist_vect) )

	
# ===== print into file
printH20(xyzH2O,filename)
printPt(coords,filename)
print coords
print xyzH2O

# ===== freezing coords
if len(sys.argv)>2:
  print len(sys.argv)
  if sys.argv[2]=="freeze":
		f=open(filename,"a")
		for i in range(4,N_Pt+4):
		  f.write(str(i)+" "+"F"+"\n")
		f.write("\n")
		f.close()
else:
  print "No freezing of atoms."
    
    
    
    


