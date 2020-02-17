# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import pi 
import numpy as np
#Create coordinate lists in x and z direction
xlst = []
zlst = []


ca = 0.505 #m
la = 1.611 #m

#Read data of distributed load magnitude
data = np.loadtxt('AERO.dat',delimiter = ',')  #Creating an Array from the aerodynamic load data file

#Calculate locations of X and Z coordinates
Nz = len(data)
Nx = len(data[0])
print('Nx:', Nx)
print('Nz:', Nz)
for i in range(1,Nz+1):
    
    thetazi  = (i-1)/Nz*pi
    thetazi1 = (i)/Nz*pi
    z  =  -(1/2)*((ca/2)*(1-np.cos(thetazi))+(ca/2)*(1-np.cos(thetazi1)))    
    zlst.append(z)
    
for k in range(1,Nx+1):
    thetaxk  = (k-1)/Nx*pi
    thetaxk1 = (k)/Nx*pi
    x  = (1/2)*((la/2)*(1-np.cos(thetaxk))+(la/2)*(1-np.cos(thetaxk1)))   
    xlst.append(x)

#Create list of distributed load values where the rows are along the z direction, and the columns are along the x direction 
#Like a transpose of the data list, where columns are along the z direction and rows along the x direction
qlst = []
for xi in range(0,Nx):
    qlsti = []
    for zi in range(0, Nz):  
        qvaluez = data[zi][xi]
        qlsti.append(qvaluez)
    qlst.append(qlsti)
    

  
#Define integration step dz
dz = 1/(len(zlst)-1)

#Define resultant magnitide and point of application. The list contains [Magnitude, Z coordinate Centroid] arrays
MagCentroidlst = []    


#Slice in the X direction
for xj in range(0,Nx): 
    Atot = 0
    Aztot = 0
    
    for zj in range(0,len(zlst)):
        
#Using the trapezoidal method for integration in the Z direction
        Areai = (qlst[xj][zj] + qlst[xj][zj-1])/2*dz
        Atot = Atot + Areai
        
#Find first moment of Area to compute centroid
        Az = Areai*zlst[zj]
        Aztot = Aztot + Az
        
#Find magnitude and centroid of resultant force
    Magq = Atot    
    CentroidZ = Aztot/Atot
    Matrix = [Magq, CentroidZ]  
    
#Add to list with results. The list contains [Magnitude, Z coordinate centroid] arrays for x coordinates starting at x=0 till x=1.609 
    MagCentroidlst.append(Matrix)

print(MagCentroidlst)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    