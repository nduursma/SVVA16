# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import pi 
import numpy as np
from Interpolation import patchinterpolate



#Read data of distributed load magnitude
data = np.loadtxt('AERO.dat',delimiter = ',')  #Creating an Array from the aerodynamic load data file
xlst, zlst, y = patchinterpolate(20,20,data)
qlst = y


"""
#Create coordinate lists in x and z direction
xlst = []
zlst = []
"""

ca = 0.505 #m
la = 1.611 #m

#Calculate locations of X and Z coordinates
Nz = len(qlst)
Nx = len(qlst[0])

"""
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
""" 



#Create list of distributed load values where the rows are along the z direction, and the columns are along the x direction 
#Like a transpose of the data list, where columns are along the z direction and rows along the x direction

    


#Create list with distributed load for every slice in x direction [N/m], and its centroid in Z direction
Centroid_Zlst = [] 
Maglst = []


#Define integration steps dx and dz
dz = 1/(Nz-1)
dx = 1/(Nx-1)

#Vtot = magnitude of the resultant force caused by the TOTAl distributed load [N]. 
#This resultant force has been computed in two different ways, therefore Vtot and Vtot2 have been defined.
Vtot = 0
Vtot2 = 0

#Slice in the X direction
for xj in range(0,Nx): 

#Az = Area below distributed load graph for every slice in x direction [N/m]
    Az = 0
    
#Vz = Volume / magnitude of resultant force for every slice in x direction [N]
    Vz = 0  
    
#Q_y = first moment of area in y direction 
    Q_y = 0

    for zj in range(0,Nz):
        
#Using the trapezoidal method for double integration in the X and Z direction (to verify total resultant force)
        Vxzi = (qlst[zj][xj] + qlst[zj][xj-1]+qlst[zj-1][xj] + qlst[zj-1][xj-1])/4*dz*dx
        Vtot = Vtot + Vxzi
  
#Using the trapezoidal method for single integration in the Z direction        
        Azi = (qlst[zj][xj] + qlst[zj-1][xj])/2*dz
        Vzi = Azi*dx
        Az = Az + Azi
        Vz = Vz + Vzi
        Mag = Vz
        
       
#First moment of area             
        Qi = Azi*zlst[zj]
        Q_y = Qi + Q_y

#Resultant force        
    Vtot2 = Vtot2 + Az*dx
   
#Centroid
    Centroid_Z = Q_y/Az    
    Centroid_Zlst.append(Centroid_Z)

#Magnitude
    Maglst.append(Mag)
    
    
    
print('Magnitude of resultant force for verification:', Vtot)
print('Magnitude of resultant force computed:', Vtot2)
print('Centroids:', Centroid_Zlst)     
print('Magnitudes:', Maglst)      
        

    
    
    
    