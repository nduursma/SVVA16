# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import pi 
import numpy as np
from Interpolation import patchinterpolate

#z coordinate of the hinge line
c_a = 0.505
h_a = 0.161
hl_z = -c_a + h_a/2


#Read data of distributed load magnitude
data = np.loadtxt('AERO.dat',delimiter = ',')  #Creating an Array from the aerodynamic load data file
#xlst, zlst, qlst = patchinterpolate(200,200,data)


#Outputs a list with the total magnitude of the resulant force, the magnitude per slice in x direction, and its belonging centroid.
def Magnitude_Centroid(x_mesh,z_mesh,data):
    #Calculate locations of X and Z coordinates
    xlst, zlst, qlst = patchinterpolate(x_mesh,z_mesh,data)
    Nz = len(qlst)
    Nx = len(qlst[0])
    #Create list with distributed load for every slice in x direction [N], and its centroid in Z direction
    Centroid_Zlst = []
    Maglst = []
    
    
    #Define integration steps dx and dz
    dz = 1/(Nz-1)*zlst[-1]
    dx = 1/(Nx-1)*xlst[-1]
    
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
              
    """
    print('Centroids [m]:', Centroid_Zlst)     
    print('Magnitudes Resultant Force for every slice [N]:', Maglst) 
    print('Magnitude of resultant force for verification:', Vtot)
    print('Magnitude of resultant force computed:', Vtot2)    
    """
    
    return Vtot, Maglst, Centroid_Zlst     


#Find the Torque around the x axis for a specific x location 
def Torque_x(xlst, Maglst, Centroid_Zlst, xloc):
    T = 0 
    dx = 1/(len(xlst)-1)*xlst[-1]
    n = 0 
    while dx*n <= xloc:
        n = n+1        
    for m in range(n):
        T= T+Maglst[m]*(hl_z - Centroid_Zlst[m])
    return T

#Outputs an array with torque caused by distributed load up till cut in x direction
def Torque_xarray(Maglst, Centroid_Zlst):
    Tlst = []
    T = 0
    for m in range(len(Maglst)):
        T= T + Maglst[m]*(hl_z - Centroid_Zlst[m])
        Tlst.append(T)
    return Tlst

#Find the moment around the z axis for a specific x location
def Moment_z(xlst, Maglst, xloc):
    dx = 1/(len(xlst)-1)*xlst[-1]
    n = 0 
    while dx*n <= xloc:
        n = n+1        
    Mz = 0
    for r in range(n):
        Mz = Mz + dx*r*Maglst[r]
    return Mz

#Outputs an array with the moment around the z axis caused by the distributed load up till cut in x direction.
def Moment_zarray(Maglst):
    dx = 1/(len(xlst)-1)*xlst[-1]
    Mz = 0
    Mzlst = []
    for r in range(len(Maglst)):
        Mz = Mz + dx*r*Maglst[r]
        Mzlst.append(Mz)
    return Mzlst



Vtot, Maglst, Centroid_Zlst = Magnitude_Centroid(xlst,zlst,qlst)
Tlst = Torque_xarray(Maglst, Centroid_Zlst)
Mzlst = Moment_zarray(Maglst)















    
    
