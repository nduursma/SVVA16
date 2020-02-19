# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Sun May 26 15:32:55 2019

@author: nadie
"""
import numpy as np
import matplotlib.pyplot as plt
from math import pi 
from mpl_toolkits.mplot3d import Axes3D

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

print(Nz)
for i in range(1,Nz+1):
    
    #print(i)
    
    thetazi  = (i-1)/Nz*pi
    thetazi1 = (i)/Nz*pi
    z  =  -(1/2)*((ca/2)*(1-np.cos(thetazi))+(ca/2)*(1-np.cos(thetazi1)))    
    zlst.append(z)
    
for k in range(1,Nx+1):
    thetaxk  = (k-1)/Nx*pi
    thetaxk1 = (k)/Nx*pi
    x  = (1/2)*((la/2)*(1-np.cos(thetaxk))+(la/2)*(1-np.cos(thetaxk1)))   
    xlst.append(x)
      
def CreatePlots(xlst,zlst,qlst):  
            
    #Create 2D Contour plot
    X, Z = np.meshgrid(xlst, zlst)
    Y = data
    
    plt.figure()
    cp = plt.contourf(X, Z, Y )
    plt.colorbar(cp)
    
    plt.title('Contour Plot Distributed Load on Airleron [kN/m^2]')
    plt.xlabel('X along wingspan [m]')
    plt.ylabel('Z along chord [m]')
    plt.show()
    
    
    #Create 3D plot
    ax = plt.axes(projection='3d')
    
    ax.plot_surface(X, Z, Y, rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none')
    
    
    plt.title('3D Plot Distributed Load on Aileron [kN/m^2]')
    ax.set_xlabel('X along wingspan [m]')
    ax.set_ylabel('Z along chord [m]')
    ax.set_zlabel('Distributed Load [kN/m^2]')
    ax.view_init(azim=166)
    plt.show()
        
CreatePlots(xlst,zlst,data)   
    
