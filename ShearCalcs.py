# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:43:17 2020

@author: afras
"""
import numpy as np
import matplotlib.pyplot as plt
from math import *
import cross_sect_prop


Ca    = 0.505       #Chord Length Aileron
la    = 1.611       #Span of the Aileron
x1    = 0.125       #x-location of hinge 1
x2    = 0.498       #x-location of hinge 2
x3    = 1.494       #x-location of hinge 3
xa    = 24.5        #Distance between actuator 1 and 2
h     = 0.161        #Aileron Height
tsk   = 1.1/1000         #Skin thickness
tsp   = 2.4/1000         #Spar Thickness
tst   = 1.2/1000         #Thickness of stiffener
hst   = 1.3/100         #Height of stifffener
wst   = 1.7/100         #Width of stiffeners
nst   = 11          #Number of stiffeners (equalt spaced)

x = cross_sect_prop(h,tsk,tsp,tst,hst,wst,Ca,nst)

def CalcShearTorq(Ca, h, tsk, tsp,T):
    #Defining geometrical values
   
    r = 1/2 * h
    d = sqrt((Ca-r)**2+r**2)
    A1 = 1/2 * pi * r**2
    A2 = (Ca-r)*r
   
    #Coefficient matrix
   
    a11 = 1/(2*A1)*(pi*r/tsk+h/tsp)
    a12 = 1/(2*A1)*(-h/tsp)
    a13 = -1
    a21 = 1/(2*A2)*(-h/tsp)
    a22 = 1/(2*A2)*(2*d/tsk + h/tsp)
    a23 = -1
    a31 = 2*A1
    a32 = 2*A2
    a33 = 0
   
    A = np.matrix([[a11,a12,a13],
                   [a21,a22,a23],
                   [a31,a32,a33]])
   
    b = np.matrix([[0],
                   [0],
                   [T]])
   
    x = np.linalg.inv(A)*b
     
    return x


 
def CalcTorsStif(Ca,h,tsk,tsp):
   
    T = 1
    q1, q2, Gdethadx  = CalcShearTorq(Ca, h, tsk, tsp, T)
    J = T/(Gdethadx)
    J = J.item()
    return J

def ShearInteg(x,y):
   
    znew = 0
    z = [0]
    for i in range(len(x)-1):
        dz = (x[i+1] + x[i]) *0.5
        znew = znew + dz*(y[-1])/(len(x)-1)
        z.append(znew)
   
    return np.asarray(z)

def CalcShearFlow(h,Ca,tsk,tsp,Iyy,Izz,zc,Sy,Sz,stiff,dx):
   
    #Defining distances
    r = h/2
    d = sqrt(r**2+(Ca-r)**2)
   
    #Beam 1 calculation
    #-------------------------------
   
    theta = np.linspace(0,pi,1000)
   
    y1 = np.cos(theta)
    x1 = -r - zc + np.sin(theta)
   
    q1b = Sy*tsk/Izz*ShearInteg(y1*r,theta)  + Sz*tsk/Iyy*ShearInteg(x1*r,theta)
   
    #Beam 2 calculation
    #Note that the shear flow of q1b will be added to q2b after the booms have been added
    #---------------------------------
   
    s2 = np.linspace(0,h,1000)
   
    y2 = -r + s2
    x2 = (-r - zc)*np.ones(1000)
   
    q2b = -Sy*tsk/Izz*ShearInteg(y2,s2) - Sz*tsk/Iyy*ShearInteg(x2,s2) 
   
    #Beam 3&4 calculation
    #-------------------------------------
   
    s3 = np.linspace(0,h,1000)
   
    y3 = s3*r/d
    x3 = -Ca - zc +s3*(Ca-r)/d

    q3b = -Sy*tsk/Izz*ShearInteg(y3,s3) - Sz*tsk/Iyy*ShearInteg(x3,s3)
   
    #Beam 4 calculation
    #-------------------------------------
   
    s4 = np.linspace(0,h,1000)
   
    q4b = q3b
   
    #Adding effect of relevant booms z direction
   
    BA2z = stiff[1][5]
    BA3z = stiff[2][5]
    BA4z = stiff[3][5]
    BA6z = stiff[5][5]
    BA8z = stiff[7][5]
    BA10z = stiff[9][5]
   
    #booms y direction
   
    BA2y = stiff[1][6]
    BA3y = stiff[2][6]
    BA4y = stiff[3][6]
    BA6y = stiff[5][6]
    BA8y = stiff[7][6]
    BA10y = stiff[9][6]
   
    #vector with application points booms beam 1
   
    l1 = (0.5*pi*r-dx)/(pi*r)
    l3 = 1 - l1

    L = len(q1b)

    B1 = np.zeros(int(round(l1*L)))
    B2z = np.ones(int(round((l3-l1)*L)))*BA2z
    B2y = np.ones(int(round((l3-l1)*L)))*BA2y
    B3z = np.ones(int(ceil(l1*L)))*BA3z
    B3y = np.ones(int(ceil(l1*L)))*BA3y


    Bsty = np.append(B1,B2y)
    Bsty1 = np.append(Bsty,B3y)
   
    Bstz = np.append(B1,B2z)
    Bstz1 = np.append(Bstz,B3z)
   
   
   
    #vector with application points beam 3
   
    l1 = 0.5*dx/d
    l2 = dx/d
    L = len(q3b)
   
    B1 = np.zeros(int(round(l1*L)))
    B2z = np.ones(int(round(l2*L)))*BA10z
    B2y = np.ones(int(round(l2*L)))*BA10y
    B3z = np.ones(int(round(l2*L)))*(BA10z+BA8z)
    B3y = np.ones(int(round(l2*L)))*(BA10y+BA8y)
    B4y = np.ones(int(floor(l2*L)))*(BA10y+BA8y+BA6y)
    B4z = np.ones(int(floor(l2*L)))*(BA10z+BA8z+BA6z)
    B5y = np.ones(int(round((1-l1-l2*3)*L)))*(BA10y+BA8y+BA6y+BA4y)
    B5z = np.ones(int(round((1-l1-l2*3)*L)))*(BA10z+BA8z+BA6z+BA4z)
   

    Bsty = np.append(B1,B2y)
    Bsty = np.append(Bsty,B3y)
    Bsty = np.append(Bsty,B4y)
    Bsty3 = np.append(Bsty,B5y)
   
    Bstz = np.append(B1,B2z)
    Bstz = np.append(Bstz,B3z)
    Bstz = np.append(Bstz,B4z)
    Bstz3 = np.append(Bstz,B5z)
    
    #Adding booms to base shear flows
    
    q1b = q1b + Sy/Izz*Bsty1 + Sz/Iyy*Bstz1
    q2b = q2b + q1b[-1]
    q3b = q3b - Sy/Izz*Bsty3 - Sz/Iyy*Bstz3
    q4b = -q3b
    
    return q1b,q2b,q3b,q4b,theta,s2,s3
"""    
dx = x[0]
zc = x[1]
Izz = x[3]
Iyy = x[4]
stiff = x[5]
Sy = 1
Sz = 0
   
q1b, q2b, q3b, q4b, theta, s2,s3 = CalcShearFlow(h,Ca,tsk,tsp,Iyy,Izz,zc,Sy,Sz,stiff,dx)

plt.plot(theta,q1b)
plt.show()
"""   