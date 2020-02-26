# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 16:23:31 2020

@author: afras
"""

import numpy as np
import matplotlib.pyplot as plt
from math import *
from Cross_sectional_properties import cross_sect_prop

Dx, zbar,ybar,Izz,Iyy,stiff = cross_sect_prop(h,tsk,tsp,tst,hst,wst,Ca,nst)


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
    q1, q2, GdethaDx  = CalcShearTorq(Ca, h, tsk, tsp, T)
    J = T/(GdethaDx)
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

def CalcShearFlow(h,Ca,tsk,tsp,Iyy,Izz,zbar,Sy,Sz,stiff,Dx):
   
    #Defining distances
    r = h/2
    d = sqrt(r**2+(Ca-r)**2)
   
    #Beam 1 calculation
    #-------------------------------
   
    theta = np.linspace(0,pi,1000)
   
    y1 = np.cos(theta)*r
    z1 = -r - zbar + np.sin(theta)*r
   
   
    q1b = Sy*tsk/Izz*ShearInteg(y1*r,theta)  + Sz*tsk/Iyy*ShearInteg(z1*r,theta)

   
    #Beam 2 calculation
    #Note that the shear flow of q1b will be added to q2b after the booms have been added
    #---------------------------------
   
    s2 = np.linspace(0,h,1000)
   
    y2 = -r + s2
    z2 = (-r - zbar)*np.ones(1000)
   
    q2b = -Sy*tsp/Izz*ShearInteg(y2,s2) - Sz*tsp/Iyy*ShearInteg(z2,s2)

    #Beam 3&4 calculation
    #-------------------------------------
   
    s3 = np.linspace(0,d,1000)
   
    y3 = s3*r/d
    z3 = -Ca - zbar +s3*(Ca-r)/d

    q3b = -Sy*tsk/Izz*ShearInteg(y3,s3) - Sz*tsk/Iyy*ShearInteg(z3,s3)
    #Beam 4 calculation
    #-------------------------------------
   
    s4 = s3
   
    q4b = -Sy*tsk/Izz*ShearInteg(y3,s3) + Sz*tsk/Iyy*ShearInteg(z3,s3)
   
    #Adding effect of relevant booms z direction
   
    BA2z = stiff[1][5] - stiff[1][4]*zbar
    BA3z = stiff[2][5] - stiff[2][4]*zbar
    BA4z = stiff[3][5] - stiff[3][4]*zbar
    BA6z = stiff[5][5] - stiff[5][4]*zbar
    BA8z = stiff[7][5] - stiff[7][4]*zbar
    BA10z = stiff[9][5] - stiff[9][4]*zbar
   
    #booms y direction
   
    BA2y = stiff[1][6]
    BA3y = stiff[2][6]
    BA4y = stiff[3][6]
    BA6y = stiff[5][6]
    BA8y = stiff[7][6]
    BA10y = stiff[9][6]
   
    #vector with application points booms beam 1
   
    l1 = (0.5*pi*r-Dx)/(pi*r)
    l3 = 1 - l1
     
    L = len(q1b)
      
    B1 = np.zeros(int(round(l1*L)))
    B2z = np.ones(int(round((l3-l1)*L)))*BA2z
    B2y = np.ones(int(round((l3-l1)*L)))*BA2y
    B3z = np.ones(int(ceil(l1*L)))*(BA3z+BA2z)
    B3y = np.ones(int(ceil(l1*L)))*(BA3y+BA2y)


    Bsty = np.append(B1,B2y)
    Bsty1 = np.append(Bsty,B3y)
   
    Bstz = np.append(B1,B2z)
    Bstz1 = np.append(Bstz,B3z)
   
   
    #vector with application points beam 3
   
    l1 = 0.5*Dx/d
    l2 = Dx/d
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

    q3b = q3b - Sy/Izz*Bsty3 - Sz/Iyy*Bstz3
    q4b = q4b - Sy/Izz*Bsty3 + Sz/Iyy*Bstz3
   
    q2b = q2b - q1b[-1] - q4b[-1]
    #Setting up matrix for q01 and q02 calculations
    #F1b, F2b etc are integrated function
   
   
    F1b = ShearInteg(q1b*r,theta)[-1]
    F2b = ShearInteg(q2b, s2)[-1]
    F3b = ShearInteg(q3b, s3)[-1]
    F4b = ShearInteg(q4b, s4)[-1]

    A = np.matrix([[pi*r/tsk+h/tsp,-h/tsp],
                  [-h/tsp,2*d/tsk + h/tsp]])
   
    B = np.matrix([[-F1b/tsk+F2b/tsp],
                  [(F3b+F4b)/tsk-F2b/tsp]])
   
    x = np.linalg.solve(A,B)
   
    q01, q02 = x
   

    q01 = q01.item()
    q02 = q02.item()
   
    q1 = q1b + q01
    q2 = q2b -q01 + q02
    q3 = q3b - q02
    q4 = q4b - q02
   

    return q1,q2,q3,q4,theta,s2,s3,s4

def CalcShearCenter(h,Ca,tsk,tsp,Iyy,Izz,stiff,Dx):
    Sy = 1
    Sz = 0
   
    r = h/2
   
    q1,q2,q3,q4,theta,s2,s3,s4 = CalcShearFlow(h,Ca,tsk,tsp,Iyy,Izz,zbar,Sy,Sz,stiff,Dx)
   
    F1 = ShearInteg(q1*r,theta)[-1]
    F3 = ShearInteg(q3,s3)[-1]
    F4 = ShearInteg(q4,s3)[-1]
   
    phi = atan(r/(Ca-r))
    L = cos(phi)*r
   
    Mtot = F1*r - (F3+F4)*L
    eta = -Mtot
    SCz = -r - eta
    return SCz

def CalcDirectStress(My,Mz,Izz,Iyy,zbar,theta,s2,s3,s4):
   
    #Defining useful geometry
    r = h/2
    d = sqrt(r**2+(Ca-r)**2)
   
    #Bending stress in beam 1
    y1 = np.cos(theta)*r
    z1 = -r - zbar + np.sin(theta)*r
   
    sigmax1 = My*(z1-zbar)/Iyy + Mz*y1/Izz
    
    #Bending stress in beam 2

    y2 = -r + s2
    z2 = (-r - zbar)*np.ones(1000)
   
    sigmax2 = My*(z2-zbar)/Iyy + Mz*y2/Izz
   
    #Bending stress in beam 3
   
    y3 = s3*r/d
    z3 = -Ca - zbar +s3*(Ca-r)/d
   
    sigmax3 = My*(z3-zbar)/Iyy + Mz*y3/Izz
   
    #Bending stress in beam 4
   
    y4 = -s4*r/d
    z4 = -Ca - zbar +s4*(Ca-r)/d
   
    sigmax4 = My*(z4-zbar)/Iyy + Mz*y4/Izz
   
    return sigmax1, sigmax2, sigmax3, sigmax4



def CalcAllVonMises(tsk,tsp,sigmax1,sigmax2,sigmax3,sigmax4,q1,q2,q3,q4):
   
    #Use Von Mises formula for all beams
    
    vonmises1 = np.sqrt(np.square(sigmax1)+3*(np.square(q1/tsk)))
    vonmises2 = np.sqrt(np.square(sigmax2)+3*(np.square(q2/tsp)))
    vonmises3 = np.sqrt(np.square(sigmax3)+3*(np.square(q3/tsk)))
    vonmises4 = np.sqrt(np.square(sigmax4)+3*(np.square(q4/tsk)))
    
   
    return vonmises1, vonmises2, vonmises3, vonmises4
   

def CalcCombineShear(h,Ca,tsk,tsp,Iyy,Izz,zbar,Sy,Sz,stiff,Dx,T):
    
    q1, q2, q3,q4,theta,s2,s3,s4 = CalcShearFlow(h,Ca,tsk,tsp,Iyy,Izz,zbar,Sy,Sz,stiff,Dx)
    
    qt1, qt2, Gdethetadz = CalcShearTorq(Ca, h, tsk, tsp,T)
    
    q1 = q1 + qt1.item()
    q2 = q2 -qt1.item() + qt2.item()
    q3 = q3 - qt2.item()
    q4 = q4 - qt2.item()
    
    return q1,q2,q3,q4,theta,s2,s3,s4

def CreateStressArrays(Ca,h,theta,s2,s3,s4,q1,q2,q3,q4,sigmax1,sigmax2,sigmax3,sigmax4,vonmises1,vonmises2,vonmises3,vonmises4):
    
    r = h/2
    d = sqrt((Ca-r)**2+r**2)
    
    qtot = np.array([])
    sigmaxtot = np.array([])
    vonmisestot = np.array([])
    ztot = np.array([])
    ytot = np.array([])
    
    qtot = np.append(q1,q2)
    qtot = np.append(qtot,q3)
    qtot = np.append(qtot,q4)

    sigmaxtot = np.append(sigmax1, sigmax2)
    sigmaxtot = np.append(sigmaxtot, sigmax3)
    sigmaxtot = np.append(sigmaxtot, sigmax4)

    vonmisestot = np.append(vonmises1, vonmises2)
    vonmisestot = np.append(vonmisestot, vonmises3)
    vonmisestot = np.append(vonmisestot, vonmises4)

    y1 = np.cos(theta)*r
    z1 = -r + np.sin(theta)*r

    y2 = -r + s2
    z2 = (-r)*np.ones(1000)

    y3 = s3*r/d
    z3 = -Ca +s3*(Ca-r)/d 
    
    y4 = -s4*r/d
    z4 = -Ca +s4*(Ca-r)/d

    ztot = np.append(z1,z2)
    ztot = np.append(ztot,z3)
    ztot = np.append(ztot,z4)

    ytot = np.append(y1,y2)
    ytot = np.append(ytot,y3)
    ytot = np.append(ytot,y4)
    
    return ytot,ztot,qtot,sigmaxtot,vonmisestot

def MaxStress(vonmisestot,ytot,ztot):
    i = np.argmax(vonmisestot)
    zmax = ztot[i]
    ymax = ytot[i]
    maxstress = vonmisestot[i]
    
    return zmax, ymax, maxstress
    

def AileronFullStress(Sy,Sz,My,Mz,T,Ca,h,la,tsp,tsk,Iyy,Izz,zbar,Dx):
    #define x
    x = np.linspace(0,la,step)
    
    maxstresslist = []
    zcoordlist = []
    ycoordlist = []
    xcoordlist = []
    
    for i in np.arange(len(x)):
        #Find shear flows 
        q1,q2,q3,q4,theta,s2,s3,s4 = CalcCombineShear(h,Ca,tsk,tsp,Iyy,Izz,zbar,Sy[i],Sz[i],stiff,Dx,T[i])
        
        #Find bending stresses
        sigmax1, sigmax2, sigmax3, sigmax4 = CalcDirectStress(My[i],Mz[i],Izz,Iyy,zbar,theta,s2,s3,s4)
        
        #Find vonmises
        vonmises1, vonmises2, vonmises3, vonmises4 = CalcAllVonMises(tsk,tsp,sigmax1,sigmax2,sigmax3,sigmax4,q1,q2,q3,q4)
        
        #Create stress arrays for each cut
        ytot,ztot,qtot,sigmaxtot,vonmisestot = CreateStressArrays(Ca,h,theta,s2,s3,s4,q1,q2,q3,q4,sigmax1,sigmax2,sigmax3,sigmax4,vonmises1,vonmises2,vonmises3,vonmises4) 
        
        #Store each maximum value in list
        zmax,ymax,maxstress = MaxStress(vonmisestot,ytot,ztot)
        
        #Find highes stress in entire aileron
        maxstresslist.append(maxstress)
        zcoordlist.append(zmax)
        ycoordlist.append(ymax)
        xcoordlist.append(x[i])
    
    imax = np.argmax(maxstresslist)
    
    maximumstress = maxstresslist[imax]
    zloc = zcoordlist[imax]
    yloc = ycoordlist[imax]
    xloc = xcoordlist[imax]
    
    return maximumstress,zloc,yloc,xloc

#-------------------------------------------------------------------------
#                           Load input
#-------------------------------------------------------------------------    
"""
step = 50

Sy = np.linspace(0,2000,step)
Sz = np.linspace(0,12352,step)
My = np.linspace(0,5000,step)
Mz = np.linspace(0,2342,step)
T = np.linspace(0,2424,step)

maxst,zloc,yloc,xloc = AileronFullStress(Sy,Sz,My,Mz,T,Ca,h,la,tsp,tsk,Iyy,Izz,zbar,Dx)

print('Maximum stress = ', maxst, 'N/m^2')
print('locatiion of max stress = ',xloc,yloc,zloc,'(x,y,z)')
"""