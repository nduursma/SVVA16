# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:43:17 2020

@author: afras
"""
import numpy as np
import matplotlib.pyplot as plt
from math import *



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


def cross_sect_prop(h,tsk,tsp,tst,hst,wst,Ca,nst) :
    i     = 1
    stiff = []
    Ic_z  =[]
    Ic_y  =[]
   
   
    r   = h/2                                  #Radius of the semi circular part of the airfoil in m
    sk  = np.sqrt(r**2 + (Ca-r)**2)      #Length of straight part of the airfoil in m
    phi = atan((r)/(Ca-r))                 #Angle of the non circular part of the airfoil
   
    #At = (pi*r**2)/2+ r*(Ca-r)
    #print(At)
   
   
    A_sk   = sk*tsk                          #Area of the skin straight panel in m2
    A_cir  = pi*r*tsk                       #Area of the skin circular panel in m2
    A_sp   = (h)*tsp                             #Area of the spar in m2
    Ast    =  (wst*tst + (hst-tst)*tst)                 #Area of a stiffener m2
   
    ycirc = 0
    zcirc = -(r-2*(r-tsk)/(pi))  
   
    ysk   = (sk/2)*sin(phi)
    zsk   = -(Ca-((sk/2)*cos(phi)))
   
    ysp   = 0
    zsp   = -r
   
    Skt = [0,zsk,ysk,0,A_sk,A_sk*zsk,A_sk*ysk]
    Skb = [0,zsk,-ysk,0,A_sk,A_sk*zsk,A_sk*(-ysk)]
    Ssp = [0,zsp,ysp,0,A_sp,A_sp*zsp,A_sp*ysp]
    Sci = [0,zcirc,ycirc,0,A_cir,A_cir*zcirc,A_cir*ycirc]
   
    P = np.pi*r +2*sk                           #Calculation of perimeter
   
    Dx = (P)/((nst))                            #Spacing between the stiffeners (We assume that there is one stiffener at the begining of the semi-circular part)
   
    z = 0
    y = 0
   
    stiff.append([i,z,y,0,Ast,Ast*z,Ast*y])
   
   
    loc_stiff = Dx
   
    while loc_stiff < pi*r/2 :              
        i = i + 1
       
        arc_angle = (loc_stiff)/r
       
        z = -(r-r*np.cos(arc_angle))
        y = r*np.sin(arc_angle)
       
        stiff.append([i,z,y,loc_stiff,Ast,Ast*z,Ast*y])
        i = i + 1
        stiff.append([i,z,-y,loc_stiff,Ast,Ast*z,Ast*(-y)])
       
        loc_stiff = loc_stiff + Dx
       
   
    hyp = (pi*r)/2+sk - (loc_stiff)
   
    while loc_stiff< sk+(pi*r)/2 :
        i = i + 1
        z = -(Ca-hyp*cos(phi))
        y = hyp*sin(phi)
        stiff.append([i,z,y,loc_stiff,Ast,Ast*z,Ast*y])
        i = i + 1
        stiff.append([i,z,-y,loc_stiff,Ast,Ast*z,Ast*(-y)])
        loc_stiff = loc_stiff+Dx
        hyp = hyp -Dx
       
   
    graph = np.asarray(stiff)
   
    '''
    if round(2*stiff[10][3]+Dx)== round(P):
        print("ok")
    '''
     
    """
    plt.plot(graph[:,1],graph[:,2],'o')
    plt.axis('equal')
    plt.show()
   
    stiff.append(Skt)
    stiff.append(Skb)
    stiff.append(Ssp)
    stiff.append(Sci)
    """
   
    stiff = np.asarray(stiff)
   
    #print(stiff[:,1],stiff[:,2])
   
    #print(stiff)
   
    ybar = sum(stiff[:,6])/sum(stiff[:,4])
    zbar = sum(stiff[:,5])/sum(stiff[:,4])
   
    #print(stiff[:,5])
   
    z   = (stiff[:,1] -zbar)**2
    y   = (stiff[:,2] -ybar)**2
   
   
    #print(k)
    stiff = np.insert(stiff,7,z,axis=1)
    stiff = np.insert(stiff,8,y,axis=1)
   
   
    Az2 = stiff[:,4]*stiff[:,7]
    Ay2 = stiff[:,4]*stiff[:,8]
   
   
    stiff = np.insert(stiff,9,Az2,axis=1)
    stiff = np.insert(stiff,10,Ay2,axis=1)
   
   
   
    def Icircy(angle):
        Icircy = tsk*r**3*sin(angle)**2
        return Icircy
   
    def Icircz(angle):
        Icircy = tsk*r**3*cos(angle)**2
        return Icircz
   
   
    p = pi/180
    angles = np.arange(0,pi+p,p)
   
    for k in range(len(angles)-1):
        area_z = p*Icircy((angles[k]+angles[k+1])/2)
        area_y = p*Icircy((angles[k]+angles[k+1])/2)
       
        Ic_z.append(area_z)
        Ic_y.append(area_y)
   
    Icircy = sum(Ic_z)
    Icircz = sum(Ic_y)
       
    Ispary = (h*(tsp)**3)/12
    Isparz = (tsp*h**3)/12
   
    Iskz   = (sk**3*tsk*sin(phi)**2)/12
    Isky   = (sk**3*tsk*cos(phi)**2)/12
   
    Izz = Isparz + Icircz + 2*Iskz + sum(stiff[:,10])
    Iyy = Ispary + Icircy + 2*Isky + sum(stiff[:,9]) - A_cir*(4*r/(3*pi))**2
   
    """
    print(zbar)
    print(ybar)
    print("Iyy is", Iyy)
    print("Izz is", Izz)
   
    plt.plot(zbar,ybar,'o')
    plt.axis('equal')
    plt.show()
    """

    return(Dx,zbar,ybar,Izz,Iyy,stiff)

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
   
    y1 = np.cos(theta)*r
    z1 = -r - zc + np.sin(theta)*r
    
   
    q1b = Sy*tsk/Izz*ShearInteg(y1*r,theta)  + Sz*tsk/Iyy*ShearInteg(z1*r,theta)
    
   
    #Beam 2 calculation
    #Note that the shear flow of q1b will be added to q2b after the booms have been added
    #---------------------------------
   
    s2 = np.linspace(0,h,1000)
   
    y2 = -r + s2
    z2 = (-r - zc)*np.ones(1000)
   
    q2b = -Sy*tsp/Izz*ShearInteg(y2,s2) - Sz*tsp/Iyy*ShearInteg(z2,s2) 
   
    #Beam 3&4 calculation
    #-------------------------------------
   
    s3 = np.linspace(0,d,1000)
   
    y3 = s3*r/d
    z3 = -Ca - zc +s3*(Ca-r)/d

    q3b = -Sy*tsk/Izz*ShearInteg(y3,s3) - Sz*tsk/Iyy*ShearInteg(z3,s3)
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
    B3z = np.ones(int(ceil(l1*L)))*(BA3z+BA2z)
    B3y = np.ones(int(ceil(l1*L)))*(BA3y+BA2y)


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
    q4b = q3b
    
    
    q2b = q2b + q1b[-1]

    
    #Setting up matrix for q01 and q02 calculations
    #F1b, F2b etc are integrated function
    
    F1b = ShearInteg(q1b*r,theta)[-1]
    F2b = ShearInteg(q2b, s2)[-1]
    F3b = ShearInteg(q3b, s3)[-1]
    
    A = np.matrix([[(pi*r/tsk+h/tsp),-h/tsp],
                  [(-h/tsp),(2*d/tsk + h/tsp)]])
    
    B = np.matrix([[-F1b/tsk+F2b/tsp],
                  [2*F3b/tsk-F2b/tsp]])
    
    x = np.linalg.solve(A,B)
    
    q01, q02 = x
    
    
    q01 = q01.item()
    q02 = q02.item()
    
    q1 = q1b + q01
    q2 = q2b -q01 + q02
    q3 = q3b - q02
    
    return q1,q2,q3,theta,s2,s3

def CalcShearCenter(h,Ca,tsk,tsp,Iyy,Izz,stiff,dx):
    Sy = 1
    Sz = 0
    
    r = h/2
    d = sqrt(r**2+(Ca-r)**2)
    
    q1,q2,q3,theta,s2,s3 = CalcShearFlow(h,Ca,tsk,tsp,Iyy,Izz,zc,Sy,Sz,stiff,dx)
    
    F1 = ShearInteg(q1*r,theta)[-1]
    F3 = ShearInteg(q3,s3)[-1]
    
    Mtot = F1*r - 2*F3*d
    eta = -Mtot
    Scz = -r - eta
    
    return SCz

def CalcDirectStress(My,Mz,Izz,Iyy,zc,theta,s2,s3,s4):
    
    #Defining useful geometry
    r = h/2
    d = sqrt(r**2+(Ca-r)**2)
   
    #Bending stress in beam 1
    y1 = np.cos(theta)*r
    z1 = -r - zc + np.sin(theta)*r
    
    simgax1 = My*(z1-zc)/Iyy + Mz*y1/Izz
    
    #Bending stress in beam 2

    y2 = -r + s2
    z2 = (-r - zc)*np.ones(1000)
    
    sigmax2 = My*(z2-zc)/Iyy + Mz*y2/Izz
   
    #Bending stress in beam 3
    
    y3 = s3*r/d
    z3 = -Ca - zc +s3*(Ca-r)/d
    
    sigmax3 = My*(z3-zc)/Iyy + Mz*y3/Izz
    
    #Bending stress in beam 4
    
    y4 = -s4*r/d
    z4 = -Ca - zc +s4*(Ca-r)/d
    
    sigmax4 = My*(z4-zc)/Iyy + Mz*y4/Izz
    
    return sigmax1

#def CalcVonMises():
    
"""
x = cross_sect_prop(h,tsk,tsp,tst,hst,wst,Ca,nst)

dx = x[0]
zc = x[1]
Izz = x[3]
Iyy = x[4]
stiff = x[5]
Sy = 1
Sz = 0

q1, q2, q3, theta,s2,s3 = CalcShearFlow(h,Ca,tsk,tsp,Iyy,Izz,zc,Sy,Sz,stiff,dx)

#plt.plot(s3,q3)
#plt.show()

print(CalcShearCenter(h,Ca,tsk,tsp,Iyy,Izz,stiff,dx))
plt.plot(s3,q3)
plt.show()
"""

