# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from Interpolation import patchinterpolate



#Shear Center 
sc = -0.085

data = np.loadtxt('AERO.dat',delimiter = ',')
#xlst, zlst, qlst = patchinterpolate(600,600,data)

def output(xlst, zlst, qlst, sc):
    
    Nz = len(qlst)
    Nx = len(qlst[0])
    dz = 1/(Nz-1)*(zlst[-1]-zlst[0])
    dx = 1/(Nx-1)*(xlst[-1]-xlst[0])
    
    qavglst = [] 
    taulst = []
    
    qtot = 0
    
    for xj in range(0,Nx): 
        Az = 0
        tau = 0
        for zj in range(0,Nz):
            Azi = (qlst[zj][xj] + qlst[zj-1][xj])/2*dz
            Az = Az + Azi
            
            taui = (qlst[zj][xj]*(zlst[zj] - sc) + qlst[zj-1][xj]*(zlst[zj-1] - sc))/2*dz
            tau = tau + taui
            
            
            qi = (qlst[zj][xj] + qlst[zj][xj-1]+qlst[zj-1][xj] + qlst[zj-1][xj-1])/4*dz*dx
            qtot = qi + qtot
            
            
        qavglst.append(Az)
        taulst.append(tau)
    
    defl = 0
    defllst = []
    
    V = 0
    Vlst = []
    
    M = 0
    Mlst = []
    
    d1 = 0
    d1lst = []
    
    Tlst = []
    T = 0
    
    thetalst = []
    theta = 0
    
    for p in range(len(qavglst)):
        V = V + (qavglst[p-1]+qavglst[p])*dx/2
        Vlst.append(V)
        
        #M = M + (qavglst[p-1]*((p-1)*dx)+qavglst[p]*(p*dx))*dx/2
        M = M + (Vlst[p-1]+Vlst[p])*dx/2
        Mlst.append(M)
  
        d1 =  d1 + (Mlst[p-1]+Mlst[p])*dx/2
        d1lst.append(d1)
        
        defl = defl + (d1lst[p-1]+d1lst[p])*dx/2
        defllst.append(defl)
        
        T = T + (taulst[p-1]+taulst[p])*dx/2
        Tlst.append(T)
        
        theta = theta +  (Tlst[p-1]+Tlst[p])*dx/2
        thetalst.append(theta)

    Vlst = np.array(Vlst)*1000        # [N]
    Mlst = np.array(Mlst)*1000        # [Nm]
    d1lst = np.array(d1lst)*1000      # [/m]
    defllst = np.array(defllst)*1000  # [m]
    taulst = np.array(taulst)*1000    # [Nm/m]
    Tlst = np.array(Tlst)*1000        # [Nm]

    return Vlst, Mlst, d1lst, defllst, taulst, Tlst, thetalst

#Vlst, Mlst, defllst, Tlst, thetalst = output(sc)
#print(' V', Vlst[-1], '\n M_z', Mlst[-1], '\n Deflection', defllst[-1],'/EI \n Torque', Tlst[-1], '\n Theta', thetalst[-1],'/GJ')

def locationvalue(xlst, xloc, reflst):
    dx = 1/(len(xlst)-1)*xlst[-1]
    n = 0 
    while dx*n <= xloc:
        n = n+1
    if n>=len(reflst):
        value = reflst[-1]
    else:
        value = reflst[n]
    return value





 











    
    
