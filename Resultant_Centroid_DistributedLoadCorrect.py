# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from Interpolation import patchinterpolate



#Shear Center 
sc = -0.22 


def output(xlst, zlst, qlst, sc):
    
    Nz = len(qlst)
    Nx = len(qlst[0])
    dz = 1/(Nz-1)*zlst[-1]
    dx = 1/(Nx-1)*xlst[-1]
    
    qavglst = [] 
    taulst = []
    
    
    for xj in range(0,Nx): 
        Az = 0
        tau = 0
        for zj in range(0,Nz):
            Azi = (qlst[zj][xj] + qlst[zj-1][xj])/2*dz
            Az = Az + Azi
            
            taui = (qlst[zj][xj]*(zlst[zj] - sc) + qlst[zj-1][xj]*(zlst[zj-1] - sc))/2*dz
            tau = tau + taui
            
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
        
        M = M + (qavglst[p-1]*((p-1)*dx)+qavglst[p]*(p*dx))*dx/2
        Mlst.append(M)
  
        d1 =  d1 + (qavglst[p-1]*(((p-1)*dx)**2)/2+qavglst[p]*((p*dx)**2)/2)*dx/2
        d1lst.append(d1)
        
        defl = defl + (qavglst[p-1]*(((p-1)*dx)**3)/3+qavglst[p]*((p*dx)**3)/3)*dx/2
        defllst.append(defl)
        
        T = T + (taulst[p-1]+taulst[p])*dx/2
        Tlst.append(T)
        
        theta = theta +  (taulst[p-1]*((p-1)*dx)+taulst[p]*(p*dx))*dx/2
        thetalst.append(theta)
        
    return Vlst, Mlst, defllst, Tlst, thetalst


data = np.loadtxt('AERO.dat',delimiter = ',') 
xlst, zlst, qlst = patchinterpolate(100,100,data)
Vlst, Mlst, defllst, Tlst, thetalst = output(xlst, zlst, qlst, sc)
print(' V', Vlst[-1], '\n M_z', Mlst[-1], '\n Deflection', defllst[-1],'/EI \n Torque', Tlst[-1], '\n Theta', thetalst[-1],'/GJ')    




 











    
    