# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 11:57:00 2020

@author: Dimitris
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:16:02 2020

@author: Dimitris
"""
from math import * 
import numpy as np  
from matplotlib import pyplot as plt

Ca    = 0.505       #Chord Length Aileron 
la    = 1.611       #Span of the Aileron
x1    = 0.125       #x-location of hinge 1
x2    = 0.498       #x-location of hinge 2 
x3    = 1.494       #x-location of hinge 3
xa    = 24.5        #Distance between actuator 1 and 2
h     = 16.1        #Aileron Height 
tsk   = 1.1         #Skin thickness
tsp   = 2.4         #Spar Thickness
tst   = 1.2         #Thickness of stiffener
hst   = 1.3         #Height of stifffener
wst   = 1.7         #Width of stiffeners 
nst   = 11          #Number of stiffeners (equalt spaced)
d1    = 0.389       #Vertical displacement of hinge 1
d3    = 1.245       #Vertical displacement of hinge 2
theta = 30          #Maximum upward deflection 
P     = 49.2        #Load in actuator 2 
i     = 1


stiff = []


h     = h/100                                  #h in m 
tsk   = tsk/1000         #in m 
tsp   = tsp/1000        #in m 
tst   = tst/1000         #in m 
hst   = hst/100        #in m 
wst   = wst/100        #in m 

r   = h/2                                  #Radius of the semi circular part of the airfoil in m 
sk  = np.sqrt(r**2 + (Ca-r)**2)      #Length of straight part of the airfoil in m 
phi = atan((r)/(Ca-r))                 #Angle of the non circular part of the airfoil 

At = (pi*r**2)/2+ r*(Ca-r)
print(At)


A_sk   = sk*tsk                          #Area of the skin straight panel in m2
A_cir  = pi*r*tsk                       #Area of the skin circular panel in m2
A_sp   = (h)*tsp                             #Area of the spar in m2
Ast    =  (wst*tst + (hst-tst)*tst)                 #Area of a stiffener m2

ycirc = 0 
zcirc = -(r-2*r/(pi))  

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
 

plt.plot(graph[:,1],graph[:,2],'o')
plt.axis('equal')
plt.show()

stiff.append(Skt)
stiff.append(Skb)
stiff.append(Ssp)
stiff.append(Sci)


stiff = np.asarray(stiff)

print(stiff[:,1],stiff[:,2])

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


Ispary = (h*(tsp)**3)/12
Isparz = (tsp*h**3)/12 

Icirc  = pi*r**3*tsk


Iskz   = (sk**3*tsk*sin(phi)**2)/12
Isky   = (sk**3*tsk*cos(phi)**2)/12

Izz = Isparz + Icirc + 2*Iskz + sum(stiff[:,10]) + A_cir*(4*r/(3*pi))**2
Iyy = Ispary + Icirc + 2*Isky + sum(stiff[:,9]) 

print(zbar)
print(ybar)
print("Iyy is", Iyy)
print("Izz is", Izz)
#print(stiff)