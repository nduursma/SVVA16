import numpy as np
from Interpolation import patchinterpolate
from NEW_Forces_Deflections import output, locationvalue
from Reaction_Forces import reaction_forces
from Cross_sectional_properties import cross_sect_prop
from ShearCalcs import CalcTorsStif


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
E     = 73.1E9      # [Pa]
G     = 28E9        # [Pa]
J     = CalcTorsStif(Ca,h,tsk,tsp)


data = np.loadtxt('AERO.dat',delimiter = ',')

zsc = -0.085
xlst, zlst, qlst = patchinterpolate(600,600,data)
Vlst, Mlst, defllst, Tlst, thetalst = output(xlst, zlst, qlst, zsc)
Dx, zbar,ybar,Izz,Iyy, stiff = cross_sect_prop(h,tsk,tsp,tst,hst,wst,Ca,nst)
Ay,Az,By,Bz,Cy,Cz,Fy,Fz,C1,C2,C3,C4,C5 = reaction_forces(la,x1,x2,x3,xa,h,d1,d3,theta,P,E,G,zsc,Iyy,
                                                         Izz, J,xlst,Vlst,Mlst,defllst,Tlst,thetalst)
