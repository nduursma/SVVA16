# Imports
import numpy as np
from Interpolation import patchinterpolate
from NEW_Forces_Deflections import output, locationvalue
from Reaction_Forces import reaction_forces
from Cross_sectional_properties import cross_sect_prop
from ShearCalcs import CalcTorsStif
import matplotlib.pyplot as plt
from Deflections import deflections
from AllStressCalcs import CalcShearCenter,AileronFullStress


# Input parameters
Ca    = 0.505       #Chord Length Aileron 
la    = 1.611       #Span of the Aileron
x1    = 0.125       #x-location of hinge 1
x2    = 0.498       #x-location of hinge 2 
x3    = 1.494       #x-location of hinge 3
xa    = 24.5E-2     #Distance between actuator 1 and 2
h     = 16.1E-2     #Aileron Height
tsk   = 1.1E-3      #Skin thickness
tsp   = 2.4E-3      #Spar Thickness
tst   = 1.2E-3      #Thickness of stiffener
hst   = 1.3E-2      #Height of stifffener
wst   = 1.7E-2      #Width of stiffeners
nst   = 11          #Number of stiffeners (equalt spaced)
d1    = 0.389E-2    #Vertical displacement of hinge 1
d3    = 1.245E-2    #Vertical displacement of hinge 2
theta = 30          #Maximum upward deflection 
P     = 49.2E3      #Load in actuator 2
E     = 73.1E9      # [Pa]
G     = 28E9        # [Pa]


# Calculated parameters
J     = CalcTorsStif(Ca,h,tsk,tsp)
Dx, zbar,ybar,Izz,Iyy, stiff = cross_sect_prop(h,tsk,tsp,tst,hst,wst,Ca,nst)
#zsc = -0.085       # [m]
zsc = CalcShearCenter(h,Ca,tsk,tsp,Iyy,Izz,stiff,Dx,zbar)


# Reading the aerodynamic load
data = np.loadtxt('AERO.dat',delimiter = ',')


# Transforming the coordinates
xl = []
zl = []
Nz = len(data)
Nx = len(data[0])

for i in range(1,Nz+1):  
    thetazi  = (i-1)/Nz*np.pi
    thetazi1 = (i)/Nz*np.pi
    z  =  -(1/2)*((Ca/2)*(1-np.cos(thetazi))+(Ca/2)*(1-np.cos(thetazi1)))
    zl.append(z)

for k in range(1,Nx+1):
    thetaxk  = (k-1)/Nx*np.pi
    thetaxk1 = (k)/Nx*np.pi
    x  = (1/2)*((la/2)*(1-np.cos(thetaxk))+(la/2)*(1-np.cos(thetaxk1)))   
    xl.append(x)


# Interpolate the aerodynamic load
xlst, zlst, qlst = patchinterpolate(600,600,xl,zl,data)


# Integrate the aerodynamic load
Vlst, Mlst, d1lst, defllst, taulst, Tlst, thetalst = output(xlst, zlst, qlst, zsc)

# Calculate the reaction forces
Ay,Az,By,Bz,Cy,Cz,Fy,Fz,C1,C2,C3,C4,C5 = reaction_forces(la,x1,x2,x3,xa,h,d1,d3,theta,P,E,G,zsc,Iyy,
                                                         Izz, J,xlst,Vlst,Mlst,defllst,Tlst,thetalst)


# Calculate displacements, moments, etc. for every x location
x = xlst #np.arange(0,la+dx,dx)

y = []
slopey = []
momentz = []
sheary = []

z = []
slopez = []
momenty = []
shearz = []

twist = []
torque = []

for xi in x:
    Myx,Mzx,Tx,Syx,Szx,vsx,vx,wsx,wx,thetax = deflections(xi,x1,x2,x3,xa,h,theta,P,E,G,zsc,Iyy,Izz,J,xlst,Vlst,Mlst,d1lst,defllst,Tlst,thetalst,Ay,Az,By,Bz,Cy,Cz,Fy,Fz,C1,C2,C3,C4,C5)
    y.append(vx)
    slopey.append(vsx)
    momentz.append(Mzx)
    sheary.append(Syx)

    z.append(wx)
    slopez.append(wsx)
    momenty.append(Myx)
    shearz.append(Szx)

    twist.append(thetax)
    torque.append(Tx)


# Calculate the maximum stress and its location
maximumstress,zloc,yloc,xloc = AileronFullStress(x,sheary,shearz,momenty,momentz,torque,Ca,h,la,tsp,tsk,Iyy,Izz,zbar,Dx,stiff)


# Creating the plots
# y
plt.subplot(221)
plt.plot(x,y)
plt.title('Deflection in y direction')
plt.xlabel('x [m]')
plt.ylabel('v(x) [m]')
plt.grid()

plt.subplot(222)
plt.plot(x,slopey)
plt.title('Slope in y direction')
plt.xlabel('x [m]')
plt.ylabel("v'(x) [/m]")
plt.grid()

plt.subplot(223)
plt.plot(x,momentz)
plt.title('Moment around the z axis')
plt.xlabel('x [m]')
plt.ylabel('Mz(x) [Nm]')
plt.grid()

plt.subplot(224)
plt.plot(x,sheary)
plt.title('Shear force in y direction')
plt.xlabel('x [m]')
plt.ylabel('Sy(x) [N]')
plt.grid()

plt.show()

# z
plt.subplot(221)
plt.plot(x,z)
plt.title('Deflection in z direction')
plt.xlabel('x [m]')
plt.ylabel('w(x) [m]')
plt.grid()

plt.subplot(222)
plt.plot(x,slopez)
plt.title('Slope in z direction')
plt.xlabel('x [m]')
plt.ylabel("w'(x) [/m]")
plt.grid()

plt.subplot(223)
plt.plot(x,momenty)
plt.title('Moment around the y axis')
plt.xlabel('x [m]')
plt.ylabel('My(x) [Nm]')
plt.grid()

plt.subplot(224)
plt.plot(x,shearz)
plt.title('Shear force in z direction')
plt.xlabel('x [m]')
plt.ylabel('Sz(x) [N]')
plt.grid()

plt.show()

# x
plt.subplot(221)
plt.plot(x,twist)
plt.title('Twist')
plt.xlabel('x [m]')
plt.ylabel('theta(x) [rad]')
plt.grid()

plt.subplot(222)
plt.plot(x,torque)
plt.title('Torque')
plt.xlabel('x [m]')
plt.ylabel("T(x) [/m]")
plt.grid()

plt.subplot(223)
plt.plot(x,taulst)
plt.title('Distributed torque')
plt.xlabel('x [m]')
plt.ylabel('tau(x) [Nm]')
plt.grid()

plt.show()