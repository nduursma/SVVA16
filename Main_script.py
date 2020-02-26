import numpy as np
from Interpolation import patchinterpolate
from NEW_Forces_Deflections import output, locationvalue
from Reaction_Forces import reaction_forces
from Cross_sectional_properties import cross_sect_prop
from ShearCalcs import CalcTorsStif
import matplotlib.pyplot as plt
from Deflections import deflections


Ca    = 0.505       #Chord Length Aileron 
la    = 1.611       #Span of the Aileron
x1    = 0.125       #x-location of hinge 1
x2    = 0.498       #x-location of hinge 2 
x3    = 1.494       #x-location of hinge 3
xa    = 24.5E-2        #Distance between actuator 1 and 2
h     = 16.1E-2        #Aileron Height 
tsk   = 1.1E-3         #Skin thickness
tsp   = 2.4E-3         #Spar Thickness
tst   = 1.2E-3         #Thickness of stiffener
hst   = 1.3E-2         #Height of stifffener
wst   = 1.7E-2         #Width of stiffeners 
nst   = 11          #Number of stiffeners (equalt spaced)
d1    = 0.389E-2       #Vertical displacement of hinge 1
d3    = 1.245E-2       #Vertical displacement of hinge 2
theta = 30          #Maximum upward deflection 
P     = 49.2E3        #Load in actuator 2 
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

x = xlst #np.arange(0,la+dx,dx)
y = []
z = []
twist = []
Sy = []
Sz = []
momenty = []
momentz = []
torque = []

for xi in x:
    Myx,Mzx,Tx,Syx,Szx,vx,wx,thetax = deflections(xi,x1,x2,x3,xa,h,theta,P,E,G,zsc,Iyy,Izz,J,xlst,Vlst,Mlst,defllst,Tlst,thetalst,Ay,Az,By,Bz,Cy,Cz,Fy,Fz,C1,C2,C3,C4,C5)
    y.append(vx+thetax*zsc)
    z.append(wx)
    twist.append(thetax)
    momenty.append(Myx)
    momentz.append(Mzx)
    torque.append(Tx)
    Sy.append(Syx)
    Sz.append(Szx)
momentz = np.array(momentz)


plt.plot(x,y)
plt.title('Vertical Displacement')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.grid()
plt.savefig('Vertical_Displacement.jpg')
plt.show()

plt.plot(x,z)
plt.title('Horizontal Displacement')
plt.xlabel('x [m]')
plt.ylabel('z [m]')
plt.grid()
plt.savefig('Horizontal_Displacement.jpg')
plt.show()

plt.plot(x,twist)
plt.title('Twist')
plt.xlabel('x [m]')
plt.ylabel('Theta [rad]')
plt.grid()
plt.savefig('Twist.jpg')
plt.show()

plt.plot(x,momenty)
plt.title('Moment around the y axis')
plt.xlabel('x [m]')
plt.ylabel('My [Nm]')
plt.grid()
plt.savefig('Moment_y.jpg')
plt.show()

plt.plot(x,momentz)
plt.title('Moment around the z axis')
plt.xlabel('x [m]')
plt.ylabel('Mz [Nm]')
plt.grid()
plt.savefig('Moment_z.jpg')
plt.show()

plt.plot(x,torque)
plt.title('Torque')
plt.xlabel('x [m]')
plt.ylabel('T [Nm]')
plt.grid()
plt.savefig('Torque.jpg')
plt.show()
