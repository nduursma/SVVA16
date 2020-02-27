import numpy as np
from NEW_Forces_Deflections import locationvalue

def deflections(x,x1,x2,x3,xa,h,theta,P,E,G,zsc,Iyy,Izz,J,xlst,Vlst,Mlst,defllst,Tlst,thetalst,Ay,Az,By,Bz,Cy,Cz,Fy,Fz,C1,C2,C3,C4,C5):

    # Transformation of some input parameters

    # Determine the total aerodynamic load
    Vtot = Vlst[-1]

    # Decompose force P along the y and z axis
    Py = P*np.sin(theta*np.pi/180)
    Pz = P*np.cos(theta*np.pi/180)

    # Determine the distance to the actuators
    xf = x2-xa/2
    xp = x2+xa/2

    # Calculate constants used in the displacement functions
    EIzz = E*Izz
    EIyy = E*Iyy
    GJ = G*J


    # Moment around the y axis
    def My(x):

        M = 0

        # Apply the Macaulay step function
        if x-x1>=0:
            M += -Az*(x-x1)
        if x-xf>=0:
            M += -Fz*(x-xf)
        if x-x2>=0:
            M += -Bz*(x-x2)
        if x-x3>=0:
            M += -Cz*(x-x3)

        if x-xp>=0:
            M += Pz*(x-xp)

        return M

    # Moment around the z axis
    def Mz(x):
        M = (-locationvalue(xlst, x, Mlst))  # Mlst is negative, needs to be +
        if x-x1>=0:
            M += -Ay*(x-x1)
        if x-xf>=0:
            M += -Fy*(x-xf)
        if x-x2>=0:
            M += -By*(x-x2)
        if x-x3>=0:
            M += -Cy*(x-x3)
        if x-xp>=0:
            M += Py*(x-xp)
        return M

    # Torque
    def T(x):
        Torque = locationvalue(xlst,x,Tlst)
        if x-x1>=0:
            Torque += Ay*(zsc-h/2)
        if x-xf>=0:
            Torque += Fy*(-zsc) - Fz*(h/2)
        if x-x2>=0:
            Torque += By*(zsc-h/2)
        if x-x3>=0:
            Torque += Cy*(zsc-h/2)
        if x-xp>=0:
            Torque += -Py*(-zsc)+Pz*h/2
        return Torque

    # Displacement in y direction
    def v(x):
        defl = locationvalue(xlst,x,defllst)/EIzz + C1*x + C2
        if x-x1>=0:
            defl += Ay/(6*EIzz)*(x-x1)**3
        if x-xf>=0:
            defl += Fy/(6*EIzz)*(x-xf)**3
        if x-x2>=0:
            defl += By/(6*EIzz)*(x-x2)**3
        if x-x3>=0:
            defl += Cy/(6*EIzz)*(x-x3)**3
        if x-xp>=0:
            defl += -Py/(6*EIzz)*(x-xp)**3
        return defl

    # Displacement in z direction
    def w(x):
        defl = C3*x + C4
        if x-x1>=0:
            defl += Az/(6*EIyy)*(x-x1)**3
        if x-xf>=0:
            defl += Fz/(6*EIyy)*(x-xf)**3
        if x-x2>=0:
            defl += Bz/(6*EIyy)*(x-x2)**3
        if x-x3>=0:
            defl += Cz/(6*EIyy)*(x-x3)**3
        if x-xp>=0:
            defl += -Py/(6*EIyy)*(x-xp)**3
        return defl

    # Twist
    def ftheta(x):
        twist = locationvalue(xlst,x,thetalst)/GJ + C5
        if x-x1>=0:
            twist += Ay/GJ*(zsc-h/2)*(x-x1)
        if x-xf>=0:
            twist += Fy/GJ*(-zsc)*(x-xf) - Fz/GJ*h/2*(x-xf)
        if x-x2>=0:
            twist += By/GJ*(zsc-h/2)*(x-x2)
        if x-x3>=0:
            twist += Cy/GJ*(zsc-h/2)*(x-x3)
        if x-xp>=0:
            twist += -Py/GJ*(-zsc)*(x-xp)+Pz/GJ*h/2*(x-xp)
        return twist

    # Sum of forces in z direction
    def Sz(x):
        Sumz = 0
        if x-x1>=0:
            Sumz += Az
        if x-xf>=0:
            Sumz += Fz
        if x-x2>=0:
            Sumz += Bz
        if x-x3>=0:
            Sumz += Cz
        if x-xp>=0:
            Sumz += -Pz
        return Sumz

    # Sum of forces in y direction
    def Sy(x):
        Sumy = Vtot
        if x-x1>=0:
            Sumy += Ay
        if x-xf>=0:
            Sumy += Fy
        if x-x2>=0:
            Sumy += By
        if x-x3>=0:
            Sumy += Cy
        if x-xp>=0:
            Sumy += -Py
        return Sumy


    Myx = My(x)
    Mzx = Mz(x)
    Tx = T(x)
    Syx = Sy(x)
    Szx = Sz(x)
    vx = v(x)
    wx = w(x)
    thetax = ftheta(x)

    return Myx,Mzx,Tx,Syx,Szx,vx,wx,thetax


'''
from Interpolation import patchinterpolate
from NEW_Forces_Deflections import output
from Reaction_Forces import reaction_forces
import matplotlib.pyplot as plt

ca    = 0.505       # [m]
la    = 1.611       # [m]
h     = 16.1E-2     # [m]
x1    = 0.125       # [m]
x2    = 0.498       # [m]
x3    = 1.494       # [m]
xa    = 24.5E-2     # [m]
d1    = 0.389E-2    # [m]
d3    = 1.245E-2    # [m]
theta = 30          # [deg]
E     = 73.1E9      # [Pa]
G     = 28E9        # [Pa]
Iyy   = 4.59E-5     # [m4]
Izz   = 4.75E-6     # [m4]
J     = 7.749E-6    # [m4]
zsc   = -0.0855     # [m]
P     = 49.2E3      # [N]

data = np.loadtxt('AERO.dat',delimiter = ',')
xlst, zlst, qlst = patchinterpolate(600,600,data)
Vlst,Mlst,defllst,Tlst,thetalst = output(xlst,zlst,qlst,zsc)

Ay,Az,By,Bz,Cy,Cz,Fy,Fz,C1,C2,C3,C4,C5 = reaction_forces(la,x1,x2,x3,xa,h,d1,d3,theta,P,E,G,zsc,Iyy,Izz,J,xlst,Vlst,Mlst,defllst,Tlst,thetalst)

dx = 0.001
x = np.arange(0,la+dx,dx)
y = []
z = []
twist = []
momenty = []
momentz = []
torque = []
sheary = []

for xi in x:
    Myx,Mzx,Tx,Syx,Szx,vx,wx,thetax = deflections(xi,x1,x2,x3,xa,h,theta,P,E,G,zsc,Iyy,Izz,J,xlst,Vlst,Mlst,defllst,Tlst,thetalst,Ay,Az,By,Bz,Cy,Cz,Fy,Fz,C1,C2,C3,C4,C5)
    y.append(vx+thetax*zsc)
    z.append(wx)
    twist.append(thetax)
    momenty.append(Myx)
    momentz.append(Mzx)
    torque.append(Tx)
    sheary.append(Syx)

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
'''