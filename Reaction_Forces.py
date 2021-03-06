import numpy as np
from NEW_Forces_Deflections import locationvalue,output
from Interpolation import patchinterpolate


def reaction_forces(la,x1,x2,x3,xa,h,d1,d3,theta,P,E,G,zsc,Iyy,Izz,J,xlst,Vlst,Mlst,defllst,Tlst,thetalst):

    # Transformation of some input parameters

    # Determine the distance to the actuators
    xf = x2-xa/2
    xp = x2+xa/2

    # Determine the total aerodynamic load
    Vtot = Vlst[-1]

    # Decompose d1 and d3 along the y and z axis
    d1y = d1*np.cos(theta*np.pi/180)  # [m]
    d1z = -d1*np.sin(theta*np.pi/180)  # [m]

    d3y = d3*np.cos(theta*np.pi/180)  # [m]
    d3z = -d3*np.sin(theta*np.pi/180)  # [m]

    # Decompose force P along the y and z axis
    Py = P*np.sin(theta*np.pi/180)
    Pz = P*np.cos(theta*np.pi/180)

    # Determine Components of the F force
    CFy = np.sin(theta*np.pi/180)
    CFz = np.cos(theta*np.pi/180)

    # Calculate constants used in the displacement functions
    EIzz = E*Izz
    EIyy = E*Iyy
    GJ = G*J

    '''
    To be able to solve for the unknown reaction forces, a matrix equation will be set up
    
    Ab = c
    where:
        A = The moment and displacement equations
        b = [Ay,Az,By,Bz,Cy,Cz,F ,C1,C2,C3,C4,C5].T
        b = [0 ,1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,11].T
        c = The answers of the moment and displacement equations
    '''

    # Locations of the unknowns
    iAy = 0
    iAz = 1
    iBy = 2
    iBz = 3
    iCy = 4
    iCz = 5
    iF  = 6
    iC1 = 7
    iC2 = 8
    iC3 = 9
    iC4 = 10
    iC5 = 11


    # Setting up functions to calculate the moments and displacements for location x

    # Moment around the y axis
    def My(x):
        # Create the row for the A matrix
        row = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

        # Create the sum with parameters independent of the unknowns
        additional_sum = 0

        # Apply the Macaulay step function
        if x-x1 >= 0:
            row[iAz] = -1*(x-x1)
        if x-xf >= 0:
            row[iF] = -CFz*(x-xf)
        if x - x2 >= 0:
            row[iBz] = -1*(x-x2)
        if x - x3 >= 0:
            row[iCz] = -1*(x-x3)

        if x - xp >= 0:
            additional_sum += Pz*(x-xp)

        return row,additional_sum


    # Moment around the z axis
    def Mz(x):
        row = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
        additional_sum = (-locationvalue(xlst,x,Mlst))  # Mlst is negative, needs to be +
        if x-x1>=0:
            row[iAy] = -1*(x-x1)
        if x-xf>=0:
            row[iF] = -CFy*(x-xf)
        if x-x2>= 0:
            row[iBy] = -1*(x-x2)
        if x-x3>= 0:
            row[iCy] = -1*(x-x3)
        if x-xp>= 0:
            additional_sum += Py*(x-xp)
        return row,additional_sum


    # Torque
    def T(x):
        row = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
        additional_sum = locationvalue(xlst,x,Tlst)
        if x-x1>=0:
            row[iAy] = -(zsc+h/2)
        if x-xf>=0:
            row[iF] = CFy*(-zsc)-CFz*h/2
        if x-x2>= 0:
            row[iBy] = -(zsc+h/2)
        if x-x3>= 0:
            row[iCy] = -(zsc+h/2)
        if x-xp>=0:
            additional_sum += -Py*(-zsc)+Pz*h/2
        return row,additional_sum


    # Displacement in y direction
    def v(x):
        row = np.array([0., 0., 0., 0., 0., 0., 0., x, 1., 0., 0., 0.])
        additional_sum = locationvalue(xlst,x,defllst)/EIzz
        if x-x1>=0:
            row[iAy] = 1/(6*EIzz)*(x-x1)**3
        if x-xf>=0:
            row[iF] = CFy/(6*EIzz)*(x-xf)**3
        if x - x2 >= 0:
            row[iBy] = 1/(6*EIzz)*(x-x2)**3
        if x - x3 >= 0:
            row[iCy] = 1/(6*EIzz)*(x-x3)**3
        if x - xp >= 0:
            additional_sum += -Py/(6*EIzz)*(x-xp)**3
        return row,additional_sum


    # Displacement in z direction
    def w(x):
        row = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., x, 1., 0.])
        additional_sum = 0
        if x-x1>=0:
            row[iAz] = 1/(6*EIyy)*(x-x1)**3
        if x-xf>=0:
            row[iF] = CFz/(6*EIyy)*(x-xf)**3
        if x-x2>=0:
            row[iBz] = 1/(6*EIyy)*(x-x2)**3
        if x-x3>=0:
            row[iCz] = 1/(6*EIyy)*(x-x3)**3
        if x-xp>=0:
            additional_sum += -Pz/(6*EIyy)*(x-xp)**3
        return row,additional_sum


    # Twist
    def ftheta(x):
        row = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.])
        additional_sum = locationvalue(xlst,x,thetalst)/GJ
        if x-x1>=0:
            row[iAy] = 1/GJ*(-(zsc+h/2))*(x-x1)
        if x-xf>=0:
            row[iF] = CFy/GJ*(-zsc)*(x-xf)-CFz/GJ*h/2*(x-xf)
        if x-x2>=0:
            row[iBy] = 1/GJ*(-(zsc+h/2))*(x-x2)
        if x-x3>=0:
            row[iCy] = 1/GJ*(-(zsc+h/2))*(x-x3)
        if x-xp>=0:
            additional_sum += -Py/GJ*(-zsc)*(x-xp)+Pz/GJ*h/2*(x-xp)
        return row,additional_sum


    # Sum of forces in z direction
    def Sz(x):
        row = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
        additional_sum = 0
        if x-x1>=0:
            row[iAz] = 1
        if x-xf>=0:
            row[iF] = CFz
        if x-x2>=0:
            row[iBz] = 1
        if x-x3>=0:
            row[iCz] = 1
        if x-xp>=0:
            additional_sum += -Pz
        return row,additional_sum


    # Sum of forces in y direction
    def Sy(x):
        row = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
        additional_sum = Vtot
        if x-x1>=0:
            row[iAy] = 1
        if x-xf>=0:
            row[iF] = CFy
        if x-x2>=0:
            row[iBy] = 1
        if x-x3>=0:
            row[iCy] = 1
        if x-xp>=0:
            additional_sum += -Py
        return row,additional_sum


    # Solving for the unknowns

    # Creating the A matrix with output rows from the functions
    A = np.array([My(la)[0],
                  Mz(la)[0],
                  T(la)[0],
                  Sz(la)[0],
                  Sy(la)[0],
                  v(x1)[0] + ftheta(x1)[0]*zsc,
                  v(x2)[0] + ftheta(x2)[0]*zsc,
                  v(x3)[0] + ftheta(x3)[0]*zsc,
                  w(x1)[0],
                  (w(xf)[0] + ftheta(xf)[0]*h/2)/np.cos(theta*np.pi/180),
                  w(x2)[0],
                  w(x3)[0]])

    # Creating the c vector and subtracting the additional sum from the boundary condition
    c = np.array([[-My(la)[1]],
                  [-Mz(la)[1]],
                  [-T(la)[1]],
                  [-Sz(la)[1]],
                  [-Sy(la)[1]],
                  [d1y - (v(x1)[1] - ftheta(x1)[1]*zsc)],
                  [-(v(x2)[1] - ftheta(x2)[1]*zsc)],
                  [d3y - (v(x3)[1] - ftheta(x3)[1]*zsc)],
                  [d1z - w(x1)[1]],
                  [-(w(xf)[1] + ftheta(xf)[1]*h/2)/np.cos(theta*np.pi/180)],
                  [-w(x2)[1]],
                  [d3z - w(x3)[1]]])

    # Solving for the unknowns in vector b
    b = np.linalg.solve(A,c)

    Ay = b[iAy][0]
    Az = b[iAz][0]
    By = b[iBy][0]
    Bz = b[iBz][0]
    Cy = b[iCy][0]
    Cz = b[iCz][0]
    F  = b[iF][0]
    C1 = b[iC1][0]
    C2 = b[iC2][0]
    C3 = b[iC3][0]
    C4 = b[iC4][0]
    C5 = b[iC5][0]

    Fy = F*CFy
    Fz = F*CFz

    return Ay,Az,By,Bz,Cy,Cz,Fy,Fz,C1,C2,C3,C4,C5

'''
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
zsc   = -0.085      # [m]
P     = 49.2E3      # [N]

data = np.loadtxt('AERO.dat',delimiter = ',')
xlst, zlst, qlst = patchinterpolate(600,600,data)
Vlst,Mlst,defllst,Tlst,thetalst = output(xlst,zlst,qlst,zsc)

Ay,Az,By,Bz,Cy,Cz,Fy,Fz,C1,C2,C3,C4,C5,b = reaction_forces(la,x1,x2,x3,xa,h,d1,d3,theta,P,E,G,zsc,Iyy,Izz,J,xlst,Vlst,Mlst,defllst,Tlst,thetalst)
'''