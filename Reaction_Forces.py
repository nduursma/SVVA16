'''
Calculation of the reaction forces
'''

# Imports
import numpy as np
from Resultant_Centroid_DistributedLoadCorrect import Moment_z,Torque_zsc,Magnitude_Centroid


# Necessary Inputs
ca    = 0.505       # [m]
la    = 1.611       # [m]
h     = 16.1        # [cm]
x1    = 0.125       # [m]
x2    = 0.498       # [m]
x3    = 1.494       # [m]
xa    = 24.5/100    # [m]
d1    = 0.389       # [m]
d3    = 1.245       # [m]
theta = 30          # [deg]
E     = 20          # [GPa]
G     = 20          # [GPa]
Iyy   = 20          # [m4]
Izz   = 20          # [m4]
J     = 20          # [m4]
zsc   = -0.2        # [m]
P     = 49.2        # [kN]
qtot, Maglst, Centroid_Zlst  = Magnitude_Centroid()

# Transform of some input parameters

# Transform h from [cm] to [m]
h = h/100 # [m]

# Decompose d1 and d3 along the y and z axis
d1y = d1*np.cos(theta*np.pi/180) # [m]
d1z = d1*np.sin(theta*np.pi/180) # [m]

d3y = d3*np.cos(theta*np.pi/180) # [m]
d3z = d3*np.sin(theta*np.pi/180) # [m]

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


'''
To be able to solve for the unknown reaction forces, a matrix equation will be set up

Ab = c
where:
    A = The moment and displacement equations
    b = [Ay,Az,By,Bz,Cy,Cz,Fy,Fz,C1,C2,C3,C4,C5].T
    b = [0 ,1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,11,12].T
    c = The answers of the moment and displacement equations
'''


# Locations of the unknowns
Ay = 0
Az = 1
By = 2
Bz = 3
Cy = 4
Cz = 5
Fy = 6
Fz = 7
C1 = 8
C2 = 9
C3 = 10
C4 = 11
C5 = 12


# Setting up functions to calculate the moments and displacements for location x

# Moment around the y axis
def My(x):
    # Create the row for the A matrix
    row = np.array([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.])

    # Create the sum with parameters independent of the unknowns
    additional_sum = 0

    # Apply the Macaulay step function
    if x-x1>=0:
        row[Az] = -1*(x-x1)
    if x-xf>=0:
        row[Fz] = -1*(x-xf)
    if x-x2>=0:
        row[Bz] = -1*(x-x2)
    if x-x3>=0:
        row[Cz] = -1*(x-x3)

    if x-xp>=0:
        additional_sum = Pz*(x-xp)

    return row,additional_sum

# Moment around the z axis
def Mz(x):
    row = np.array([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.])
    additional_sum = -Moment_z(x)
    if x-x1>=0:
        row[Ay] = -1*(x-x1)
    if x-xf>=0:
        row[Fy] = -1*(x-xf)
    if x-x2>=0:
        row[By] = -1*(x-x2)
    if x-x3>=0:
        row[Cy] = -1*(x-x3)
    if x-xp>=0:
        additional_sum = Pz*(x-xp)
    return row,additional_sum

# Torque
def T(x):
    row = np.array([1/GJ*(0-zsc),0.,1/GJ*(0-zsc),0.,1/GJ*(0-zsc),0.,1/GJ*(0-zsc),0.,0.,0.,0.,0.,0.])
    additional_sum = Py/GJ*(0-zsc)+Torque_zsc(x,zsc)
    return row, additional_sum

# Displacement in y direction
def v(x):
    row = np.array([0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,0.,0.,0.])
    additional_sum = 0
    if x-x1>=0:
        row[Ay] = 1/(6*EIzz)*(x-x1)**3
    if x-xf>=0:
        row[Fy] = 1/(6*EIzz)*(x-xf)**3
    if x-x2>=0:
        row[By] = 1/(6*EIzz)*(x-x2)**3
    if x-x3>=0:
        row[Cy] = 1/(6*EIzz)*(x-x3)**3
    if x-xp>=0:
        additional_sum = -Py/(6*EIzz)*(x-xp)**3
    return row,additional_sum

# Displacement in z direction
def w(x):
    row = np.array([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,0.])
    additional_sum = 0
    if x-x1>=0:
        row[Az] = 1/(6*EIzz)*(x-x1)**3
    if x-xf>=0:
        row[Fz] = 1/(6*EIzz)*(x-xf)**3
    if x-x2>=0:
        row[Bz] = 1/(6*EIzz)*(x-x2)**3
    if x-x3>=0:
        row[Cz] = 1/(6*EIzz)*(x-x3)**3
    if x-xp>=0:
        additional_sum = -Py/(6*EIzz)*(x-xp)**3
    return row,additional_sum

# Twist
def ftheta(x):
    row = np.array([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.])
    additional_sum = 0
    if x-x1>=0:
        row[Ay] = 1/GJ*(0-zsc)*(x-x1)
    if x-xf>=0:
        row[Fy] = 1/GJ*(0-zsc)*(x-xf)
    if x-x2>=0:
        row[By] = 1/GJ*(0-zsc)*(x-x2)
    if x-x3>=0:
        row[Cy] = 1/GJ*(0-zsc)*(x-x3)
    if x-xp>=0:
        additional_sum = -Py/GJ*(0-zsc)*(x-xp)
    return row,additional_sum

# Sum of forces in z direction
def Sz():
    row = np.array([0.,-1.,0.,-1.,0.,-1.,0.,-1.,0.,0.,0.,0.,0.])
    additional_sum = Pz
    return row,additional_sum

# Sum of forces in y direction
def Sy():
    row = np.array([-1.,0.,-1.,0.,-1.,0.,-1.,0.,0.,0.,0.,0.,0.])
    additional_sum = Py + qtot
    return row,additional_sum


# Creating the A matrix with output rows from the functions
A = np.array([My(la)[0],
              Mz(la)[0],
              T(la)[0],
              Sz()[0],
              Sy()[0],
              v(x1)[0]+ftheta(x1)[0]*zsc,
              v(xf)[0]+ftheta(xf)[0]*zsc,
              v(x2)[0]+ftheta(x2)[0]*zsc,
              v(x3)[0]+ftheta(x3)[0]*zsc,
              w(x1)[0],
              w(xf)[0],
              w(x2)[0],
              w(x3)[0]])

# Creating the c vector and subtracting the additional sum from the boundary condition
c = np.array([[-My(la)[1]],
              [-Mz(la)[1]],
              [-T(la)[1]],
              [-Sz()[1]],
              [-Sy()[1]],
              [d1y-v(x1)[1]-ftheta(x1)[1]*zsc],
              [0-v(xf)[1]-ftheta(xf)[1]*zsc],
              [0-v(x2)[1]-ftheta(x2)[1]*zsc],
              [d3y-v(x3)[1]-ftheta(x3)[1]*zsc],
              [d1z-w(x1)[1]],
              [0-w(xf)[1]],
              [0-w(x2)[1]],
              [d3z-w(x3)[1]]])

# Solving for the unknowns in vector b
#b = np.linalg.solve(A,c)

#print(b)