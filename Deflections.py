import numpy as np
from Reaction_Forces import reaction_forces
from NEW_Forces_Deflections import output,locationvalue
import matplotlib.pyplot as plt

ca    = 0.505       # [m]
la    = 1.611       # [m]
h     = 16.1E-2     # [m]
x1    = 0.125       # [m]
x2    = 0.498       # [m]
x3    = 1.494       # [m]
xa    = 24.5/100    # [m]
theta = 30          # [deg]
E     = 73.1E9      # [Pa]
G     = 28E9        # [Pa]
Iyy   = 4.59E-5     # [m4]
Izz   = 4.75E-6     # [m4]
J     = 7.749E-06   # [m4]
zsc   = -0.085      # [m]
P     = 49.2E3      # [N]

rf = reaction_forces()
Rqlst,Mzqlst,vqlst,Tqlst,fthetaqlst = output(zsc)

Py = P*np.sin(theta*np.pi/180)
Pz = P*np.cos(theta*np.pi/180)

xf = x2-xa/2
xp = x2+xa/2

EIzz = E*Izz
EIyy = E*Iyy
GJ = G*J


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


# Moment around the y axis
def My(x):
    # Create the row for the A matrix
    row = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

    # Create the sum with parameters independent of the unknowns
    additional_sum = 0

    # Apply the Macaulay step function
    if x - x1 >= 0:
        row[Az] = -1 * (x - x1)
    if x - xf >= 0:
        row[Fz] = -1 * (x - xf)
    if x - x2 >= 0:
        row[Bz] = -1 * (x - x2)
    if x - x3 >= 0:
        row[Cz] = -1 * (x - x3)

    if x - xp >= 0:
        additional_sum += Pz * (x - xp)

    M = row.dot(rf)+additional_sum

    return M[0]


# Moment around the z axis
def Mz(x):
    row = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
    additional_sum = (-locationvalue(x, Mzqlst))

    if x - x1 >= 0:
        row[Ay] = -1 * (x - x1)
    if x - xf >= 0:
        row[Fy] = -1 * (x - xf)
    if x - x2 >= 0:
        row[By] = -1 * (x - x2)
    if x - x3 >= 0:
        row[Cy] = -1 * (x - x3)
    if x - xp >= 0:
        additional_sum += Pz * (x - xp)

    M = row.dot(rf)+additional_sum

    return M[0]


# Torque
def T(x):
    row = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
    additional_sum = locationvalue(x, Tqlst)

    if x - x1 >= 0:
        row[Ay] = -zsc
    if x - xf >= 0:
        row[Fy] = -zsc
    if x - x2 >= 0:
        row[By] = -zsc
    if x - x3 >= 0:
        row[Cy] = -zsc
    if x - xp >= 0:
        additional_sum += -Py * (-zsc)

    T = row.dot(rf) + additional_sum

    return T[0]


# Displacement in y direction
def v(x):
    row = np.array([0., 0., 0., 0., 0., 0., 0., 0., x, 1., 0., 0., 0.])
    additional_sum = locationvalue(x, vqlst) / EIzz

    if x - x1 >= 0:
        row[Ay] = 1 / (6 * EIzz) * (x - x1) ** 3
    if x - xf >= 0:
        row[Fy] = 1 / (6 * EIzz) * (x - xf) ** 3
    if x - x2 >= 0:
        row[By] = 1 / (6 * EIzz) * (x - x2) ** 3
    if x - x3 >= 0:
        row[Cy] = 1 / (6 * EIzz) * (x - x3) ** 3
    if x - xp >= 0:
        additional_sum += -Py / (6 * EIzz) * (x - xp) ** 3

    defl = row.dot(rf) + additional_sum

    return defl[0]


# Displacement in z direction
def w(x):
    row = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., x, 1., 0.])
    additional_sum = 0

    if x - x1 >= 0:
        row[Az] = 1 / (6 * EIyy) * (x - x1) ** 3
    if x - xf >= 0:
        row[Fz] = 1 / (6 * EIyy) * (x - xf) ** 3
    if x - x2 >= 0:
        row[Bz] = 1 / (6 * EIyy) * (x - x2) ** 3
    if x - x3 >= 0:
        row[Cz] = 1 / (6 * EIyy) * (x - x3) ** 3
    if x - xp >= 0:
        additional_sum += -Py / (6 * EIyy) * (x - xp) ** 3

    defl = row.dot(rf) + additional_sum

    return defl

# Twist
def ftheta(x):
    row = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.])
    additional_sum = locationvalue(x, fthetaqlst) / GJ

    if x - x1 >= 0:
        row[Ay] = 1 / GJ * (-zsc) * (x - x1)
    if x - xf >= 0:
        row[Fy] = 1 / GJ * (-zsc) * (x - xf)
    if x - x2 >= 0:
        row[By] = 1 / GJ * (-zsc) * (x - x2)
    if x - x3 >= 0:
        row[Cy] = 1 / GJ * (-zsc) * (x - x3)
    if x - xp >= 0:
        additional_sum += -Py / GJ * (-zsc) * (x - xp)

    twist = row.dot(rf) + additional_sum

    return twist