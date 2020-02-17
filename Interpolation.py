"""
Interpolation of Aerodynamic load

Created by B. van Dillen
"""

#Imports
import numpy as np


# Necessary inputs
ca = 0.505 #m
la = 1.611 #m


# Reading the aerodynamic data
data = np.loadtxt('AERO.dat',delimiter = ',')


# Calculate locations of x and z coordinates
Nz = len(data)
Nx = len(data[0])

zlst = []
xlst = []

# Calculate the z coordinates
for i in range(1, Nz + 1):
    thetazi = (i - 1) / Nz * np.pi
    thetazi1 = (i) / Nz * np.pi
    zi = -(1 / 2) * ((ca / 2) * (1 - np.cos(thetazi)) + (ca / 2) * (1 - np.cos(thetazi1)))
    zlst.append(zi)

# Calculate the x coordinates
for i in range(1, Nx + 1):
    thetaxi = (i - 1) / Nx * np.pi
    thetaxi1 = (i) / Nx * np.pi
    xi = (1 / 2) * ((la / 2) * (1 - np.cos(thetaxi)) + (la / 2) * (1 - np.cos(thetaxi1)))
    xlst.append(xi)


# Creating the 2D mesh
X, Z = np.meshgrid(xlst, zlst)

# Create lists to save the coordinates and A matrices
rectangle_lst = []
A_lst = []

# For loop going through every row
for i in range(Nz-1):
    rectangle_row = []
    A_row = []
    # For loop going through every column
    for j in range(Nx-1):
        # Determine x1, x2, z1, z2
        x1 = X[i][j]
        x2 = X[i][j+1]
        z1 = Z[i][j]
        z2 = Z[i+1][j]

        # Determine points of the rectangle
        points_matrix = [[x1,z1],
                         [x2,z2]]
        # Append to the row list
        rectangle_row.append(points_matrix)

        # Determine A matrix
        A = [[1,x1,z1,x1*z1],
             [1,x2,z1,x2*z1],
             [1,x1,z2,x1*z2],
             [1,x2,z2,x2*z2]]
        # Append to the row list
        A_row.append(A)

    # Append the rows to the lists
    rectangle_lst.append(rectangle_row)
    A_lst.append(A_row)

# Transform into a numpy array
rectangle_lst = np.array(rectangle_lst)
A_lst = np.array(A_lst)
