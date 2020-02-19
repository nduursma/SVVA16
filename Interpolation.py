import numpy as np
import matplotlib.pyplot as plt
from math import pi 
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt('AERO.dat',delimiter = ',')

def patchinterpolate(x_mesh, z_mesh, data):
    
    #GETTING X AND Z VALUES OF DATA
    
    xlst = []
    zlst = []
    Nz = len(data)
    Nx = len(data[0])
    ca = 0.505
    la = 1.611
    
    for i in range(1,Nz+1):  
        thetazi  = (i-1)/Nz*pi
        thetazi1 = (i)/Nz*pi
        z  =  -(1/2)*((ca/2)*(1-np.cos(thetazi))+(ca/2)*(1-np.cos(thetazi1)))    
        zlst.append(z)
    
    for k in range(1,Nx+1):
        thetaxk  = (k-1)/Nx*pi
        thetaxk1 = (k)/Nx*pi
        x  = (1/2)*((la/2)*(1-np.cos(thetaxk))+(la/2)*(1-np.cos(thetaxk1)))   
        xlst.append(x)
      
    #INTERPOLATION

    #GENERATING LIST OF ALL PATCHES
    squarelst =[]
    
    for j in range(len(data)-1):
        rowlst = []
        for i in range(len(data[0])-1):
            mat = np.array([[xlst[i],zlst[j]], [xlst[i+1],zlst[j+1]]])
            rowlst.append(mat)
        squarelst.append(rowlst)
    
    
    #GENERATING A MATRIX   
    
    Alst = []
    
    for j in range(len(data)-1):
        Arowlst = []
        for i in range(len(data[0])-1):
            x0 = squarelst[j][i][0][0]
            z0 = squarelst[j][i][0][1]
            x1 = squarelst[j][i][1][0]
            z1 = squarelst[j][i][1][1]
            mat = np.array([[1, x0, z0, x0*z0],
                            [1, x1, z0, x1*z0],
                            [1, x0, z1, x0*z1],
                            [1, x1, z1, x1*z1]])
            Arowlst.append(mat)
        Alst.append(Arowlst)
    
    #GENERATING LIST OF ALL FUNCTION VALUES
    f_lst = []
    for j in range(len(data)-1):
        f_row = []
        for i in range(len(data[0])-1):
            f_11 = data[j][i]
            f_21 = data[j][i+1]
            f_12 = data[j+1][i]
            f_22 = data[j+1][i+1]
            vec = np.array([[f_11],[f_21],[f_12],[f_22]])
            f_row.append(vec)
        f_lst.append(f_row)
    
    #GENERATING LIST OF a:
    a_lst = []
    for j in range(len(data)-1):
        a_row = []
        for i in range(len(data[0])-1):
            a_vec = np.linalg.solve(Alst[j][i], f_lst[j][i])
            a_row.append(a_vec)
        a_lst.append(a_row)
    
    #GENERATING POINTS TO BE INTERPOLATED
    x_lst = np.linspace(xlst[0],xlst[-1],x_mesh)
    z_lst = np.linspace(zlst[0], zlst[-1],z_mesh)
    
    #GENERATING INTERPOLATED VALUES
    
    pxz_lst = []
               
    
    for j in range(len(zlst)-1):
        z1 = zlst[j]

        z2 = zlst[j+1]

        for k in range(len(z_lst)):

            pxz_row = []
            
            z_coordinate = z_lst[k]

            for i in range(len(xlst)-1):

                x1 = xlst[i]

                x2 = xlst[i+1]

                

                for l in range(len(x_lst)):

                    x_coordinate = x_lst[l]

                    if x1 <= x_coordinate <= x2 and z2 <= z_coordinate <= z1:

                        pxz = a_lst[j][i][0] + a_lst[j][i][1]*x_coordinate + a_lst[j][i][2]*z_coordinate + a_lst[j][i][3]*x_coordinate*z_coordinate

                        pxz_row.append(float(pxz))

            if len(pxz_row) != 0:        

                pxz_lst.append(pxz_row)
                
    return x_lst, z_lst, pxz_lst
