# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 12:01:11 2020

@author: Raven
"""

"""This script moves the forces to the hinge line, and calculates the corresponding moment"""

import numpy as np

from Resultant_Centroid_DistributedLoadCorrect import integrate

#Read data of distributed load magnitude
data = np.loadtxt('AERO.dat',delimiter = ',')  #Creating an Array from the aerodynamic load data file

Vtot, Centroid_Zlst, Maglst  = integrate(16,16,data)

def integrate_and_move_forces(x_mesh, z_mesh, data):
    Vtot, Centroid_Zlst, Maglst  = integrate(x_mesh, z_mesh, data)
    
    Momentmaglst = []
    
    for num1, num2 in zip(Centroid_Zlst, Maglst):
        Momentmaglst.append(num1 * num2)
    
    return Vtot, Centroid_Zlst, Maglst, Momentmaglst
    