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


def integrate_and_move_forces(x_mesh, z_mesh, data):
    #Retrieving total resultant force, list if centroids and list of magnitudes
    Vtot, Centroid_Zlst, Maglst  = integrate(x_mesh, z_mesh, data)
    
    #creating list to store moments
    Momentmaglst = []
    
    for num1, num2 in zip(Centroid_Zlst, Maglst):
        #calculating moments, equal to each force in Maglst multiplied by their arms --> centroids
        Momentmaglst.append(num1 * num2)
    
    return Vtot, Centroid_Zlst, Maglst, Momentmaglst
    