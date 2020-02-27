#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:33:38 2020

@author: Yann

Interpolation Verification
"""

import numpy as np
from Interpolation import patchinterpolate
from GraphsDistributedLoad import CreatePlots

def f(x,z):
    y = 10+5*x**(1/2)+6*z+16*x*z**2
    return y

def error_analysis(x, z, yint, f):
    actual_values = []
    for zi in z:
        actual_values_row = []
        for xi in x:
            actual_values_row.append(f(xi,zi))
        actual_values.append(actual_values_row)
    
    MSE = np.sum(abs(yint-actual_values)**2)/(np.shape(yint)[0]*np.shape(yint)[1])
    
    return MSE

xlst = np.array(np.linspace(0,1.611,41))
zlst = np.array(np.linspace(0,0.505,81))

data =[]
for zi in zlst:
    data_row = []
    for xi in xlst:
        y = f(xi,zi)
        data_row.append(y)
    data.append(data_row)
data = np.array(data)

x, z, yint = patchinterpolate(600,600,xlst,zlst,data)

MSE = error_analysis(x,z,yint,f)

print(MSE)



