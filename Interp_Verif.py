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
    y = 10+5*x+6*z+16*x*z
    return y

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

x, z, yint = patchinterpolate(10,10,xlst,zlst,data)

#CreatePlots(xlst,zlst,data)

#CreatePlots(x,z,yint)

