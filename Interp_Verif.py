#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:33:38 2020

@author: Yann

Interpolation Verification
"""

import numpy as np

def f(x,z):
    y = 10+5*x+6*y+16*x*y
    return y

xlst = np.linspace(0,1.611,41)
zlst = np.linspace(0,0.505,81)

data =[]
for zi in zlst:
    data_row = []
    for xi in xlst:
        y = f(xi,zi)
        data_row.append(y)
    data.append(data_row)
