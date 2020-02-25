#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 14:00:35 2020

@author: Yann

Main Simulation Script
"""

import numpy as np
from Interpolation import patchinterpolate
from NEW_Forces_Deflections import output, locationvalue

data = np.loadtxt('AERO.dat',delimiter = ',')

sc = -0.085
xlst, zlst, qlst = patchinterpolate(600,600,data)

Vlst, Mlst, defllst, Tlst, thetalst = output(xlst, zlst, qlst, sc)

