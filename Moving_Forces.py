# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 12:01:11 2020

@author: Raven
"""

import numpy as np

from Resultant_Centroid_DistributedLoadCorrect import integrate
data = np.loadtxt('AERO.dat',delimiter = ',')
Vtot, Centroid_Zlst, Maglst  = integrate(16,16,data)