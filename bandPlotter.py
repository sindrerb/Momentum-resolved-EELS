# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 09:59:53 2017

@author: sindrerb
"""
import numpy as np
import matplotlib.pyplot as plt

basis = np.loadtxt("./Kronig-Penney-Model/Build/BASISFILE", dtype = 'float')
waves = np.loadtxt("./Kronig-Penney-Model/Build/WAVEFILE", dtype = 'float')
potential = np.loadtxt("./Kronig-Penney-Model/Build/POTENTIALFILE", dtype = 'float', skiprows = 1)[0]

Q = waves[:,0:3]
E = waves[:,3]


for i,q in enumerate(Q):
    plt.scatter(q[0],E[i])
    
plt.show()