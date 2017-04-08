# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 14:16:14 2017

@author: sindrerb
"""

import numpy as np
import matplotlib.pyplot as plt

energy = np.loadtxt("./Dielectric-function/build/ENERGYFILE", dtype = 'float')
spectra = np.loadtxt("./Dielectric-function/build/SPECTRUMFILE", dtype = 'float')

plt.plot(energy,spectra[0,:])

"""
total = np.zeros((1,len(energy)))

for i in range(0,len(spectra[:,0])):
    total = spectra[i,:]
    plt.plot(energy,1.0/total)
"""    
    
    
    
    
    