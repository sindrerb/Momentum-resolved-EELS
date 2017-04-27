# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 16:30:41 2017

@author: sindrerb
"""

import numpy as np
import matplotlib.pyplot as plt

spectra = np.loadtxt("./MonteCarloSEELS/Build/SPECTRA", dtype = 'float')

energy = np.linspace(0,10,49)


for i, spectrum in enumerate(spectra):
    
    print i
    plt.plot(energy,spectrum[:])
    
plt.show()