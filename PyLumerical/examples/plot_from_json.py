# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 18:28:31 2024

@author: mans4209
"""

import numpy as np
import matplotlib.pyplot as plt
import json

qlist = [4, 8, 12, 16, 20, 28, 32]
# qlist = [28]

max_ang = []
Q = []

for i, q in enumerate(qlist):      
    file = open('coarse_sweep/sweep_data_q{:d}.json'.format(q), 'r')
    datafile = json.load(file)
        
    theta = np.array(datafile['theta'])
    P = np.array(datafile['P_v_theta'])
    
    max_ang.append(theta[np.where(P == max(P))][0])

    Q.append(max(np.array(datafile['Q'])))
    
    
plt.plot(qlist, max_ang, 's', color='k')
plt.xlabel('q')
plt.ylabel('Emission Angle of Peak Intensity (Deg)')
plt.show()

plt.plot(qlist, Q, 's', color='k')
plt.xlabel('q')
plt.ylabel('Q Factor')
plt.show()