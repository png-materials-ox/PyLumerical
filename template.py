# -*- coding: utf-8 -*-
"""
Created on Wed May  8 17:40:03 2024

@author: mans4209
"""

import numpy as np

class Template:
    def __init(self, roc=2, wlen=637e-09, fdtd=None):
        self.L_depth = .5 * self.wlen 
        self.feature_width = 2*np.sqrt(self.roc**2-(self.roc-self.L_depth)**2)
        self.xy_span = 1.5*self.feature_width
        self.xy_span_bleed = self.xy_span + 300e-09
        
        self.fdtd = fdtd
        
        if not self.fdtd:
            raise IOError ("A lumerical FDTD object must be defined in the constructor")
            