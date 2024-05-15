# -*- coding: utf-8 -*-
"""
Created on Mon May 13 11:28:21 2024

@author: mans4209
"""
import numpy as np
import scipy.constants as sc

class AnalysisBuilder:
    def __init__(self, fdtd=None, pdir='', monitor_names={}, min_wlen=637e-09, max_wlen=638e-09):
        
        self.fdtd = fdtd

        self.fdtd.load(pdir + 'fsp.fsp')
        self.fdtd.loaddata(pdir + 'lnp.ldf')
        
        self.monitor_names = monitor_names

        self.xy_span = fdtd.getnamed("structure::dielectric mediu", "x span")
        self.wlen = fdtd.getnamed("source", "center wavelength")

        self.x = self.fdtd.getdata(self.monitor_names['m2'],"x") 
        self.y = self.fdtd.getdata(self.monitor_names['m2'],"y")
        self.z = self.fdtd.getdata(self.monitor_names['m4'],"z") 

        self.x_pts = np.size(self.x)     # Points in the axis
        self.z_pts = np.size(self.z)

        self.x_res = self.xy_span/self.x_pts # resolution of simulation region
        self.midpoint = np.floor(self.x_pts/2)+1

        self.f = np.squeeze(self.fdtd.getdata(self.monitor_names['m2'],"f")) #Gets frequency data and removes singleton dimensions
        self.lam = sc.c/self.f

        ### Define the freq ROI for the peak interested ####
        
        self.min_wlen = min_wlen
        self.max_wlen = max_wlen
        
        self.lam_min_i = self.fdtd.find(self.lam, min_wlen)
        self.lam_max_i = self.fdtd.find(self.lam, max_wlen)
        
        # self.min_wlen = self.fdtd.find(self.lam, min_wlen)
        # self.max_wlen = self.fdtd.find(self.lam, max_wlen)
        
        self.w_range_max = (2*np.pi*sc.c)/self.min_wlen
        self.w_range_min = (2*np.pi*sc.c)/self.max_wlen
        self.f_range_max = sc.c/self.min_wlen
        self.f_range_min = sc.c/self.max_wlen
        