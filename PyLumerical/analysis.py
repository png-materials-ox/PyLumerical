# -*- coding: utf-8 -*-

import numpy as np
import imp
# import os
import scipy.constants as sc
import cavityanalysis

# os.add_dll_directory("C:\\Program Files\\Lumerical\\v232\\api\\python\\")
lumapi = imp.load_source("lumapi","C:\\Program Files\\Lumerical\\v232\\api\\python\\lumapi.py")

fdtd = lumapi.FDTD()

pdir = 'C:\\Users\\mans4209\\Desktop\\testsim\\'

fdtd.load(pdir + 'fsp.fsp')
fdtd.loaddata(pdir + 'lnp.ldf')

# Monitor names
m0 = 'n'
m1 = 'xy_exoplanar'
m2 = 'xy_middle'
m3 = 'xy_exofeatured'
m4 = 'xz_middle'
m5 = 'xz_edge'
m6 = 'yz_middle'
m7 = 'yz_edge'

xy_span = fdtd.getnamed("structure::dielectric mediu", "x span")
wlen = fdtd.getnamed("source", "center wavelength")

x = fdtd.getdata(m2,"x") 
y = fdtd.getdata(m2,"y")
z = fdtd.getdata(m4,"z") 

x_pts = fdtd.size(x)     # Points in the axis
z_pts = fdtd.size(z)

x_res = int(xy_span/x_pts) # resolution of simulation region
midpoint = np.floor(x_pts[0]/2)+1

f = np.squeeze(fdtd.getdata(m2,"f")) #Gets frequency data and removes singleton dimensions
lam = sc.c/f

### Define the freq ROI for the peak interested ####

max_wlen = 638.0*1e-09
min_wlen = 630.0*1e-09

max_wlen = fdtd.find(lam, max_wlen)
min_wlen = fdtd.find(lam, min_wlen)
w_range_max = (2*np.pi*sc.c)/min_wlen
w_range_min = (2*np.pi*sc.c)/max_wlen
f_range_max = sc.c/min_wlen
f_range_min = sc.c/max_wlen

cavana = cavityanalysis.CavityAnalysis(fdtd=fdtd)
cavana.Qfactor()


# ## Decide whether to run the convergence test
# Q_test=0;

# fullrun=0;
# savedatafiles=1;
# savemonitor=1;

# lam_max=638.0; #nm
# lam_min=630.0; #nm


# ###################################################
# ### Define the freq ROI for the peak interested ####
# #lamoff = 0.2;
# lam_max_i=find(ax_lam,lam_max*1e-9);
# lam_min_i=find(ax_lam,lam_min*1e-9);
# w_range_max=2*pi*c/lam_min/1e-9;
# w_range_min=2*pi*c/lam_max/1e-9;
# f_range_max=c/lam_min/1e-9;
# f_range_min=c/lam_max/1e-9;
# ###################################################