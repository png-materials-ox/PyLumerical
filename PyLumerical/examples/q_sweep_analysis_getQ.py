# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:50:34 2024

@author: mans4209
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 18:23:45 2024

@author: mans4209
"""

# -*- coding: utf-8 -*-

import numpy as np
import imp
# import os
import scipy.constants as sc
from PyLumerical import cavityanalysis, lumericalanalysis
import matplotlib.pyplot as plt
import progressbar
import json
# import analysisbuilder

# os.add_dll_directory("C:\\Program Files\\Lumerical\\v232\\api\\python\\")
lumapi = imp.load_source("lumapi","C:\\Program Files\\Lumerical\\v232\\api\\python\\lumapi.py")

fdtd = lumapi.FDTD()

basedir = 'C://Users//mans4209//Documents//LumericalFiles//monolithic_cavity_sweep//sweep_q_mesh1//'

qlist = [4, 5]

Qmax = []
err = []
Vlam = []
collection_eff = []

with progressbar.ProgressBar(max_value=len(qlist)) as bar:
    for i, q in enumerate(qlist):
    
        pdir = basedir + '20240605_monolithic_diamond_cavity_roc_q{:d}_Lsub20nm_theta90//'.format(q)
        fdtd.load(pdir + 'fsp.fsp')
        analysis = lumericalanalysis.LumericalAnalysis(fdtd=fdtd)
        
        res = analysis.resonances()
        Qmax.append(max(res['Q'].values))
        err.append(res['error'].iloc[res['Q'].argmax()])
        
        mode_vol = analysis.mode_volume_2D()
        
        Vlam.append(np.mean(mode_vol['Vol_lam_avg']))
        
        farfield = analysis.farfield_analysis(NA=0.9)
        
        collection_eff.append(np.mean((farfield['T34'])))
    
        with open('sweep_data_q{:d}.json'.format(q), 'w', encoding='utf-8') as f:
            json.dump({'f0': res['f0'].values.tolist(),
                        'decay': res['decay'].values.tolist(),
                        'Q': res['Q'].values.tolist(),
                        'amp': res['amp'].values.tolist(),
                        'phase': res['phase'].values.tolist(),
                        'error': res['error'].values.tolist(),
                        
                        'Vol_abs_xz':mode_vol['Vol_abs_xz'].values.tolist(),
                        'Vol_lam_xz':mode_vol['Vol_lam_xz'].values.tolist(),
                        'Vol_abs_yz':mode_vol['Vol_abs_yz'].values.tolist(),
                        'Vol_lam_yz':mode_vol['Vol_lam_yz'].values.tolist(),
                        'Vol_abs_avg':mode_vol['Vol_abs_avg'].values.tolist(),
                        'Vol_lam_avg':mode_vol['Vol_lam_avg'].values.tolist(),
                        
                        'f':farfield['f'].tolist(),
                        'sp':farfield['sp'].tolist(),
                        'dp':farfield['dp'].tolist(),
                        'wlen':farfield['wlen'].tolist(),
                        'purcell':farfield['Fp'].tolist(),
                        'cone_angle':farfield['cone_angle'].tolist(),
                        'collection_eff':farfield['T34'].tolist(),
                        'theta':farfield['theta'].tolist(),
                        'P_v_theta':farfield['P_v_theta'].tolist(),
                        'Twlen':farfield['Twlen'].tolist(),
                        'trans':farfield['trans'].tolist(),
                        
                        
                    }, f, ensure_ascii=False, indent=4)
        
        bar.update(i)
    
plt.plot(qlist, Qmax, 's', color='k')
plt.xlabel('q')
plt.ylabel('Q Factor')
plt.show()

plt.plot(qlist, Vlam, 's', color='k')
plt.xlabel('q')
plt.ylabel('Normalized Mode Volume')
plt.show()

plt.plot(qlist, collection_eff, 's', color='k')
plt.xlabel('q')
plt.ylabel('Collection Efficiency')
plt.show()

