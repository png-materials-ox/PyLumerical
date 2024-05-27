# -*- coding: utf-8 -*-
"""
Created on Sun May 26 16:30:57 2024

@author: mans4209
"""
import matplotlib.pyplot as plt
from  matplotlib import patches
from matplotlib.figure import Figure
from matplotlib import rcParams
import numpy as np
import pandas as pd
import scipy.constants as sc

class LumericalAnalysis:
    
    def __init__(self, fdtd=None):
        self.fdtd = fdtd
        if not self.fdtd:
            raise IOError ("A lumerical FDTD object must be defined in the constructor")
    
    def electric_field_magnitude(self) -> pd.DataFrame:
        t = self.fdtd.getdata('Qanalysis::t1', 't').flatten()
        Ex = np.squeeze(self.fdtd.getdata('Qanalysis::t1', 'Ex'))
        Ey = np.squeeze(self.fdtd.getdata('Qanalysis::t1', 'Ey'))
        Ez = np.squeeze(self.fdtd.getdata('Qanalysis::t1', 'Ez'))
        
        E = np.sqrt(Ex**2 + Ey**2 + Ez**2)
        return pd.DataFrame(data={'t':t, 'Ex':Ex, 'Ey':Ey, 'Ez':Ez, 'E':E})    
    
    def f_spectrum(self, plotting=False, saveplot=False) -> list:
        if not self.fdtd.havedata('Qanalysis'):
            self.fdtd.runanalysis
            
        # wlen = self.fdtd.getresult('Qanalysis', 'w').flatten()
        f_spectrum = self.fdtd.getresult('Qanalysis', 'f_spectrum').flatten()
        
        if plotting:
            plt.plot(f_spectrum)
            plt.ylabel('Response (AU)')
            if saveplot:
                plt.savefig('f_spectrum.jpg', bbox_inches="tight")
            plt.show()
            
        return f_spectrum
    
    def omega(self):
        if not self.fdtd.havedata('Qanalysis'):
            self.fdtd.runanalysis     
        return self.fdtd.getresult('Qanalysis', 'w').flatten()
        
    def resonances(self) -> pd.DataFrame:
        if not self.fdtd.havedata('Qanalysis'):
            self.fdtd.runanalysis     
        
        res = self.fdtd.getresult('Qanalysis', 'resonances')
        df = pd.DataFrame(data={'f0':res[:,0], 'decay':res[:,1], 'Q':res[:,2],'amp':res[:,3], 'phase':res[:,4], 'error':res[:,5]})
        return df
    
    def complex_eigenfrequency(self, omega0=1, Q=100):
        omega_imag = 2*Q
        return omega0 - 1j * omega_imag
    
    # def transfer_function(self, eigenfreqs=[], w=[], wr=[], Q=[]):
    #     r = 1
    #     b = 0
        
    #     return (r)/(1j* (w - wr) - (wr/(2*Q)))
        
        
    
    def farfield_analysis(self):
        # start_unix_time = now;
        if not self.fdtd.havedata('farfield'):
            self.fdtd.runanalysis
        
        # savebool = 'no';
        
        # s=splitstring(pwd, '/');
        # if(s{1}=='C:'){    
        #     fname = s{length(s)};    
        #     ?fname;
        # }
        # else{
        #     save == "no";
        #     ?"Save directory is not in the C drive. File not saved";
        # }
        #loaddata('C:\Users\mans4209\Documents\Lumerical Files\SIL_Simulations\fsf\workspace_data\rect_test_2_monitor_0p1um_below.ldf');
        # start_farfield_unix_time = now;
        
        ##############################################
        # Plot results from single simulation
        
        ## Scaling factors for plots
        um = 1e6;
        nm = 1e9;
        
        ##############################################
        # Get field from y-normal monitor
        #x = getdata("xz_profile", "x");
        #z = getdata("xz_profile", "z");
        #E2 = getelectric("xz_profile");
        ##x = getdata("xz_mid", "x");
        ##z = getdata("xz_mid", "z");
        ##E2 = getelectric("xz_mid");
        #image(x*um, z*um, pinch(E2, 4, 10), "x (um)", "z (um)", "", "logplot");
        
        ##############################################
        # Get field from z-normal monitor (farfield::z2)
        x = self.fdtd.getdata("farfield::z2", "x")
        y = self.fdtd.getdata("farfield::z2", "y")
        E2 = self.fdtd.getelectric("farfield::z2")
        self.fdtd.image(x*um, y*um, self.fdtd.pinch(E2, 4, 1), "x (um)", "y (um)", "")
        
        ##############################################
        # Get results from "farfield" analysis group
        ## Total transmission
        T = self.fdtd.getresult("farfield", "T");
        plt.plot(T['lambda']*nm, T['T']) 
        plt.xlabel('$\lambda$ (nm)')
        plt.ylabel('Total transmission')
        plt.show()
        
        ## Poynting vector in far field
        S = self.fdtd.getresult("farfield", "farfield")
        self.fdtd.visualize(S)
        
        ## P_vs_theta
        P_vs_theta = self.fdtd.getresult("farfield", "P_vs_theta");
        plt.plot(P_vs_theta['theta_degrees'], self.fdtd.pinch(P_vs_theta['P'], 2, 1))
        plt.xlabel('theta (degrees)')
        plt.ylabel('Normalized Power')
        plt.show()
        
        ## Purcell factor
        Purcell = self.fdtd.getresult("farfield", "Purcell")
        plt.plot(Purcell['lambda']*nm, Purcell['purcell'])
        plt.xlabel('lambda (nm)')
        plt.ylabel('Purcell factor')
        
        ## Power in angular cone
        # create 2D mesh of theta
        NA = 0.55
        cone_angle = np.arcsin(0.55)*180/np.pi
        Theta = self.fdtd.meshgridx(P_vs_theta['theta_radians'],P_vs_theta['f'])
        # integrate over specified cone
        T34 = 0.5*2*np.pi*np.real(self.fdtd.integrate(P_vs_theta['P']*self.fdtd.sin(Theta)*(Theta<=cone_angle*np.pi/180),1,P_vs_theta['theta_radians']))
        
        # plot final results
        plt.plot(P_vs_theta['lambda']*nm,T34, label="Normalized power")
        plt.plot(P_vs_theta['lambda']*nm,T34/Purcell['purcell'], label="Optical extraction efficiency")
        plt.xlabel('wavelength (nm)')
        plt.legend()
        
        # cone_angle = 0:1:70;
        # T34_m = matrix(length(P_vs_theta.lambda),length(cone_angle));
        
        # for(i=1:length(cone_angle)){ #degrees
        #     # integrate over specified cone
        #     #print(cone_angle(i));
        #     T34_m(:,i) = 0.5*2*pi*real(integrate(P_vs_theta.P*sin(Theta)*(Theta<=cone_angle(i)*pi/180),1,P_vs_theta.theta_radians));
        # }
        
        # end_unix_time = now;
        
        #  ?"Farfield simulation took " + num2str((end_unix_time - start_farfield_unix_time)/60) + " minutes";
        #  ?"Total simulation took " + num2str((end_unix_time - start_unix_time)/60) + " minutes";
        
        # ##### SAVE ####
        # #get_datetime;
        # ?fname;
        # select("::model::hemisph_surf");
        # r = get("radius");
        
        # dx = getdata('dipole_source','x');
        # dy = getdata('dipole_source','y');
        # dz = getdata('dipole_source','z');
        
        # dx = pinch(dx,1);
        # dy = pinch(dy,1);
        # dz = pinch(dz,1);
        
        # #>>>> Save Data <<<<
        # if (savebool=="yes"){
        #     ff = getresult('farfield', 'farfield');
        #     f = getdata('xy_profile','f');
        #     #f = getdata('xy_mid','f');
        #     wlen = P_vs_theta.lambda;
        #     Fp = Purcell.purcell;
        #     theta = P_vs_theta.theta_degrees;
        #     P_v_theta = pinch(P_vs_theta.P, 2, 1);
        #     Twlen = T.lambda*nm;
        #     trans = T.T;
        #     sp = sourcepower(f);
        #     dp = dipolepower(f);
        #     plot(wlen, trans*(sourcepower(f)/dipolepower(f)));
            
        #     vtksave(fname + ".vtu",ff);
        #     matlabsave(fname + ".mat", f, sp, dp, wlen, Fp, cone_angle, T34, T34_m, theta, P_v_theta, Twlen, trans);
        # }else if (savebool == "no"){
        #     ?"Data not saved";
        # }else{
        #     ?"Save toggle incorrect. Use yes or no";
        # }
        
    
    def mode_volume_2D_QNM(self, dipole_shift=0):
    
        x = self.fdtd.getdata('mode volume 2D::xy_middle',"x").flatten()    
        y = self.fdtd.getdata('mode volume 2D::xy_middle',"y").flatten()
        z = self.fdtd.getdata('mode volume 2D::xz_middle',"z").flatten()    
    
        x_pt = np.shape(x)           # Points in the axis
        z_pt = np.shape(z)           # Points in the z axis            
    
        f = self.fdtd.getresult('mode volume 2D::xy_middle', 'f');
        wlen = sc.c/f;
        midpoint_z_ind = int(self.fdtd.round(.5*self.fdtd.size(x , 1)))    # find the z midpoint index value
        midpoint_x_ind = int(self.fdtd.floor(x_pt[0]/2)+1)

        # Integration ranges (0 -> rmax, zmin -> zmax)
        z_int_range = z                                               
        r_int_range = x[0:midpoint_x_ind]
        
        # Get electric fields
        E_xz = self.fdtd.pinch(self.fdtd.getelectric('mode volume 2D::xz_middle'))
        E_yz = self.fdtd.pinch(self.fdtd.getelectric('mode volume 2D::yz_middle'))
        
        # Get refractive indices
        n_xz = self.fdtd.pinch(self.fdtd.getdata('mode volume 2D::n',"index_z"))                       
        
        nv_zpos_ind = int(self.fdtd.find(z, dipole_shift))           # Position of the NV
        n_nv_xz = np.real(n_xz[midpoint_x_ind, nv_zpos_ind]) # Refractive index at the position of the NV
        n_nv_yz = np.real(n_xz[midpoint_x_ind, nv_zpos_ind])                                  
        
        # Predefine arrays for looping
        Vol_abs_xz = np.ones(len(f))
        Vol_lam_xz = np.ones(len(f))
        
        Vol_abs_yz = np.ones(len(f))
        Vol_lam_yz = np.ones(len(f))
        
        # Loop over each frequency
        for i in range(len(f)):
            n = n_xz[0:midpoint_x_ind,0:z_pt[0]]
        
            E_xz_res = self.fdtd.pinch(E_xz, 3, 1)
            E_xz_mid = E_xz_res[0:midpoint_x_ind,0:z_pt[0]]                      # half of the x-z cut
            eps_E_xz = np.real(n**2)*E_xz_mid

            # Calculate epsilon * E at the NV position 
            eps_E_xz_at_nv = eps_E_xz[midpoint_x_ind-1, nv_zpos_ind]                          
            eps_E_xz_max = max(eps_E_xz.flatten())    # Max value                                
            
            E_yz_res = self.fdtd.pinch(E_yz, 3, 1)
            E_yz_mid = E_yz_res[0:midpoint_x_ind,0:z_pt[0]]                      # half of the x-z cut
            eps_E_yz = np.real(n**2)*E_yz_mid
            eps_E_yz_at_nv = eps_E_yz[midpoint_x_ind-1, nv_zpos_ind]                         
            eps_E_yz_max = max(eps_E_yz.flatten())
                
            
            V0_xz = E_xz_mid*np.real(n**2)*2*np.pi                                  # n^2 
            Vol_raw1_xz = self.fdtd.integrate(V0_xz,2,z_int_range)
            Vol_raw2_xz = 1e18*abs(self.fdtd.integrate(Vol_raw1_xz*r_int_range,1,r_int_range)) # Unit in um^3.
            
            Vol_abs_xz[i] = Vol_raw2_xz/eps_E_xz_max                         # unit:um^3, for air-like mode
            Vol_lam_xz[i] = Vol_abs_xz[i]/(wlen[i]^3*1e18)  
            
            V0_yz = E_yz_mid*np.real(n**2)*2*np.i                                  # n^2 
            Vol_raw1_yz = self.fdtd.integrate(V0_yz,2,z_int_range)
            Vol_raw2_yz = 1e18*abs(self.fdtd.integrate(Vol_raw1_yz*r_int_range,1,r_int_range)) # Unit in um^3.
            
            Vol_abs_yz[i] = Vol_raw2_yz/(eps_E_yz_max)                         # unit:um^3, for air-like mode
            Vol_lam_yz[i] = Vol_abs_yz(i)/(wlen(i)^3*1e18)   
        
        
        ## Average mode volume
        Vol_abs_avg = (Vol_abs_xz+Vol_abs_yz)/2
        Vol_lam_avg = (Vol_lam_xz+Vol_lam_yz)/2
        n_nv_avg = (n_nv_xz+n_nv_yz)/2
        In2_nv_avg = (eps_E_xz_at_nv+eps_E_yz_at_nv)/2
        In2_max_avg = (eps_E_xz_max+eps_E_yz_max)/2          # unit:lam^3
        
        return Vol_abs_avg
