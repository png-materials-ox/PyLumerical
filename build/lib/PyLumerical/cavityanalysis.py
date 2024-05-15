# -*- coding: utf-8 -*-
import numpy as np
import scipy.constants as sc

from scipy.fft import fft, ifft
from scipy.signal import find_peaks

class CavityAnalysis:
    def __init__(self, builder=None, fdtd=None):
        self.fdtd = fdtd
        if not self.fdtd:
            raise IOError ("A lumerical FDTD object must be defined in the constructor")
            
        self.builder = builder
    
    def Qfactor(self):    
        
        # Input properties
        number_resonances = 2 # Fill in the value
        make_plots = 0 # Fill in the value
        
        min_filter_width = 1  # min width of filter in units of resonance FWHM
        max_filter_width = 6  # min width of filter in units of resonance FWHM
        filter_width_test_points = 20
        zero_pad = 2 ** 16  # fft zero padding
        
        self.fdtd.groupscope("::model::Q monitors");
        
        # Load monitor data for the first monitor
        t = self.fdtd.getdata("t_h_1","t");		## Getting the time axis.        
        field0_t_Ex = np.squeeze(self.fdtd.getdata("t_h_1","Ex"));	##Getting the field evolution
        field0_t_Ey = np.squeeze(self.fdtd.getdata("t_h_1","Ey"));
        field0_t_Ez = np.squeeze(self.fdtd.getdata("t_h_1","Ez"));
        
        N = 4
        
        tp1_lim = np.linspace(0.55,0,10);    # lower limit of time signal
        tp2_lim = np.linspace(0.65,1,10);    # upper limit of time signal
        
        # Perform fft to frequency domain for all monitors
        # w = np.fft.fftfreq(len(t), d=t[1] - t[0])
        w = self.fdtd.fftw(t, 1, zero_pad).flatten()
        field_w = np.zeros((len(w), 6 * N), dtype=complex)
        for i in range(N):
            mname = "t_h_" + str(i + 1)
            for j in range(6):
                if j == 0:
                    component = "Ex"
                elif j == 1:
                    component = "Ey"
                elif j == 2:
                    component = "Ez"
                elif j == 3:
                    component = "Hx"
                elif j == 4:
                    component = "Hy"
                elif j == 5:
                    component = "Hz"
                    
                if j > 3.5: 
                    extra_factor = np.sqrt(sc.mu_0/sc.epsilon_0)
                else: 
                    extra_factor = 1
                if self.fdtd.havedata(mname, component):  # Assuming havedata() checks if data exists
                    extra = 2*extra_factor*(np.arange(0,len(w), 1) <= (len(w)/2+0.1))
                    ft = self.fdtd.fft(self.fdtd.pinch(self.fdtd.getdata(mname,component)),1,zero_pad)
                    field_w[:, (6*i)+j] = extra * ft.flatten()
                    
        
        #################################################
        # find resonant peaks, including all monitors
        #################################################

        w_i = np.where((self.builder.w_range_min  < w) & (self.builder.w_range_max > w))[0]
        w_zm = w[w_i]

        f_source = (abs(self.fdtd.sourcenorm(w/(2*np.pi)))**2).flatten() *( (np.arange(0, len(w), 1)) <= (len(w)/2+0.1))
        f_spectrum = np.sum(abs(field_w)**2,1)/f_source
        
        f_spectrum_zm = f_spectrum[w_i] #zm - zoom
        f_spectrum_zm = f_spectrum_zm - min(f_spectrum_zm)
        
        p_zm = self.fdtd.findpeaks(f_spectrum_zm, number_resonances).flatten()
        p_zm = [int(pz) for pz in p_zm]
        p = self.fdtd.find(w, w_zm[p_zm]).flatten()
        p = [int(pa) for pa in p]
        f0_zm = w_zm[p_zm]/(2*np.pi)
        f0 = w[p]/(2*np.pi)
        
        # reserve memory for results
        peak_spectra = self.fdtd.matrix(len(w),number_resonances)
        peak_filters2 = self.fdtd.matrix(len(w),number_resonances)
        
        ranges = 3
        # calculate slope of decay using 40-60% of time signal
        tp1 = int(np.round(tp1_lim[ranges]*len(t))+1)
        tp2 = int(np.round(tp2_lim[ranges]*len(t)))
        t2 = t[tp1:tp2]
        log_field_all = self.fdtd.matrix(tp2-tp1+1,number_resonances); ## Selecting 40-60% of time axis.
        
        Q = self.fdtd.matrix(number_resonances)
        delta_Q = self.fdtd.matrix(number_resonances)+1e300        ## Starting value of error in Q
        slope_mean0 = self.fdtd.matrix(number_resonances)
        slope_delta0 = self.fdtd.matrix(number_resonances)
            
        # loop over each peak 
        for i in range(number_resonances):
            # find FWHM of peak
            peak_val = f_spectrum_zm[p_zm[i]]
            continue_search = 1
        #     for(p1_zm=p_zm(i)-1; (p1_zm>1) & continue_search ; 1) {        ###Not sure what '&' means?
        #         if(f_spectrum_zm(p1_zm)<=peak_val/2) {         ### peak_val/2 -->Half maximum
        #             continue_search = 0; 
        #         } else {
        #             p1_zm = p1_zm-1;                ###p1 is the lower bond index of FWHM
        #         }
        #     }
        #     continue_search = 1;
        #     for(p2_zm=p_zm(i)+1; (p2_zm<length(w_zm))& continue_search; 1) {
        #         if(f_spectrum_zm(p2_zm)<=peak_val/2) { 
        #             continue_search = 0; 
        #         } else {
        #             p2_zm = p2_zm+1;            ### p2 is the upper bond index of FWHM
        #         }
        #     }
        #     if(p1_zm < 1) { p1_zm = 1; }
        #     if(p2_zm > length(w)) { p2_zm = length(w); }
        #     FWHM = w_zm(p2_zm)-w_zm(p1_zm);           
            
        #     for(filter_width=linspace(min_filter_width,max_filter_width,filter_width_test_points)) {
        #         # calculate the filter for the peak
        #         peak_filter = exp( -0.5*(w-w(p(i)))^2/(filter_width*FWHM)^2 ); ###Gaussian filter
        #                          ###This line does not work btw. w has different size than the filter_width
        #         # inverse fft to get data in time domain
        #         field2_t = 0;
        #         for(mcount=1:6*N) { 
        #             field2_t = field2_t + abs(invfft(pinch(field_w,2,mcount)*peak_filter))^2;
        #         }                        ###Summing all monitors in field_w to get time signal
        #                                  ###why not using the original data?
        #         field2_t = field2_t(tp1:tp2);        ### crop it with tp limits
        #         log_field = log10(abs(field2_t));
        
        #         # calculate slope and Q from the slope of the decay
        #         # estimate error from the slope
        #         slope = (log_field(2:length(t2))-log_field(1:length(t2)-1))/
        #             (t(2:length(t2))-t(1:length(t2)-1));
        #         slope_mean = sum(slope)/length(slope);
        #         slope_delta = sqrt( sum((slope-slope_mean)^2)/length(slope) );
        #         Q_test = -w(p(i))*log10(exp(1))/(slope_mean);
        #         delta_Q_test = abs(slope_delta/slope_mean*Q_test);
        #         if(delta_Q_test < delta_Q(i)) {
        #             Q(i) = Q_test;
        #             delta_Q(i) = delta_Q_test;
        #             slope_mean0(i)=slope_mean;
        #             slope_delta0(i)=slope_delta;
                    
        #             # collect data for final plot
        #             peak_spectra(1:length(w),i) = f_spectrum * peak_filter^2;
        #             peak_filters2(1:length(w),i) = peak_filter^2;
        #             log_field_all(1:length(t2),i) = log_field;
        
        #         }
        #     }
        #     # output summary of peak results to script window
        #     ?"Resonance " + num2str(i) + ":";
        #     ?"    frequency = " + num2str(w(p(i))/(2*pi)*1e-12) + "THz, or "+num2str(2*pi*c/w(p(i))*1e9)+" nm";
        #     ?"    Q = " + num2str(Q(i)) +" +/- " + num2str(delta_Q(i));
        #     ?"    FWHM = " + num2str(FWHM/(2*pi)) + "Hz, or " + num2str(((2*pi*c)/FWHM)*1e09) +" nm";     
        
        # #?'Centre = ' + num2str(c/f0);
        # #?'FWHM = ' + num2str((2*pi*c)/FWHM);
        # #?Q;
        # ?delta_Q;
        # #?slope_mean0;
        # #?slope_delta0;
        
    # groupscope("::model");
            