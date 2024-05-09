# -*- coding: utf-8 -*-
import numpy as np
import scipy.constants as sc

class CavityAnalysis:
    def __init__(self, fdtd=None):
        self.fdtd = fdtd
        if not self.fdtd:
            raise IOError ("A lumerical FDTD object must be defined in the constructor")
    
    def Qfactor(self):    
        tlower = np.linspace(0.55,0,10);    # lower limit of time signal
        tupper = np.linspace(0.65,1,10);    # upper limit of time signal
    
        # simplify input variable names by removing spaces
        number_resonances = 1;
        
        min_filter_width = 1; # min width of filter in units of resonance FWHM
        max_filter_width = 6; # max width of filter in units of resonance FWHM
        filter_width_test_points = 20;
        zero_pad = 2**16; # fft zero padding
                         # Note fft zero pad should be a power of 2, 
                         # and larger gives more resolution in the 
                         # frequency domain.
        self.fdtd.groupscope("::model::Q monitors");
        
        
        # for(N=0; havedata("t_h_"+num2str(N+1)); N=N+1) { 1; } ##Determine how many monitors are there. 
        # 														## When havedata() is not true, it out puts 'Warning: Q_analysis_high_2.lsf line 33: in havedata, no data structures named t_h_19 were found'
        # 														##Then it records current N as the total number. 
        # 														##{1;} does not do anything. 
        
        num_mons = 5
        
        #################################################
        # get the monitor data for the first monitor
        #################################################
        t = self.fdtd.getdata("t_h_1","t");		## Getting the time axis.
        field0_t_Ex = self.fdtd.pinch(self.fdtd.getdata("t_h_1","Ex"));	##Getting the field evolution
        field0_t_Ey = self.fdtd.pinch(self.fdtd.getdata("t_h_1","Ey"));
        field0_t_Ez = self.fdtd.pinch(self.fdtd.getdata("t_h_1","Ez"));
        
        #################################################
        # do fft to frequency domain for all monitors
        #################################################
        w = self.fdtd.fftw(t,1,zero_pad);  #Angular frequency	axis					###Returns angular frequency. Value of zero_pad specifies the resolution of the frequency domain result.
        field_w = self.fdtd.matrix(len(w) ,6*num_mons);	

        ###Initialize a 2x2 matrix of length(w) and 6xNumber of monitors.
        for i in range(num_mons):
          mname = "t_h_" + self.fdtd.num2str(i)
          for j in range(5):
            if self.fdtd.almostequal(j,1): 
                component = "Ex"
            elif self.fdtd.almostequal(j,2): 
                component = "Ey"
            elif self.fdtd.almostequal(j,3): 
                component = "Ez"
            elif self.fdtd.almostequal(j,4): 
                component = "Hx"
            elif self.fdtd.almostequal(j,5): 
                component = "Hy"
            elif self.fdtd.almostequal(j,6): 
                component = "Hz"
            
            if j > 3.5: 
                extra_factor = np.sqrt(sc.mu0/sc.e0)
            else: 
                extra_factor = 1
            
            print('1')
            # if self.fdtd.havedata(mname, component):
              # field_w[0:len(w), 6*(i-1)+j:] = 2*extra_factor*((0:len(w)) <= (len(w)/2+0.1)) * self.fdtd.fft(self.fdtd.pinch(self.fdtd.getdata(mname,component)),1,zero_pad)
               
        
        
        # #################################################
        # # find resonant peaks, including all monitors
        # #################################################
        # w_i=find((w_range_min<w)&(w<w_range_max));
        # w_zm=w(w_i);
        
        # f_source=abs(sourcenorm(w/(2*pi)))^2*( (1:length(w)) <= (length(w)/2+0.1));
        # f_spectrum= sum(abs(field_w)^2,2)/f_source;
        
        # f_spectrum_zm=f_spectrum(w_i); #zm - zoom
        # f_spectrum_zm=f_spectrum_zm-min(f_spectrum_zm);
        
        # p_zm=findpeaks(f_spectrum_zm,number_resonances);
        # p=find(w,w_zm(p_zm));
        # f0_zm=w_zm(p_zm)/(2*pi);
        # f0=w(p)/(2*pi);
        
        # #w_i=find((spec_minangf<w)&(w<spec_maxangf));
        # #p = findpeaks(f_spectrum(w_i),number_resonances); 	## Find the largest peaks, number = number_resonances.
        # #f0 = w(w_i(p))/(2*pi);	## Find the central frequency
        # ##################################################
        # # find quality factors
        # #################################################
        
        # # reserve memory for results
        # peak_spectra = matrix(length(w),number_resonances);
        # peak_filters2 = matrix(length(w),number_resonances);
        
        # # calculate slope of decay using 40-60% of time signal
        # tp1 = round(tp1_lim(ranges)*length(t))+1;
        # tp2 = round(tp2_lim(ranges)*length(t));
        # t2 = t(tp1:tp2);
        # log_field_all = matrix(tp2-tp1+1,number_resonances); ## Selecting 40-60% of time axis.
        
        # Q = matrix(number_resonances);
        # delta_Q = matrix(number_resonances)+1e300;        ## Starting value of error in Q
        # slope_mean0=matrix(number_resonances);
        # slope_delta0=matrix(number_resonances);
        # # loop over each peak 
        # for(i=1:number_resonances) {
        #     # find FWHM of peak
        #     peak_val = f_spectrum_zm(p_zm(i));
        #     continue_search = 1;
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
        # }
        # #?'Centre = ' + num2str(c/f0);
        # #?'FWHM = ' + num2str((2*pi*c)/FWHM);
        # #?Q;
        # ?delta_Q;
        # #?slope_mean0;
        # #?slope_delta0;
        
        # groupscope("::model");
        
        # cd(filepath_data);##############################
        # #closeall;
        #     #titlebuff="Q at ranges=4 - Field decay vs t";
        #     #plot(t2,log_field,"time","log10(I) a.u.",titlebuff);
        #     #outbuff=titlebuff+".jpg";
        #     #selectfigure(1);
        #     #exportfigure(outbuff);
        
        #     #closeall;
        #     #titlebuff="Q at ranges=4 - slope (gradient) vs t";
        #     #plot(t2(1:(length(t2)-1)),slope,"time","a.u.",titlebuff);
        #     #outbuff=titlebuff+".jpg";
        #     #selectfigure(1);
        #     #exportfigure(outbuff);
            
        #     ## Plot Q spectrum and filtering
        #     p1 = find(w,0.8*min(w(p)));
        #     p2 = find(w,1.2*max(w(p)));
        #     f_Q = w/(2*pi);
        #     #closeall;
        #     titlebuff="Q at ranges=4 - spectrum filtering";
        #     plot(c/f_Q(p1:p2)*1e9,f_spectrum(p1:p2)/max(f_spectrum(p1:p2)),peak_filters2(p1:p2,1:number_resonances)
        #                  ,"Wavelength (nm)","Arbitrary units","Q Spectrum and filters");
        #     outbuff=titlebuff+".jpg";                     
        #     selectfigure(1);
        #     exportfigure(outbuff);   