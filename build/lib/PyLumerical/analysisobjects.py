# -*- coding: utf-8 -*

class AnalysisObjects:
    
    def __init__(self, fdtd=None, mon=None):
        self.fdtd = fdtd
        # self.sim = simulation.Simulation(fdtd=fdtd, xy_span_bleed=xy_span, z_min=0, z_max=1e-06)
        self.mon = mon
        
    def farfield(self, xy_span=1e-06, z_span=1e-06, x=0, y=0, z=0, theta_max=90,
                 N_theta=180, Nphi=91):
        
        self.fdtd.unselectall()
        self.fdtd.addobject('scat_ff_radiation')
        self.fdtd.select('ff_radiation_plot') 
        self.fdtd.set('name','farfield')
        self.fdtd.set('x span', xy_span)
        self.fdtd.set('y span', xy_span)
        self.fdtd.set('z span', z_span)
        self.fdtd.set('x', x);
        self.fdtd.set('y', y);
        self.fdtd.set('z', z) # Need to work out how to find the direct centre of the overall structure

        # with open('./farfield_script.txt', 'r', encoding="utf8") as f:
        #     script = f.read().rstrip()
        
        script = """            
            # setup downsampling for far field projections
            farfieldsettings("override near field mesh",true);
            near_points = %near field points per wavelength%;
            if(near_points < 2) {
                ?"Warning, you can't set near field points per wavelength to less than 2";
                near_points = 2;
            }
            if(near_points > 10) {
                ?"Warning, you can't set near field points per wavelength to more than 10";
                near_points = 10;
            }
            farfieldsettings("near field samples per wavelength",near_points);
            
            ##############################################
            # automatically unfold field data if symmetry BC is applied
            # collect power transmission data
            if (havedata("x1", "f")) {
              symm_x = 0;
              Tx1 = -transmission("x1");
            } else {
              xtemp = getdata("y2", "x");
              ztemp = getdata("y2", "z");
              Eztemp = pinch(getdata("y2", "Ez"));
              
              Ez2mid = sum(Eztemp(round(length(xtemp)/2), 1:length(ztemp))^2);
              if (Ez2mid != 0) {
                symm_x = 1;
              } else {
                symm_x = -1;
              }
              Tx1 = transmission("x2");
            }
            
            if (havedata("y1", "f")) {
              symm_y = 0;
              Ty1 = -transmission("y1");
            } else {
              ytemp = getdata("x2", "y");
              ztemp = getdata("x2", "z");
              Eztemp = pinch(getdata("x2", "Ez"));
              
              Ez2mid = sum(Eztemp(round(length(ytemp)/2), 1:length(ztemp))^2);
              if (Ez2mid != 0) {
                symm_y = 1;
              } else {
                symm_y = -1;
              }
              Ty1 = transmission("y2");
            }
            if(include_z_min_boundary) {
                if (havedata("z1", "f")) {
                    symm_z = 0;
                    Tz1 = -transmission("z1");
                } else {
                    xtemp = getdata("y2", "x");
                    ztemp = getdata("y2", "z");
                    Eytemp = pinch(getdata("y2", "Ey"));
                    
                    Ey2mid = sum(Eytemp(1:length(xtemp), round(length(ztemp)/2))^2);
                    if (Ey2mid != 0) {
                        symm_z = 1;
                    } else {
                        symm_z = -1;
                    }
                    Tz1 = transmission("z2");
                }
            } else {
                Tz1 = 0;
            }
            f = getdata("x2","f");   # get freqency data
            Tx2 = transmission("x2"); # transmission data
            Ty2 = transmission("y2");
            Tz2 = transmission("z2");
            Ttotal = Tx1+Tx2+Ty1+Ty2+Tz1+Tz2;
            T = matrixdataset("T");
            T.addparameter("lambda",c/f,"f",f);
            T.addattribute("T",Ttotal);
            ##############################################
            
            ##############################################
            # Identify the closest wavelength to target:
            if(!%project all wavelengths%) {
                target_wavelength = %target wavelength%;
                i_target = find(f,c/target_wavelength);
                ?"Target wavelength = " + num2str(target_wavelength);
                ?"Wavelength used = " + num2str(c/f(i_target));
                f = f(i_target);
            } else {
                i_target = 1:length(f);
                ?"Calculating far field of all wavelengths";
            }
            ##############################################
            
            ##############################################
            # define the theta and phi grid
            theta = linspace(%theta min%,%theta max%,%N theta%)*pi/180; # user-modifiable in the Variables tab
            phi = linspace(0,180,%N phi per half%)*pi/180; # user-modifiable in the Variables tab
            if(almostequal(phi(1)+2*pi,phi(end))) {
                phi = phi(1:end-1);
                wrap_phi = true;
            } else {
                wrap_phi = false;
            }
            nt = length(theta);
            np = length(phi);
            nf = length(f);
            Theta = meshgridx(theta,phi);
            Phi = meshgridy(theta,phi);
            
            x = sin(Theta)*cos(Phi);
            y = sin(Theta)*sin(Phi);
            z = cos(Theta);
            
            ##############################################
            # project from all monitors, using symmetry/anti-symmetry where appropriate
            if (havedata("z2","Ex")) { # have z data, 3D simulation
            ##############################################
            # Angular distribution calculation for a 3D simulation begins
            
                npts = length(Phi);
            
                ######## x-y plane (phi=0 corresponds to the direction (0,1,0))
                ?"    Projecting";
            
                # Calculate far field by summing contribution from each monitor
                temp = farfieldexact("x2",x,y,z,i_target) + farfieldexact("y2",x,y,z,i_target) + farfieldexact("z2",x,y,z,i_target);
                shape = [npts,3,nf]; # deal with extra dimension not included if nf=1    
                temp = reshape(temp,shape);  
                if(havedata("x1")){
                  temp = temp - reshape(farfieldexact("x1",x,y,z,i_target),shape);
                }else{
                  temp2 = reshape(farfieldexact("x2",-x,y,z,i_target),shape);
                  s = symm_x*[1,-1,-1];
                  temp2(1:npts,1,:) = s(1)*temp2(1:npts,1,:);
                  temp2(1:npts,2,:) = s(2)*temp2(1:npts,2,:);
                  temp2(1:npts,3,:) = s(3)*temp2(1:npts,3,:);
                  temp = temp - temp2;
                }
                if(havedata("y1")){
                  temp = temp - reshape(farfieldexact("y1",x,y,z,i_target),shape);
                }else{
                  temp2 = reshape(farfieldexact("y2",x,-y,z,i_target),shape);
                  s = symm_y*[-1,1,-1];  	
                  temp2(1:npts,1,:) = s(1)*temp2(1:npts,1,:);
                  temp2(1:npts,2,:) = s(2)*temp2(1:npts,2,:);
                  temp2(1:npts,3,:) = s(3)*temp2(1:npts,3,:);
                  temp = temp - temp2;
                }
                if(include_z_min_boundary) {
                    if(havedata("z1")){
                      temp = temp - reshape(farfieldexact("z1",x,y,z,i_target),shape);
                    }else{
                        temp2 = reshape(farfieldexact("z2",x,y,-z,i_target),shape);
                        s = symm_z*[-1,-1,1];
                        temp2(1:npts,1,:) = s(1)*temp2(1:npts,1,:);
                        temp2(1:npts,2,:) = s(2)*temp2(1:npts,2,:);
                        temp2(1:npts,3,:) = s(3)*temp2(1:npts,3,:);
                        temp = temp - temp2;
                    }
                }
                Ex = temp(1:length(Phi),1,:);
                Ey = temp(1:length(Phi),2,:);
                Ez = temp(1:length(Phi),3,:);
                E2temp = reshape( abs(Ex)^2+abs(Ey)^2+abs(Ez)^2, [nt, np, nf] );
                
                # reconstruct other quadrants
                phi = linspace(0,360,2*np-1)*pi/180;
                E2 = matrix(nt,2*np-1,nf);
                E2(:,1:np,:) = E2temp;
                E2(:,(np+1):(2*np-1),:) = E2temp(:,(np-1):-1:1,:);
                #image(theta,phi,E2);
                index = getresult("index","index_x");
                index = real(pinch(index(1,1,1,i_target)));
                # calculate Poynting vector, normalize to sourcepower
                P = E2 * sqrt(eps0/mu0) * meshgrid3dz(theta,phi,index/sourcepower(f));
                
                Theta = meshgridx(theta,phi);
                Phi = meshgridy(theta,phi);
                np = length(phi);
                x = sin(Theta)*cos(Phi);
                y = sin(Theta)*sin(Phi);
                z = cos(Theta);
                
            } else {
                ?"Error, projection not available for 2D simulations";
                break;
            }
            
            ##############################################
            # construct the connectivity matrix for an unstructured mesh
            pos = find( (Theta != theta(end)) & (Phi != phi(end)) );
            C = matrix(2*(nt-1)*(np-!wrap_phi),3);
            cPos = 1:length(pos);
            C(cPos,1) = pos;
            C(cPos,2) = pos+1;
            C(cPos,3) = pos+1+nt;
            if(wrap_phi) {
                posEnd = find( (Theta != theta(end)) & (Phi == phi(end)) );
                cPosEnd = cPos(end) + (1:length(posEnd));
                C(cPosEnd,1) = posEnd;
                C(cPosEnd,2) = posEnd+1;
                C(cPosEnd,3) = posEnd+1+nt-(np*nt);
            } 
            cPos = cPos + length(pos) + wrap_phi*(nt-1);
            C(cPos,1) = pos;
            C(cPos,2) = pos+1+nt;
            C(cPos,3) = pos+nt;
            if(wrap_phi) {
                cPosEnd = cPos(end) + (1:length(posEnd));
                posEnd;
                C(cPosEnd,1) = posEnd;
                C(cPosEnd,2) = posEnd+1+nt-(np*nt);
                C(cPosEnd,3) = posEnd+nt-(np*nt);
            } 
            
            # construct a dataset and visualize it
            farfield = unstructureddataset("ff",x,y,z,C);
            farfield.addparameter("lambda",c/f,"f",f);
            farfield.addattribute("P",reshape(P,[nt*np,nf]));
                
            # integrate over all phi for symmetric result.
            P_vs_theta = matrixdataset("P");
            P_vs_theta.addparameter("theta_degrees",theta*180/pi,"theta_radians",theta);
            P_vs_theta.addparameter("lambda",c/f,"f",f);
            P_vs_theta.addattribute("P",integrate(P,2,phi)/(2*pi));
            
            # add the Purcell factor as a result for the same frequency grid, using dipolepower script command
            Purcell = matrixdataset("Purcell");
            Purcell.addparameter("lambda",c/f,"f",f);
            Purcell.addattribute("purcell",dipolepower(f)/sourcepower(f));
            

        """
    
        self.fdtd.set('analysis script', script)
        self.fdtd.addanalysisresult('T')
        self.fdtd.addanalysisresult('Purcell')
        self.fdtd.addanalysisresult('P_vs_theta')
        self.fdtd.addanalysisprop("near field points per wavelength", 0, 3)
        self.fdtd.addanalysisprop("N phi per half", 0, Nphi)
        self.fdtd.addanalysisprop("project all wavelengths", 0, 1)
        self.fdtd.addanalysisprop('include_z_min_boundary', 0, 1);
        self.fdtd.set('theta max', theta_max)
        self.fdtd.set('N theta', N_theta)
        self.fdtd.unselectall()

    def Qfactor(self, x=0, y=0, z=0):
        
        self.fdtd.unselectall()
        self.fdtd.select('source')
        fstart = self.fdtd.get('frequency start')
        fstop = self.fdtd.get('frequency stop')
        
        self.fdtd.addobject('Qanalysis')
        self.fdtd.select('Qanalysis')
        self.fdtd.set('x span', 0.1e-06)
        self.fdtd.set('y span', 0.1e-06)
        self.fdtd.set('f min', fstart)
        self.fdtd.set('f max', fstop)
        self.fdtd.set('t start', 100e-15)
        self.fdtd.set('x', x)
        self.fdtd.set('y', y)
        self.fdtd.set('z', z)
        self.fdtd.addanalysisresult('f_spectrum')
        self.fdtd.unselectall()

    def mode_volume_3D(self, xy_span=1e-06, z_span=1e-06, x=0, y=0, z=0):
        self.fdtd.addobject('mode_volume')
        self.fdtd.select('mode_volume')
        self.fdtd.set('x span', xy_span)
        self.fdtd.set('y span', xy_span)
        self.fdtd.set('z span', z_span)
        self.fdtd.set('x', x)
        self.fdtd.set('y', y)
        self.fdtd.set('z', z)
        self.fdtd.unselectall
        
    def mode_volume_2D(self, xy_span_pml=1e-06, apod="Start", apod_center=0, apod_start_w=100e-09,
                       z_min=0, z_max=1e-06, dipole_shift=0, lambda_res=640e-09):
        # self.fdtd.unselectall
        self.fdtd.addanalysisgroup()
        self.fdtd.set('name', 'mode volume 2D')
        
        self.mon.index_monitor(name="n", monitor_type="2D Y-normal",
                        x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                        z_min=z_min, z_max=z_max)
        
        self.fdtd.select('n')
        self.fdtd.addtogroup('mode volume 2D')
        
        self.mon.power_monitor(name='xy_middle', montype="2D Z-normal", plane="xy",
                              x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                              y_min=-.5*xy_span_pml, y_max=.5*xy_span_pml, 
                              z=0,
                              apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)        
        
        self.mon.power_monitor(name='xz_middle', montype="2D Y-normal", plane="xz",
                              x_min=-.5*xy_span_pml, x_max=.5*xy_span_pml, 
                              z_min=z_min, z_max=z_max, 
                              y=0,
                              apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)    
        
        self.mon.power_monitor(name='yz_middle', montype="2D X-normal", plane="yz",
                              y_min=-.5*xy_span_pml, y_max=.5*xy_span_pml, 
                              z_min=z_min, z_max=z_max, 
                              x=0,
                              apod="Start", apod_center=apod_center, apod_time_width=apod_start_w)    
        
        self.fdtd.select('xy_middle')
        self.fdtd.addtogroup('mode volume 2D')        
        
        self.fdtd.select('xz_middle')
        self.fdtd.addtogroup('mode volume 2D')        
        
        self.fdtd.select('yz_middle')
        self.fdtd.addtogroup('mode volume 2D')  
        
        # with open('./modeVol_2D.txt', 'r', encoding="utf8") as f:
        #     script = f.read().rstrip()
        
        script = """
            
            dipole_shift = 0; 
            
            x = getdata('xy_middle',"x");    
            y = getdata('xy_middle',"y");
            z = getdata('xz_middle',"z");    
            
            x_pt=size(x);           # Points in the axis
            z_pt=size(z);           # Points in the z axis            
            
            f = getresult('xy_middle', 'f');
            wlen = c/f;
            midpoint_z_ind = round(.5*size(x , 1));    # find the z midpoint index value
            midpoint_x_ind = floor(x_pt(1)/2)+1;
            
            
            
            # Integration ranges (0 -> rmax, zmin -> zmax)
            z_int_range = z;                                               
            r_int_range = x(1:midpoint_x_ind);                                      
            
            # Get electric fields
            E_xz = pinch(getelectric('xz_middle'));
            E_yz = pinch(getelectric('yz_middle'));
            
            # Get refractive indices
            n_xz = pinch(getdata('n',"index_z"));                       
            
            nv_zpos_ind = find(z, dipole_shift);            # Position of the NV
            n_nv_xz = real(n_xz(midpoint_x_ind, nv_zpos_ind)); # Refractive index at the position of the NV
            n_nv_yz = real(n_xz(midpoint_x_ind, nv_zpos_ind));                                   
            
            # Predefine arrays for looping
            Vol_abs_xz = ones(length(f));
            Vol_lam_xz = ones(length(f));
            
            Vol_abs_yz = ones(length(f));
            Vol_lam_yz = ones(length(f));
            
            # Loop over each frequency
            for(i=1:length(f)){
                n = n_xz(1:midpoint_x_ind,1:z_pt(1));
            
                E_xz_res = pinch(E_xz, 3, 1);
                E_xz_mid = E_xz_res(1:midpoint_x_ind,1:z_pt(1));                      # half of the x-z cut
                eps_E_xz = real(n^2)*E_xz_mid;
                # Calculate epsilon * E at the NV position 
                eps_E_xz_at_nv = eps_E_xz(midpoint_x_ind, nv_zpos_ind);                          
                eps_E_xz_max = max(eps_E_xz);    # Max value                                
                
                E_yz_res = pinch(E_yz, 3, 1);
                E_yz_mid = E_yz_res(1:midpoint_x_ind,1:z_pt(1));                      # half of the x-z cut
                eps_E_yz = real(n^2)*E_yz_mid;
                eps_E_yz_at_nv = eps_E_yz(midpoint_x_ind, nv_zpos_ind);                          
                eps_E_yz_max = max(eps_E_yz); 
                    
                
                V0_xz = E_xz_mid*real(n^2)*2*pi;                                  # n^2 
                Vol_raw1_xz = integrate(V0_xz,2,z_int_range);
                Vol_raw2_xz = 1e18*abs(integrate(Vol_raw1_xz*r_int_range,1,r_int_range)); # Unit in um^3.
                
                Vol_abs_xz(i) = Vol_raw2_xz/(eps_E_xz_max);                         # unit:um^3, for air-like mode
                Vol_lam_xz(i) = Vol_abs_xz(i)/(wlen(i)^3*1e18);   
                
                V0_yz = E_yz_mid*real(n^2)*2*pi;                                  # n^2 
                Vol_raw1_yz = integrate(V0_yz,2,z_int_range);
                Vol_raw2_yz = 1e18*abs(integrate(Vol_raw1_yz*r_int_range,1,r_int_range)); # Unit in um^3.
                
                Vol_abs_yz(i) = Vol_raw2_yz/(eps_E_yz_max);                         # unit:um^3, for air-like mode
                Vol_lam_yz(i) = Vol_abs_yz(i)/(wlen(i)^3*1e18);   
            }
            
            ## Average mode volume
            Vol_abs_avg = (Vol_abs_xz+Vol_abs_yz)/2;
            Vol_lam_avg = (Vol_lam_xz+Vol_lam_yz)/2;
            n_nv_avg = (n_nv_xz+n_nv_yz)/2;
            In2_nv_avg = (eps_E_xz_at_nv+eps_E_yz_at_nv)/2;
            In2_max_avg = (eps_E_xz_max+eps_E_yz_max)/2;             # unit:lam^3
        """
        
        self.fdtd.select('mode volume 2D')
        self.fdtd.addanalysisprop('dipole_shift', 0, dipole_shift)
        self.fdtd.addanalysisprop('lambda_res', 0, lambda_res)
        self.fdtd.addanalysisprop('XY_span', 0, 5e-06) # CHANGE THIS
        self.fdtd.addanalysisprop('res_i', 0, 10)
        
        self.fdtd.set('analysis script', script)
        self.fdtd.addanalysisresult('Vol_abs_xz')
        self.fdtd.addanalysisresult('Vol_lam_xz')
        self.fdtd.addanalysisresult('Vol_abs_yz')
        self.fdtd.addanalysisresult('Vol_lam_yz')
        self.fdtd.addanalysisresult('Vol_abs_avg')
        self.fdtd.addanalysisresult('Vol_lam_avg')
        
        
        self.fdtd.unselectall()
    
