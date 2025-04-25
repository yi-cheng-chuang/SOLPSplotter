# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 00:37:51 2025

@author: ychuang
"""

"""

    #calculate the poloidal angle function    
            
    def rad_polangle_method(self, RadLoc, VertLoc, gfile, shift, pol_list, plot_angle, rn):
            
            
            mag_axis_z = gfile['zmaxis']
            # print(mag_axis_z)
            
            mag_axis_r = gfile['rmaxis'] + shift
            
                             
            xprloc = RadLoc[:, 72][19]
            xpzloc = VertLoc[:, 72][19]
            
            xpoint_R = xprloc - mag_axis_r
            xpoint_Z = xpzloc - mag_axis_z
            xpoint_angle = np.arctan2(xpoint_Z, xpoint_R)*180 / np.pi
            

            angle_list = []
            
            for pol_index in pol_list:
                rloc = RadLoc[:, int(pol_index)][rn]
                zloc = VertLoc[:, int(pol_index)][rn]
                
                x = rloc - mag_axis_r
                y = zloc - mag_axis_z
                
                
                angle = np.arctan2(y, x)*180 / np.pi
                if angle <= xpoint_angle:
                    angle = angle + 360
                angle_list.append(angle)
                
            pol_loc = np.asarray(list(map(int, pol_list)))
            
            if plot_angle: 
                plt.figure(figsize=(7,7))
                plt.plot(pol_loc, angle_list, 'o', color = 'r', label= 'poloidal angle')
                plt.xlabel('poloidal index')
                plt.title('poloidal angle verses poloidal index')
                plt.legend()
                    
            return angle_list, xpoint_angle
        
        
        
    def calc_polangle_rad(self, pol_list, plot_angle, input_rloc):
            
            withshift = self.DF.withshift
            withseries = self.DF.withseries
        
        
            if withshift == False and withseries == False:
                
                rad_grid = self.data['grid']['RadLoc']
                vert_grid = self.data['grid']['VertLoc']     
                g_data = self.data['gfile']['g']
                shift_val = self.data['dircomp']['shift_value']
                
                            
                angle_list, xpoint_angle = self.rad_polangle_method(RadLoc = rad_grid, VertLoc = vert_grid, 
                    gfile = g_data, shift = shift_val, plot_angle = plot_angle, 
                    pol_list = pol_list, rn = input_rloc)
                
                print('xpoint angle is {:.2f}'.format(xpoint_angle))
                
                angle_dic = {'angle_list': angle_list, 'xpoint_angle': xpoint_angle}
                
                self.data['angle'] = angle_dic
                
            elif withshift == True and withseries == False:        
                angle_dic = {}
                angle_list_dic = {}
                xpoint_angle_dic = {}
                index_angle_dic = {}
                for aa in self.data['dircomp']['multi_shift']:
                    
                    rad_grid = self.data['grid']['RadLoc'][aa]
                    vert_grid = self.data['grid']['VertLoc'][aa]    
                    g_data = self.data['gfile']['g']
                    shift_val = self.data['dircomp']['shift_dic'][aa]
                    
                    angle_list, xpoint_angle = self.rad_polangle_method(RadLoc = rad_grid, VertLoc = vert_grid, 
                        gfile = g_data, shift = shift_val, plot_angle = plot_angle, 
                        pol_list = pol_list, rn = input_rloc)
                    
                    angle_list_dic[aa] = angle_list
                    xpoint_angle_dic[aa] = xpoint_angle
                    
                    
                    print('xpoint angle is {:.2f} for {} case'.format(xpoint_angle, aa))
                    
                       
                angle_dic = {'angle_list': angle_list_dic, 'xpoint_angle': xpoint_angle_dic}
                
                self.data['angle'] = angle_dic
                
               
            elif withshift == False and withseries == True:
                
                den_list = list(self.data['dircomp']['Attempt'].keys())
                aa = den_list[0]
                
                rad_grid = self.data['grid']['RadLoc'][aa]
                vert_grid = self.data['grid']['VertLoc'][aa]    
                g_data = self.data['gfile']['g']
                shift_val = self.data['dircomp']['shift_value']
                
                
                angle_list, xpoint_angle = self.calc_pol_angle_method(RadLoc = rad_grid, VertLoc = vert_grid, 
                    gfile = g_data, shift = shift_val, plot_angle = plot_angle, pol_list = pol_list)
                
                angle_dic = {'angle_list': angle_list, 'xpoint_angle': xpoint_angle}
                
                self.data['angle'] = angle_dic
            
            elif withshift == True and withseries == True:
                print('calc_pol_angle function is not there yet!')
            
            else:
                
                             
                
                print('calc_pol_angle function has a bug')



"""



"""


Old RZ coordinate function load with OutputGen command


#-----------------------------------------------------------------------------
    

    def calcpsi_method(self, rad_range, pol_range, psi_RBS, attempt_number, outputdir):
               
        psival = np.zeros((rad_range, pol_range))
       
        RadLoc = np.loadtxt('{}/RadLoc{}'.format(outputdir, str(attempt_number)),
                    usecols = (3)).reshape((rad_range, pol_range))
        VertLoc = np.loadtxt('{}/VertLoc{}'.format(outputdir, str(attempt_number)), 
                      usecols = (3)).reshape((rad_range, pol_range))
        
        # coord_dic = {'RadLoc': RadLoc, 'VertLoc': VertLoc}
            
        for pol_loc in range(pol_range):
            for i in range(rad_range):
                # print(i)
                psival[i, pol_loc] = psi_RBS(RadLoc[i, pol_loc], 
                                                      VertLoc[i, pol_loc])
        return RadLoc, VertLoc, psival, pol_range, rad_range


                
    def calcpsi(self):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
            
        if withshift == False and withseries == False:
            n_pol = int(self.data['b2fgeo']['nx'] + 2)
            # print('xdim is {}'.format(pol_range))
            n_rad = int(self.data['b2fgeo']['ny'] + 2)
            # print('ydim is {}'.format(rad_range))
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            # psiNinterp_RGI = self.data['gfile']['gcomp']['interp_dic']['RGI'] 
            # psiNinterp_2d = self.data['gfile']['gcomp']['interp_dic']['2d']
            Attempt = self.data['dircomp']['Attempt']
            DRT = self.data['dirdata']['outputdir']['Output']
            

            RadLoc, VertLoc, psival, pol_range, rad_range = self.calcpsi_method(rad_range = n_rad, 
                                pol_range = n_pol, psi_RBS = psiNinterp_RBS, 
                                attempt_number = Attempt, outputdir = DRT)
            
            
            coord_dic = {'RadLoc': RadLoc, 'VertLoc': VertLoc}
            self.data['grid'] = coord_dic
            self.data['psi']['psival'] = psival
            self.data['DefaultSettings']['XDIM'] = pol_range
            self.data['DefaultSettings']['YDIM'] = rad_range
        
        elif withshift == True and withseries == False:
            psival_dic = {}
            RadLoc_dic = {}
            VertLoc_dic = {}
            xdim_dic = {}
            ydim_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                
                n_pol = int(self.data['b2fgeo'][aa]['nx'] + 2)
                # print('xdim is {}'.format(str(pol_range)))
                n_rad = int(self.data['b2fgeo'][aa]['ny'] + 2)
                # print('ydim is {}'.format(str(rad_range)))
                psiNinterp_RBS = self.data['gfile']['gcomp'][aa]['interp_dic']['RBS']
                Attempt = self.data['dircomp']['Attempt'][aa]
                DRT = self.data['dirdata']['outputdir'][aa]['Output']
                
                RadLoc, VertLoc, psival, pol_range, rad_range = self.calcpsi_method(rad_range = n_rad, 
                                    pol_range = n_pol, psi_RBS = psiNinterp_RBS, 
                                    attempt_number = Attempt, outputdir = DRT)
                
                psival_dic[aa] = psival
                RadLoc_dic[aa] = RadLoc
                VertLoc_dic[aa] = VertLoc
                xdim_dic[aa] = pol_range
                ydim_dic[aa] = rad_range
            
            coord_dic = {'RadLoc': RadLoc_dic, 'VertLoc': VertLoc_dic}
            self.data['grid'] = coord_dic
            self.data['psi']['psival'] = psival_dic
            self.data['DefaultSettings']['XDIM'] = xdim_dic
            self.data['DefaultSettings']['YDIM'] = ydim_dic
                
            
        elif withshift == False and withseries == True:
            # print('we are working on calcpsi for series case')
            
            series_rep = list(self.data['dircomp']['Attempt'].keys())[0]
            n_pol = int(self.data['b2fgeo']['nx'] + 2)
            # print('xdim is {}'.format(pol_range))
            n_rad = int(self.data['b2fgeo']['ny'] + 2)
            # print('ydim is {}'.format(rad_range))
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            Attempt = self.data['dircomp']['Attempt'][series_rep]
            DRT = self.data['dirdata']['outputdir'][series_rep]['Output']
            
            
            RadLoc, VertLoc, psival, pol_range, rad_range = self.calcpsi_method(rad_range = n_rad, 
                                pol_range = n_pol, psi_RBS = psiNinterp_RBS, 
                                attempt_number = Attempt, outputdir = DRT)
            
            coord_dic = {'RadLoc': RadLoc_dic, 'VertLoc': VertLoc_dic}
            self.data['grid'] = coord_dic
            self.data['psi']['psival'] = psival       
            self.data['DefaultSettings']['XDIM'] = pol_range
            self.data['DefaultSettings']['YDIM'] = rad_range
        
        elif withshift == True and withseries == True:
            print('calcpsi function is not there yet!')
        
        else:
            print('calcpsi function has a bug!')
        
    

#----------------------------------------------------------------------------


    def calcpsi_block_method(self, file_loc, shift):
        b2mn = self.lcm.scrape_b2mn(file_loc + '/b2mn.dat')
        print(file_loc)
        print(b2mn)
        jxa = b2mn['jxa']
        pol_index = int(jxa)
        
        simu_dir = file_loc.rsplit("/",1)[0]
        fname = file_loc.rsplit("/",1)[1]
        print('plasmf input: {}'.format(simu_dir))
        print('filename input: {}'.format(fname))
        
        geo = self.lcm.read_b2fgmtry(simu_dir + '/baserun/b2fgmtry')
        rad_range = int(geo['ny'] + 2)
        
        crLowerLeft = geo['crx'][pol_index,:,0]
        crLowerRight = geo['crx'][pol_index,:,1]
        crUpperLeft = geo['crx'][pol_index,:,2]
        crUpperRight = geo['crx'][pol_index,:,3]
        czLowerLeft = geo['cry'][pol_index,:,0]
        czLowerRight = geo['cry'][pol_index,:,1]
        czUpperLeft = geo['cry'][pol_index,:,2]
        czUpperRight = geo['cry'][pol_index,:,3]
        
        
        g = self.lcm.loadg(self.data['dirdata']['gbase'] 
                               + '/MAST__RMP_results/g027205.00275_efitpp')
        # g = loadg(self.data['dirdata']['gdir'][0])
        psiN = (g['psirz'] - g['simag']) / (g['sibry'] - g['simag'])

        dR = g['rdim'] / (g['nw'] - 1)
        dZ = g['zdim'] / (g['nh'] - 1)
        
        gZ = np.zeros(g['nh'])
        for i in range(g['nh']):
            gZ[i] = g['zmid'] - 0.5 * g['zdim'] + i * dZ
            
        
        gR = np.zeros(g['nw'])
        for i in range(g['nw']):
            gR[i] = g['rleft'] + i * dR + float(shift)
            

        psiNinterp_RBS = interpolate.RectBivariateSpline(gR, gZ, np.transpose(psiN))
        
        
        psi_solps_RBS = np.zeros(rad_range)
        for i in range(rad_range):
            psi_LL_RBS = psiNinterp_RBS(crLowerLeft[i], czLowerLeft[i])
            psi_UL_RBS = psiNinterp_RBS(crUpperLeft[i], czUpperLeft[i])
            psi_LR_RBS = psiNinterp_RBS(crLowerRight[i], czLowerRight[i])
            psi_UR_RBS = psiNinterp_RBS(crUpperRight[i], czUpperRight[i])
            psi_solps_RBS[i] = np.mean([psi_LL_RBS, psi_UL_RBS, 
                                        psi_LR_RBS, psi_UR_RBS])
        
        return psi_solps_RBS



"""



