# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 16:45:24 2025

@author: ychuang
"""

import matplotlib.pyplot as plt
from scipy import interpolate
from twscan_module.twinscan_prepare import twscan_assist


class midnd_plot:
    
    def __init__(self, DF, data, twa: twscan_assist):
        
        self.DF = DF
        self.data = data
        self.twa = twa
    
    
    def neudenplot_method(self, iterlist, cl_dic, A_dic, scan_style, scandetail, xcoord_type):
        
        
        midplane_psi = self.data['midplane_calc']['psi_solps_mid']
        r_rsep = self.data['midplane_calc']['R_Rsep']
        
        
        psi_to_dsa_func = interpolate.interp1d(midplane_psi, r_rsep, fill_value = 'extrapolate')
        
       
        fig, axs = plt.subplots()
        
        
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        dat_struc = {'nx': nx, 'ny': ny}

        
        for aa in iterlist:
            
            
            """
            label= 'core density {} $10^{19}$'.format(aa)
            
            """
            
            axs.legend(loc= 'lower left', fontsize=10)
            
            if self.DF.series_flag == 'twin_scan':
                
                if scan_style == 'tempscan':
                    
                    ad = aa[1]
                
                elif scan_style == 'denscan':
                    
                    ad = aa[0]
                
                else:
                    print('neteTSplot_method, please check scan_style')
                    
                
                psi_coord = self.data['midplane_profile'][aa[0]][aa[1]]['psiN']
                mid_neu_pro = self.data['midplane_profile'][aa[0]][aa[1]]['mid_nd']
                rrsep_solps = psi_to_dsa_func(psi_coord)
            
            else:
                ad = aa
                
                psi_coord = self.data['midplane_profile'][aa[0]][aa[1]]['psiN']
                mid_neu_pro = self.data['midplane_profile'][aa[0]][aa[1]]['mid_nd']
                rrsep_solps = psi_to_dsa_func(psi_coord)
                
            

            if scan_style == 'denscan':
                
                
                if xcoord_type == 'psi':
                    
                    axs.plot(psi_coord, mid_neu_pro, color = cl_dic[ad], 
                             label= '{}'.format(A_dic[ad]))
                
                elif xcoord_type == 'rrsep':
                    
                    axs.plot(rrsep_solps, mid_neu_pro, color = cl_dic[ad], 
                             label= '{}'.format(A_dic[ad]))
                    
                    
                axs.set_title('Density scan with Te = {} eV'.format(scandetail))
                axs.legend()
            
            elif scan_style == 'tempscan':
                
                
                if xcoord_type == 'psi':
                    
                    axs.plot(psi_coord, mid_neu_pro, color = cl_dic[ad], 
                                label= '{}'.format(A_dic[ad]))
                
                elif xcoord_type == 'rrsep':
                    
                    axs.plot(rrsep_solps, mid_neu_pro, color = cl_dic[ad], 
                                label= '{}'.format(A_dic[ad]))
                    
                    
                    
                axs.set_title('Temperature scan with Ne = {}'.format(scandetail))
                axs.legend()

                
            
            else:
                print('neudenplot_method, please check the scan parameter')
    
    
    
    def midnd_plot(self, scan_style, xcoord_type):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == True and withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            label_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            asp_ch = self.data['dircomp']['multi_shift']
            
            print('I am still working on it')
        
        elif withshift == False and withseries == True:
            
            # series_flag = self.DefaultSettings['series_flag']
            
            
            if self.DF.series_flag == 'twin_scan':
                
                dircomp = self.data['dircomp']
                
                if scan_style == 'tempscan':
                    
                    key_a = 'denscan_list'
                    key_b = 'tempscan_list'
                
                elif scan_style == 'denscan':
                    
                    key_a = 'tempscan_list'
                    key_b = 'denscan_list'
                
                else:
                    print('twinscan_plot_method, please check the scan_style!')
                
                keylist_a = []
                
                
                for x in dircomp[key_a]:
                    keylist_a.append('{:.3f}'.format(x))
                
                for ta in keylist_a:
                    
                    keylist_b = []
                    
                    for x in dircomp[key_b]:
                        keylist_b.append('{:.3f}'.format(x))
                    
                    
                    iter_key, color_dic, scan_title, label_dic = self.twa.twinscan_prep(ta = ta, 
                    keylist_b = keylist_b, scan_style = scan_style)
                    
                    
                    print('check:')
                    print(iter_key)
                    print(color_dic)
                    print(label_dic)
                    
                    self.neudenplot_method(iterlist = iter_key, cl_dic = color_dic, 
                                A_dic = label_dic, scan_style = scan_style, xcoord_type= xcoord_type,
                                scandetail = scan_title)
                    
             
            else:
                print('neteTS_plot, please check the series flag')