# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 15:30:16 2024

@author: ychuang
"""

import matplotlib.pyplot as plt
from SOLPSplotter_ioutflux import iout_flux
from SOLPSplotter_NTplot import NT_plot
import numpy as np
from matplotlib.offsetbox import AnchoredText


class twinscan_ioutflux(iout_flux, NT_plot):
    def __init__(self, DefaultSettings, loadDS):
        NT_plot.__init__(self, DefaultSettings, loadDS)
        iout_flux.__init__(self, DefaultSettings, loadDS)
    
    
    def set_plot(self):
        if self.Publish == 'b2plottersetting':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 3, markersize= 5)
            plt.rcParams.update({'font.size': 14})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
            # plt.rcParams["text.usetex"] = True
  
        else:
            print('Publish setting is incorrect or add another setting')


    def twinscaniout_method(self, iterlist, cl_dic, scan_style, ang_list, plot_option, format_option):
        
        if format_option == '4x1':
            
            fig, axs = plt.subplots(4, 1)
        
        elif format_option == '2x2':
            
            fig, axs = plt.subplots(2, 2)
        
              
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
                
        # if dat_size == 'full':

        #     dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        # elif dat_size == 'small':
        #     dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        plt.subplots_adjust(hspace=.0)
        anchored_text_1 = AnchoredText('{}'.format('(a) Core particle flux [$s^{-1}*m^{-2}$]'), 
                                              loc='lower center')
        anchored_text_2 = AnchoredText('{}'.format('(b) Separatrix particle flux [$s^{-1}*m^{-2}$]'), 
                                              loc='lower center')
        anchored_text_3 = AnchoredText('{}'.format('(a) Core heat flux [$s^{-1}*m^{-2}$]'), 
                                              loc='upper right')
        anchored_text_4 = AnchoredText('{}'.format('(b) Separatrix heat flux [$s^{-1}*m^{-2}$]'), 
                                             loc='lower center')
        
        
        if scan_style == 'denscan':
            
            if format_option == '4x1':
                
                axs[3].axvline(x= 0, color='black', lw=3, ls='-', label= 'LFS')
                axs[3].axvline(x= 180, color='brown', lw=3, ls='-', label= 'HFS')
                axs[3].legend(loc = 'upper right')
        
        
                 
        
        for aa in iterlist:
            
            # psi_coord, mid_ne_pro, mid_te_pro, mid_neu_pro = self.nete_midprof(itername = aa, 
            #                                         data_struc = dat_struc)
            if plot_option == 'radial flux poloidal plot':
                fna_sol = self.data['iout_data']['flux_y_0'][aa[0]][aa[1]][19, 25:73]
                fna_core = self.data['iout_data']['flux_y_0'][aa[0]][aa[1]][0, 25:73]
                fhe_sol = self.data['iout_data']['eheatflux_y_0'][aa[0]][aa[1]][19, 25:73]
                fhe_core = self.data['iout_data']['eheatflux_y_0'][aa[0]][aa[1]][0, 25:73]
                

            """
            label= 'core density {} $10^{19}$'.format(aa)
            
            """
            
            # axs.legend(loc= 'lower left', fontsize=10)
            
            if self.series_flag == 'twin_scan':
                
                if scan_style == 'tempscan':
                    
                    ad = aa[1]
                    ap = aa[0]
                    
                
                elif scan_style == 'denscan':
                    
                    ad = aa[0]
                    ap = aa[1]
                
                else:
                    print('neteTSplot_method, please check scan_style')
            
            else:
                ad = aa
                
            

            if scan_style == 'denscan':
                
                title_ap = float(ap)*pow(10, 5)
                label_ad = float(ad)*pow(10, 20)
                
                if format_option == '4x1':
                    
                    axs[0].plot(ang_list, fna_core, color = cl_dic[ad], 
                             label= '{:.3E} (1/s)'.format(label_ad))
                    axs[1].plot(ang_list, fna_sol, color = cl_dic[ad], 
                             label= '{:.3E} (1/s)'.format(label_ad))
                    axs[2].plot(ang_list, fhe_core, color = cl_dic[ad], 
                             label= '{:.3E} (1/s)'.format(label_ad))
                    axs[3].plot(ang_list, fhe_sol, color = cl_dic[ad], 
                             label= '{:.3E} (1/s)'.format(label_ad))
                    axs[2].legend(loc = 'lower left')

                
                elif format_option == '2x2':
                    
                    axs[0, 0].plot(ang_list, fna_core, color = cl_dic[ad], 
                             label= '{:.3E} (1/s)'.format(label_ad))
                    axs[1, 0].plot(ang_list, fna_sol, color = cl_dic[ad], 
                             label= '{:.3E} (1/s)'.format(label_ad))
                    # axs[1, 0].set_yscale('log')
                    axs[0, 1].plot(ang_list, fhe_core, color = cl_dic[ad], 
                             label= '{:.3E} (1/s)'.format(label_ad))
                    axs[1, 1].plot(ang_list, fhe_sol, color = cl_dic[ad], 
                             label= '{:.3E} (1/s)'.format(label_ad))
                    
                
            elif scan_style == 'tempscan':
                
                title_ap = float(ap)*pow(10, 20)
                label_ad = float(ad)*pow(10, 5)
                
                if format_option == '4x1':
                    axs[0].plot(ang_list, fna_core, color = cl_dic[ad], 
                             label= '{:.3E} W'.format(label_ad))
                    axs[1].plot(ang_list, fna_sol, color = cl_dic[ad], 
                             label= '{:.3E} W'.format(label_ad))
                    axs[2].plot(ang_list, fhe_core, color = cl_dic[ad], 
                             label= '{:.3E} W'.format(label_ad))
                    axs[3].plot(ang_list, fhe_sol, color = cl_dic[ad], 
                             label= '{:.3E} W'.format(label_ad))


                elif format_option == '2x2':
                    
                    axs[0, 0].plot(ang_list, fna_core, color = cl_dic[ad], 
                             label= '{:.3E} W'.format(label_ad))
                    axs[1, 0].plot(ang_list, fna_sol, color = cl_dic[ad], 
                             label= '{:.3E} W'.format(label_ad))
                    # axs[1, 0].set_yscale('log')
                    axs[0, 1].plot(ang_list, fhe_core, color = cl_dic[ad], 
                             label= '{:.3E} W'.format(label_ad))
                    axs[1, 1].plot(ang_list, fhe_sol, color = cl_dic[ad], 
                             label= '{:.3E} W'.format(label_ad))
            
            else:
                print('neudenplot_method, please check the scan parameter')

                
        if scan_style == 'denscan':
            
            title_ap = float(ap)*pow(10, 5)
            label_ad = float(ad)*pow(10, 20)
            
            if format_option == '4x1':
                
                axs[0].set_title('Particle flux scan with heat flux = {:.3E} W'.format(title_ap))
                axs[0].add_artist(anchored_text_1)
                axs[1].add_artist(anchored_text_2)
                axs[2].add_artist(anchored_text_3)
                axs[3].add_artist(anchored_text_4)
                axs[0].axvline(x= 0, color='black', lw=3, ls='-')
                axs[1].axvline(x= 0, color='black', lw=3, ls='-')
                axs[2].axvline(x= 0, color='black', lw=3, ls='-')
                axs[3].axvline(x= 0, color='black', lw=3, ls='-', label= 'LFS')
                axs[0].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[1].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[2].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[3].axvline(x= 180, color='brown', lw=3, ls='-', label= 'HFS')
                axs[3].set_xlabel('poloidal angle')
                axs[2].legend(loc = 'lower left')
                # axs[3].legend(loc = 'upper right')
            
            elif format_option == '2x2':
                
                fig.suptitle('Particle flux scan with heat flux = {:.3E} W'.format(title_ap))
                axs[0, 0].add_artist(anchored_text_1)
                axs[1, 0].add_artist(anchored_text_2)
                axs[0, 1].add_artist(anchored_text_3)
                axs[1, 1].add_artist(anchored_text_4)
                axs[0, 0].axvline(x= 0, color='black', lw=3, ls='-', label= 'LFS')
                axs[1, 0].axvline(x= 0, color='black', lw=3, ls='-')
                axs[0, 1].axvline(x= 0, color='black', lw=3, ls='-')
                axs[1, 1].axvline(x= 0, color='black', lw=3, ls='-')
                axs[0, 0].axvline(x= 180, color='brown', lw=3, ls='-', label= 'HFS')
                axs[1, 0].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[0, 1].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[1, 1].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[1, 0].set_xlabel('poloidal angle')
                axs[1, 1].set_xlabel('poloidal angle')
                axs[0, 0].legend(loc = 'best')
        
        
        elif scan_style == 'tempscan':
            
            title_ap = float(ap)*pow(10, 20)
            label_ad = float(ad)*pow(10, 5)
            
            if format_option == '4x1':
                                
                axs[0].set_title('Heat flux scan with particle flux = {:.3E} (1/s)'.format(title_ap))
                axs[0].add_artist(anchored_text_1)
                axs[1].add_artist(anchored_text_2)
                axs[2].add_artist(anchored_text_3)
                axs[3].add_artist(anchored_text_4)
                axs[0].axvline(x= 0, color='black', lw=3, ls='-', label= 'LFS')
                axs[1].axvline(x= 0, color='black', lw=3, ls='-')
                axs[2].axvline(x= 0, color='black', lw=3, ls='-')
                axs[3].axvline(x= 0, color='black', lw=3, ls='-')
                axs[0].axvline(x= 180, color='brown', lw=3, ls='-', label= 'HFS')
                axs[1].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[2].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[3].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[3].set_xlabel('poloidal angle')
                axs[0].legend(loc = 'best')

            elif format_option == '2x2':
                
                fig.suptitle('Heat flux scan with particle flux = {:.3E} (1/s)'.format(title_ap))
                axs[0, 0].add_artist(anchored_text_1)
                axs[1, 0].add_artist(anchored_text_2)
                axs[0, 1].add_artist(anchored_text_3)
                axs[1, 1].add_artist(anchored_text_4)
                axs[0, 0].axvline(x= 0, color='black', lw=3, ls='-', label= 'LFS')
                axs[1, 0].axvline(x= 0, color='black', lw=3, ls='-')
                axs[0, 1].axvline(x= 0, color='black', lw=3, ls='-')
                axs[1, 1].axvline(x= 0, color='black', lw=3, ls='-')
                axs[0, 0].axvline(x= 180, color='brown', lw=3, ls='-', label = 'HFS')
                axs[1, 0].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[0, 1].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[1, 1].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[1, 0].set_xlabel('poloidal angle')
                axs[1, 1].set_xlabel('poloidal angle')
                axs[0, 0].legend(loc = 'best')
       
                
                
        

       


    def iout_twinscan_prep(self, ta, keylist_b, scan_style):
     
         if self.withshift == False and self.withseries == True:
             
             if self.series_flag == 'twin_scan':
                 
                 color_list = ['red', 'orange', 'green', 'blue', 'purple']
                 
                 color_dic = self.pair_dic(keys = keylist_b, values = color_list)
                 
                 
                 scan_list = []
                 iter_key = []
                 
                 
                 for tb in keylist_b:
                     
                     if scan_style == 'tempscan':
                         
                         it_in = (ta, tb)
                     
                     elif scan_style == 'denscan':
                         
                         it_in = (tb, ta)
                     
                     else:
                         print('twinscan_plot_method, please check the scan_style!')
                     
                     
                     nx = self.data['b2fgeo']['nx']
                     ny = self.data['b2fgeo']['ny']
                     
                     iter_key.append(it_in)
             
                 return iter_key, color_dic


    def twinscan_ioutplot(self, scan_style, plot_option, format_option):
        
        if self.withshift == False and self.withseries == True:
            
            # series_flag = self.DefaultSettings['series_flag']
            
            
            if self.series_flag == 'twin_scan':
                
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
                
                pol_list_a = []
                for i in range(48):
                    pol_list_a.append('{}'.format(24 + i))
                
                
                self.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
                ang_list = self.data['angle']['angle_list']
                
                for x in dircomp[key_a]:
                    keylist_a.append('{:.3f}'.format(x))
                
                for ta in keylist_a:
                    
                    keylist_b = []
                    
                    for x in dircomp[key_b]:
                        keylist_b.append('{:.3f}'.format(x))
                    
                    
                    iter_key, color_dic= self.iout_twinscan_prep(ta = ta, 
                    keylist_b = keylist_b, scan_style = scan_style)
                    
                    
                    print('check:')
                    print(iter_key)
                    print(color_dic)
                    # print(label_dic)
                    self.twinscaniout_method(iterlist = iter_key, cl_dic = color_dic, 
                                scan_style = scan_style, ang_list = ang_list, 
                    plot_option = plot_option, format_option = format_option)
                    
             
            else:
                print('neteTS_plot, please check the series flag')