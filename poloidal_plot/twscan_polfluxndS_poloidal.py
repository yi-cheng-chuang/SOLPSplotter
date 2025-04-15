# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 01:37:50 2024

@author: ychuang
"""


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnchoredText
from twscan_module.twinscan_prepare import twscan_assist


class twscan_polfluxndS_polplot:
    
    
    def __init__(self, DF, data, twa: twscan_assist):
        
        self.DF = DF
        self.data = data
        self.twa = twa

    
    def twpolfluxndS_method(self, iterlist, cl_dic, scan_style, plot_option, ang_list, pol_list, log_flag):
        

            
        fig, axs = plt.subplots(3, 1)
        
        
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']

        
        plt.subplots_adjust(hspace=.0)
        anchored_text_1 = AnchoredText('{}'.format('Poloidal flux [$m^{-2}*s^{-1}$]'), 
                                              loc='lower center')
        anchored_text_2 = AnchoredText('{}'.format('Source [$m^{-3}*s^{-1}$]'), 
                                              loc='upper center')
        anchored_text_3 = AnchoredText('{}'.format('Neutral density [$m^{-3}$]'), 
                                              loc='upper right')
        
        for aa in iterlist:
            
            
            fnax = self.data['b2wdat'][aa[0]][aa[1]]['b2npc_fnaxs'][0][1:97, 1:37]
            hz = self.data['b2wdat'][aa[0]][aa[1]]['hz'][1:nx+1, 1:ny+1]
            hy = self.data['b2wdat'][aa[0]][aa[1]]['hy'][1:nx+1, 1:ny+1]
            tor_area = np.multiply(hz, hy)
            fnaxs = np.divide(fnax, tor_area)           
            
            
            s_term = self.data['b2wdat'][aa[0]][aa[1]]['b2npc_sna'][0][1:nx+1, 1:ny+1]
            vol = self.data['b2wdat'][aa[0]][aa[1]]['vol'][1:nx+1, 1:ny+1]
            source = np.divide(s_term, vol)
            
            neuden = self.data['ft44'][aa[0]][aa[1]]['dab2'][:, :, 0]
            
            

            # st = int(pol_list[0])
            # ed = int(pol_list[-1]) + 1
            
            
            
            if plot_option == 'core':
                
                fx_list = []
                s_list = []
                nd_list = []
                
                for ii in pol_list:
                    
                    fx_list.append(fnaxs[int(ii), 18])
                    s_list.append(source[int(ii), 18])
                    nd_list.append(neuden[int(ii), 18])
                
                
            elif plot_option == 'inner leg':
                
                index_list = np.linspace(0, 25, 26)
                
                fx_list = []
                s_list = []
                nd_list = []
                
                for ii in index_list:
                    
                    fx_list.append(abs(fnaxs[int(ii), 18]))
                    s_list.append(source[int(ii), 18])
                    nd_list.append(neuden[int(ii), 18])
                
                
            elif plot_option == 'outer leg':
                
                index_list = np.linspace(60, 95, 36)
                
                fx_list = []
                s_list = []
                nd_list = []
                
                for ii in index_list:
                    
                    fx_list.append(fnaxs[int(ii), 18])
                    s_list.append(source[int(ii), 18])
                    nd_list.append(neuden[int(ii), 18])
            
            
            
            
            
            if self.DF.series_flag == 'twin_scan':
                
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
                
                
                if log_flag:
                    axs[1].set_yscale('log')
                    axs[2].set_yscale('log')
                else:
                    pass
                    
                axs[0].set_title('Particle flux scan with heat flux = {:.3E} W'.format(title_ap))
                
                if plot_option == 'core':
                    axs[0].plot(ang_list, fx_list, color = cl_dic[ad], 
                              label= '{:.3E} (1/s)'.format(label_ad))
                    axs[1].plot(ang_list, s_list, color = cl_dic[ad])
                    axs[2].plot(ang_list, nd_list, color = cl_dic[ad])
                else:
                    axs[0].plot(index_list, fx_list, color = cl_dic[ad], 
                              label= '{:.3E} (1/s)'.format(label_ad))
                    axs[1].plot(index_list, s_list, color = cl_dic[ad])
                    axs[2].plot(index_list, nd_list, color = cl_dic[ad])
                    
                
                axs[2].set_xlabel('poloidal angle')
                axs[0].legend(loc = 'best')
                
                
            elif scan_style == 'tempscan':
                
                title_ap = float(ap)*pow(10, 20)
                label_ad = float(ad)*pow(10, 5)
                
                
                if log_flag:
                    axs[1].set_yscale('log')
                    axs[2].set_yscale('log')
                else:
                    pass
                
                
                
                axs[0].set_title('Heat flux scan with particle flux = {:.3E} (1/s)'.format(title_ap))
                
                if plot_option == 'core':
                    
                    axs[0].plot(ang_list, fx_list, color = cl_dic[ad], 
                              label= '{:.3E} W'.format(label_ad))
                    axs[1].plot(ang_list, s_list, color = cl_dic[ad])
                    axs[2].plot(ang_list, nd_list, color = cl_dic[ad])
                else:
                    axs[0].plot(index_list, fx_list, color = cl_dic[ad], 
                              label= '{:.3E} W'.format(label_ad))
                    axs[1].plot(index_list, s_list, color = cl_dic[ad])
                    axs[2].plot(index_list, nd_list, color = cl_dic[ad])
                    
                    
                axs[2].set_xlabel('poloidal angle')
                axs[0].legend(loc = 'best')
            
            else:
                print('neudenplot_method, please check the scan parameter')
        
        
        
        
        if scan_style == 'denscan':
            
            
            axs[0].add_artist(anchored_text_1)
            axs[1].add_artist(anchored_text_2)
            axs[2].add_artist(anchored_text_3)
            
            if plot_option == 'core':
                axs[0].axvline(x= 0, color='black', lw=3, ls='-')
                axs[1].axvline(x= 0, color='black', lw=3, ls='-', label= 'LFS')
                axs[2].axvline(x= 0, color='black', lw=3, ls='-')
                axs[0].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[1].axvline(x= 180, color='brown', lw=3, ls='-', label= 'HFS')
                axs[2].axvline(x= 180, color='brown', lw=3, ls='-')
            
            else:
                pass
            
            axs[1].legend(loc = 'best')
            
            
        elif scan_style == 'tempscan':

            axs[0].add_artist(anchored_text_1)
            axs[1].add_artist(anchored_text_2)
            axs[2].add_artist(anchored_text_3)
            if plot_option == 'core':
                
                axs[0].axvline(x= 0, color='black', lw=3, ls='-')
                axs[1].axvline(x= 0, color='black', lw=3, ls='-', label= 'LFS')
                axs[2].axvline(x= 0, color='black', lw=3, ls='-')
                axs[0].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[1].axvline(x= 180, color='brown', lw=3, ls='-', label= 'HFS')
                axs[2].axvline(x= 180, color='brown', lw=3, ls='-')
            
            else:
                pass
            
            axs[1].legend(loc = 'best')


    
    def twpolfluxndS_plot(self, scan_style, pol_list, log_flag):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        
        if withshift == False and withseries == True:
            
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
                
                # pol_list_a = []
                # for i in range(48):
                #     pol_list_a.append('{}'.format(25 + i))
                
                
                # self.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
                ang_list = self.data['angle']['angle_list']
                
                for x in dircomp[key_a]:
                    keylist_a.append('{:.3f}'.format(x))
                
                for ta in keylist_a:
                    
                    keylist_b = []
                    
                    for x in dircomp[key_b]:
                        keylist_b.append('{:.3f}'.format(x))
                    
                    
                    iter_key, color_dic= self.twa.twpolplot_prep(ta = ta, 
                    keylist_b = keylist_b, scan_style = scan_style)
                    
                    
                    print('check:')
                    print(iter_key)
                    print(color_dic)

                    self.twpolfluxndS_method(iterlist = iter_key, cl_dic = color_dic, 
                                scan_style = scan_style, ang_list = ang_list, 
                                pol_list = pol_list, plot_option = 'core', log_flag = log_flag)
                    
                    self.twpolfluxndS_method(iterlist = iter_key, cl_dic = color_dic, 
                                scan_style = scan_style, ang_list = ang_list, 
                                pol_list = pol_list, plot_option = 'inner leg', log_flag = log_flag)
                    
                    self.twpolfluxndS_method(iterlist = iter_key, cl_dic = color_dic, 
                                scan_style = scan_style, ang_list = ang_list, 
                                pol_list = pol_list, plot_option = 'outer leg', log_flag = log_flag)
             
            else:
                print('neteTS_plot, please check the series flag')
        