# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 18:25:23 2025

@author: ychuang
"""

# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.offsetbox import AnchoredText

from SOLPS_input.header import *
from twscan_module.twinscan_prepare import twscan_assist




class twscan_opacity_polplot:
    
    
    def __init__(self, DF, data, twa: twscan_assist):
        
        self.DF = DF
        self.data = data
        self.twa = twa
    
    
    
    def twscan_opacity_polplot_method(self, iterlist, cl_dic, scan_style, ang_list, plot_option, format_option):
        
        if format_option == '3x1':
            
            fig, axs = plt.subplots(3, 1)
        
        elif format_option == '2x1':
            
            fig, axs = plt.subplots(2, 1)
        
        elif format_option == '1x1':
            
            fig, axs = plt.subplots()
        
        
        
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']

        
        plt.subplots_adjust(hspace=.0)
        anchored_text_1 = AnchoredText('{}'.format('Density pedestal width [mm]'), 
                                              loc='upper right')
        anchored_text_2 = AnchoredText('{}'.format('Neutral penetration length [mm]'), 
                                              loc='lower right')
        anchored_text_3 = AnchoredText('{}'.format('Neutral opaqueness'), 
                                              loc='upper right')
        anchored_text_4 = AnchoredText('{}'.format('Neutral density'), 
                                              loc='upper right')
        
        for aa in iterlist:
            


            """
            label= 'core density {} $10^{19}$'.format(aa)
            
            """
            
            # axs.legend(loc= 'lower left', fontsize=10)
            
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
            
            
            efold = self.data['opacity_poloidal'][aa]['efold_length']*pow(10, 3)
            width = self.data['opacity_poloidal'][aa]['pedestal_width']*pow(10, 3)*2
            opq = self.data['opacity_poloidal'][aa]['dimensionless_opaqueness']
            nd_sep = self.data['opacity_poloidal'][aa]['neutral_density']
                
       

            if scan_style == 'denscan':
                
                title_ap = float(ap)*pow(10, 5)
                label_ad = float(ad)*pow(10, 20)
                
                if format_option == '3x1':
                    
                    axs[0].set_title('Particle flux scan with heat flux = {:.3E} W'.format(title_ap))
                    axs[0].plot(ang_list, width, color = cl_dic[ad], 
                             label= '{:.3E} (1/s)'.format(label_ad))
                    axs[1].plot(ang_list, efold, color = cl_dic[ad])
                    axs[2].plot(ang_list, opq, color = cl_dic[ad])
                    axs[0].set_ylim(0, 20)
                    axs[1].set_ylim(0, 20)
                    axs[2].set_ylim(0, 4)
                    axs[2].set_xlabel('poloidal angle')
                    axs[0].legend(loc = 'best')
                
                elif format_option == '2x1':
                    
                    axs[0].set_title('Particle flux scan with heat flux = {:.3E} W'.format(title_ap))
                    axs[0].plot(ang_list, efold, color = cl_dic[ad])
                    axs[1].plot(ang_list, opq, color = cl_dic[ad], 
                             label= '{:.3E} (1/s)'.format(label_ad))
                    axs[0].set_ylim(0, 20)
                    axs[1].set_ylim(0, 4)
                    axs[1].set_xlabel('poloidal angle')
                    axs[1].legend(loc = 'upper center')
                
                elif format_option == '1x1':
                    
                    axs.set_title('Particle flux scan with heat flux = {:.3E} W'.format(title_ap))
                    axs.plot(ang_list, nd_sep, color = cl_dic[ad], label= '{:.3E} (1/s)'.format(label_ad))
                    axs.set_xlabel('poloidal angle')
                    axs.legend(loc = 'upper center')
                    
                
            elif scan_style == 'tempscan':
                
                title_ap = float(ap)*pow(10, 20)
                label_ad = float(ad)*pow(10, 5)
                
                if format_option == '3x1':
                    
                    axs[0].set_title('Heat flux scan with particle flux = {:.3E} (1/s)'.format(title_ap))
                    axs[0].plot(ang_list, width, color = cl_dic[ad], 
                             label= '{:.3E} W'.format(label_ad))
                    axs[1].plot(ang_list, efold, color = cl_dic[ad])
                    axs[2].plot(ang_list, opq, color = cl_dic[ad])
                    axs[0].set_ylim(0, 20)
                    axs[1].set_ylim(0, 20)
                    axs[2].set_ylim(0, 4)
                    axs[2].set_xlabel('poloidal angle')
                    axs[0].legend(loc = 'best')
                    
                
                elif format_option == '2x1':
                    
                    axs[0].set_title('Heat flux scan with particle flux = {:.3E} (1/s)'.format(title_ap))
                    axs[0].plot(ang_list, efold, color = cl_dic[ad])
                    axs[1].plot(ang_list, opq, color = cl_dic[ad], 
                             label= '{:.3E} W'.format(label_ad))
                    axs[0].set_ylim(0, 20)
                    axs[1].set_ylim(0, 4)
                    axs[1].set_xlabel('poloidal angle')
                    axs[1].legend(loc = 'upper center')
                
                
                elif format_option == '1x1':
                    
                    axs.set_title('Heat flux scan with particle flux = {:.3E} (1/s)'.format(title_ap))
                    axs.plot(ang_list, nd_sep, color = cl_dic[ad], label= '{:.3E} W'.format(label_ad))
                    axs.set_xlabel('poloidal angle')
                    axs.legend(loc = 'upper center')
            
            else:
                print('neudenplot_method, please check the scan parameter')
        
        
        
        
        if scan_style == 'denscan':
            
                       
            if format_option == '3x1':
                
                
                axs[0].add_artist(anchored_text_1)
                axs[1].add_artist(anchored_text_2)
                axs[2].add_artist(anchored_text_3)
                axs[0].axvline(x= 0, color='black', lw=3, ls='-')
                axs[1].axvline(x= 0, color='black', lw=3, ls='-', label= 'LFS')
                axs[2].axvline(x= 0, color='black', lw=3, ls='-')
                axs[0].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[1].axvline(x= 180, color='brown', lw=3, ls='-', label= 'HFS')
                axs[2].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[1].legend(loc = 'best')
            
            elif format_option == '2x1':
                

                axs[0].add_artist(anchored_text_2)
                axs[1].add_artist(anchored_text_3)
                axs[0].axvline(x= 0, color='black', lw=3, ls='-')
                axs[1].axvline(x= 0, color='black', lw=3, ls='-', label= 'LFS')
                axs[0].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[1].axvline(x= 180, color='brown', lw=3, ls='-', label= 'HFS')
                axs[1].legend(loc = 'upper center')
            
            elif format_option == '1x1':
                

                axs.add_artist(anchored_text_4)
                axs.axvline(x= 0, color='black', lw=3, ls='-', label= 'LFS')
                axs.axvline(x= 180, color='brown', lw=3, ls='-', label= 'LFS')
                axs.legend(loc = 'upper center')
                
            
            
        elif scan_style == 'tempscan':
            
            
            if format_option == '3x1':


                axs[0].add_artist(anchored_text_1)
                axs[1].add_artist(anchored_text_2)
                axs[2].add_artist(anchored_text_3)
                axs[0].axvline(x= 0, color='black', lw=3, ls='-')
                axs[1].axvline(x= 0, color='black', lw=3, ls='-', label= 'LFS')
                axs[2].axvline(x= 0, color='black', lw=3, ls='-')
                axs[0].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[1].axvline(x= 180, color='brown', lw=3, ls='-', label= 'HFS')
                axs[2].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[0].set_ylim(0, 20)
                axs[1].set_ylim(0, 20)
                axs[2].set_ylim(0, 4)
                axs[1].legend(loc = 'best')
                
            
            elif format_option == '2x1':
                
                
                axs[0].add_artist(anchored_text_2)
                axs[1].add_artist(anchored_text_3)
                axs[0].axvline(x= 0, color='black', lw=3, ls='-')
                axs[1].axvline(x= 0, color='black', lw=3, ls='-', label= 'LFS')
                axs[0].axvline(x= 180, color='brown', lw=3, ls='-')
                axs[1].axvline(x= 180, color='brown', lw=3, ls='-', label= 'HFS')
                axs[0].set_ylim(0, 20)
                axs[1].set_ylim(0, 4)
                axs[1].legend(loc = 'upper center')
            
            
            elif format_option == '1x1':
                

                axs.add_artist(anchored_text_4)
                axs.axvline(x= 0, color='black', lw=3, ls='-', label= 'LFS')
                axs.axvline(x= 180, color='brown', lw=3, ls='-', label= 'HFS')
                axs.set_xlabel('poloidal angle')
                axs.legend(loc = 'upper center')
    
    
    
    
    
    
    
    def twscan_opacity_polplot(self, scan_style, plot_option, format_option):
        
        
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
                    
                    
                    iter_key, color_dic= self.twa.twscan_plot_prep(ta = ta, 
                    keylist_b = keylist_b, scan_style = scan_style)
                    
                    
                    print('check:')
                    print(iter_key)
                    print(color_dic)
                    # print(label_dic)
                    
                    self.twscan_opacity_polplot_method(iterlist = iter_key, cl_dic = color_dic, 
                                scan_style = scan_style, ang_list = ang_list, 
                    plot_option = plot_option, format_option = format_option)
                    
             
            else:
                print('neteTS_plot, please check the series flag')



