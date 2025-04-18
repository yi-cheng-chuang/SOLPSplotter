# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 14:09:51 2025

@author: ychuang
"""




import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from numpy import ma
from matplotlib.offsetbox import AnchoredText
from matplotlib import colors, cm
from twscan_module.twinscan_prepare import twscan_assist
from contour_plot.contourplot_toolbox import contour_plot_method_collect
from matplotlib.colors import LogNorm




class eirene_binfocus_contour:


    
    def __init__(self, DF, data, twa: twscan_assist, cpmc: contour_plot_method_collect):
        
        self.DF = DF
        self.data = data
        self.twa = twa
        self.cpmc = cpmc


    
    def twscan_binfocus_eirene_contourplot_method(self, itername, scan_style, plot_option, axs, cmap, norm):
        
         
        if plot_option == 'Neuden contour':
            simu_dir = self.data['dirdata']['simudir'][itername[0]][itername[1]]
            dat = self.data['ft46'][itername[0]][itername[1]]['pdena'][:, 0]
            
            
        
        if self.DF.series_flag == 'twin_scan':
            
            if scan_style == 'tempscan':
                
                ad = itername[1]
                ap = itername[0]
                
            
            elif scan_style == 'denscan':
                
                ad = itername[0]
                ap = itername[1]
            
            else:
                print('neteTSplot_method, please check scan_style')
        
        else:
            ad = itername

        if scan_style == 'denscan':
            
            title_ap = float(ap)*pow(10, 5)
            
            self.cpmc.subplot_eirene_contour_method(simu_direc = simu_dir, data = dat, 
                                               itername = itername, axs = axs, cmap= cmap, norm = norm)
               
        elif scan_style == 'tempscan':
            
            title_ap = float(ap)*pow(10, 20)
                              
            self.cpmc.subplot_eirene_contour_method(simu_direc = simu_dir, data = dat, 
                                               itername = itername, axs = axs, cmap= cmap, norm = norm)
      
        else:
            print('neudenplot_method, please check the scan parameter')
            
        
            # plt.tight_layout()    
        
               



    def twscan_binfocus_eirene_contourplot(self, scan_style, plot_option):
        
        
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
                
                for x in dircomp[key_a]:
                    keylist_a.append('{:.3f}'.format(x))
                
                
                nx = self.data['b2fgeo']['nx']
                ny = self.data['b2fgeo']['ny']
                
                
                anchored_text_1 = AnchoredText('{}'.format('density pedestal width [mm]'), 
                                                      loc='lower center')
                anchored_text_2 = AnchoredText('{}'.format('neutral penetration length [mm]'), 
                                                      loc='lower right')
                anchored_text_3 = AnchoredText('{}'.format('dimensionless opaqueness'), 
                                                      loc='lower right')
                anchored_text_4 = AnchoredText('{}'.format('neutral density'), 
                                                      loc='lower right')
                
                
                
                
                for ta in keylist_a:
                    
                    keylist_b = []
                    
                    for x in dircomp[key_b]:
                        keylist_b.append('{:.3f}'.format(x))
                    
                    
                    iter_key, color_dic= self.twa.twscan_plot_prep(ta = ta, 
                    keylist_b = keylist_b, scan_style = scan_style)
                    
                    fig, axs = plt.subplots(1, 5, sharey=True)
                    
                    plt.subplots_adjust(hspace=.0)
                    
                    print('check:')
                    print(iter_key)
                    print(color_dic)
                    # print(label_dic)
                    
                    it = 0
                    
                    CPB = cm.viridis
                    Lnorm = LogNorm(vmax = pow(10, 19), vmin = pow(10, 14))
                        
                    for aa in iter_key:
                        
                        
                        self.twscan_binfocus_eirene_contourplot_method(itername = aa, scan_style = scan_style, 
                                    plot_option = plot_option, axs = axs[it], cmap= CPB, norm = Lnorm)
                        
                        it = it + 1
                    
                    fig.subplots_adjust(right=0.8)
                    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
                    smap = cm.ScalarMappable(Lnorm, CPB)    
                    fig.colorbar(smap, cax= cbar_ax)
                                        
                    
                    
                
                # plt.tight_layout() 
             
            else:
                print('neteTS_plot, please check the series flag')