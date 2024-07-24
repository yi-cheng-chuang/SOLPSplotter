# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 16:29:09 2024

@author: ychuang
"""



import matplotlib.pyplot as plt
from SOLPSplotter_iout_polplot import Plotiout
from matplotlib.offsetbox import AnchoredText
import numpy as np



class Present_ioutplot(Plotiout):
    def __init__(self, DefaultSettings, loadDS):
        Plotiout.__init__(self, DefaultSettings, loadDS)
        
        self.Publish = DefaultSettings['Publish']
        self.data['DefaultSettings']['Publish'] = self.Publish



    def set_plot(self):
        if self.Publish == 'b2plottersetting':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 3, markersize= 6)
            plt.rcParams.update({'font.size': 10})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
    
        else:
            print('Publish setting is incorrect or add another setting')
    
      
    
    
    def graph_item(self, topic):
        
        if topic == 'all_fluxes':
    
            quant_list = ['poloidal_flux','radial_flux', 'flux_x_0', 'flux_y_0']
            
            pol_text = AnchoredText('{}'.format('Poloidal flux $\Gamma_x$ [$m^{-1} s^{-1}$]'), 
                                         loc='upper center')
            
            rad_text = AnchoredText('{}'.format('Radial flux $\Gamma_y$ [$m^{-1} s^{-1}$]'), 
                                         loc='upper center')
            
            x0_text = AnchoredText('{}'.format('Poloidal flux $\sqrt{g}/h_x*\Gamma_x$ [$m*s^{-1}$]'), 
                                         loc='lower center')
            
            y0_text = AnchoredText('{}'.format('radial flux $\sqrt{g}/h_y*\Gamma_y$ [$m*s^{-1}$]'), 
                                         loc='upper center')
            
            text_list = [pol_text, rad_text, x0_text, y0_text]
            
            graph_dic = {'xlabel': 'poloidal angle', 'title': 'Particle flux at separatrix'}
            
            return quant_list, text_list, graph_dic
        
        elif topic == 'fluxes_geo':
    
            quant_list = ['flux_x_0', 'flux_y_0']
            
            pol_text = AnchoredText('{}'.format('Poloidal flux $\Gamma_x$ [$m^{-1} s^{-1}$]'), 
                                         loc='upper center')
            
            neu_text = AnchoredText('{}'.format('Neutral density [$m^{-3}$]'), 
                                         loc='upper center')
            
            rad_text = AnchoredText('{}'.format('Radial flux $\Gamma_y$ [$m^{-1} s^{-1}$]'), 
                                         loc='upper center')
            
            text_list = [pol_text, rad_text]
            
            graph_dic = {'xlabel': 'poloidal angle', 'title': 'Particle flux at separatrix'}
            
            return quant_list, text_list, graph_dic
        
        
        elif topic == 'fluxes_no_geo':
    
            quant_list = ['poloidal_flux', 'radial_flux']
            
            pol_text = AnchoredText('{}'.format('(a) Poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                         loc='upper center')
            
            neu_text = AnchoredText('{}'.format('Neutral density [$m^{-3}$]'), 
                                         loc='upper center')
            
            rad_text = AnchoredText('{}'.format('(b) Radial flux $\Gamma_r$ [$m^{-2} s^{-1}$]'), 
                                         loc='upper center')
            
            text_list = [pol_text, rad_text]
            
            graph_dic = {'xlabel': 'poloidal angle', 'title': 'Particle flux at separatrix'}
            
            return quant_list, text_list, graph_dic
        
        elif topic == 'neu_den':
            
            neu_text = AnchoredText('{}'.format('Neutral density [$m^{-3}$]'), 
                                         loc='upper center')
            
            text_list = [neu_text]
            
            graph_dic = {'xlabel': 'poloidal angle', 'title': 'neutral density'}
            
            return text_list, graph_dic
        
        elif topic == 'coe_check':
            
            quant_list = ['hx','hx_divide_sqrt_g', 'hy','hy_divide_sqrt_g']
    
            hx_text = AnchoredText('{}'.format('$h_x$: [m]'), 
                                         loc='upper center')
            
            hxg_text = AnchoredText('{}'.format('$h_x/ \sqrt{g}$: [$m^{-2}$]'), 
                                         loc='upper center')
            
            hy_text = AnchoredText('{}'.format('$h_y$: [m]'), 
                                         loc='upper center')
            
            hyg_text = AnchoredText('{}'.format('$h_y/ \sqrt{g}$: [$m^{-2}$]'), 
                                         loc='upper center')
            
            text_list = [hx_text, hxg_text, hy_text, hyg_text]
            
            graph_dic = {'xlabel': 'poloidal angle', 'title': 'coefficients at the separatrix'}
            
            return quant_list, text_list, graph_dic
        
        
        elif topic == 'annual_coe':
                
            quant_list = ['hx', 'hy']
            
            hx_text = AnchoredText('{}'.format('$h_x$'), 
                                         loc='upper right')
            
            hxg_text = AnchoredText('{}'.format('$h_x/ \sqrt{g}$: [$m^{-2}$]'), 
                                         loc='upper center')
            
            hy_text = AnchoredText('{}'.format('$h_y$'), 
                                         loc='upper center')
            
            hyg_text = AnchoredText('{}'.format('$h_y/ \sqrt{g}$: [$m^{-2}$]'), 
                                         loc='upper center')
            
            text_list = [hx_text, hy_text]
    
            graph_dic = {'xlabel': 'poloidal angle', 'title': 'geometry coefficients at the separatrix'}
            
            return quant_list, text_list, graph_dic
        
        
        elif topic == 'mag_pol':
            
            quant_list_A = ['B', 'bx', 'bz']
            
            B_text = AnchoredText('{}'.format('B: [T]'), 
                                         loc='upper center')
            
            bx_text = AnchoredText('{}'.format('$B_{pol}$: [T]'), 
                                         loc='upper center')
            
            bz_text = AnchoredText('{}'.format('$B_{tor}$: [T]'), 
                                         loc='lower center')
            
            text_list_A = [B_text, bx_text, bz_text]
            
            
            quant_list_B = ['bx', 'bz']
            
            
            Bxratio_text = AnchoredText('{}'.format('$B_{pol}/B_{pol}^{MAST}$'), 
                                          loc='lower left')
            
            Bzratio_text = AnchoredText('{}'.format('$B_{tor}/B_{tor}^{MAST}$'), 
                                          loc='lower left')
            
            xzratio_text = AnchoredText('{}'.format('$B_{pol}/B_{tor}$'), 
                                          loc='lower center')
            
            text_list_B = [Bxratio_text, Bzratio_text, xzratio_text]
            
            graph_dic = {'xlabel': 'poloidal angle', 
                         'title_A': 'magnetic strength at the separatrix',
                         'title_B': 'magnetic strength ratio at the separatrix'}
            
            quant_list_dic = {'A': quant_list_A, 'B': quant_list_B}
            text_list_dic = {'A': text_list_A, 'B': text_list_B}
            
            return quant_list_dic, text_list_dic, graph_dic
        
        elif topic == 'source':
            
            quant_list = 
        
        
    
    
    
    
    
    
    def iout_plotcomb(self, topic_label_dic ,topic, psi_dic, 
                      pol_list, color_dic, A_dic):
        
        
        if topic_label_dic[topic] == '4':
            
            if self.withshift == True and self.withseries == False:
                
               series_dat = self.data['dircomp']['multi_shift']
              
               quant_list, text_list, graph_dic = self.graph_item(topic = topic)
            
               self.iout_4subplot_multi(psi_dic = psi_dic, pol_list = pol_list, data4_list = quant_list, color_dic = color_dic, 
                                       A_dic = A_dic, text_list = text_list, graph_dic = graph_dic, 
                                       topic = topic , series_dat = series_dat)
        
        
        elif topic_label_dic[topic] == '2':
            
            if self.withshift == True and self.withseries == False:
                
                series_dat = self.data['dircomp']['multi_shift']
                
                quant_list, text_list, graph_dic = self.graph_item(topic = topic)
              
                self.iout_2subplot_multi(psi_dic = psi_dic, pol_list = pol_list, data2_list = quant_list, color_dic = color_dic, 
                                        A_dic = A_dic, text_list = text_list, graph_dic = graph_dic, 
                                        topic = topic , series_dat = series_dat)
        
        elif topic_label_dic[topic] == 'neuden_sp':
            
            if self.withshift == True and self.withseries == False:
                
                series_dat = self.data['dircomp']['multi_shift']
                
                text_list, graph_dic = self.graph_item(topic = topic)
            
            
                self.plot_neuden_multi(psi_dic = psi_dic, pol_list = pol_list, color_dic = color_dic, 
                                        A_dic = A_dic, text_list = text_list, graph_dic =  graph_dic, 
                                        solrange = False, series_dat = series_dat)
            
    
                
                
                
    def iout_plot(self, topic_label_dic, plotall,topic, psi_dic, 
                      pol_list, color_dic, A_dic):
        
        if plotall == True:
            for kk in list(topic_label_dic.keys()):
                
                print(type(kk))
                # print(list(topic_label_dic.keys()))
                
                self.iout_plotcomb(topic_label_dic ,topic, psi_dic, 
                                  pol_list, color_dic, A_dic)