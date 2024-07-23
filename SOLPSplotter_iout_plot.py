# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 19:21:34 2024

@author: ychuang
"""

import matplotlib.pyplot as plt
from SOLPSplotter_iout_dataprocess import iout_process
from matplotlib.offsetbox import AnchoredText
import numpy as np



class Plotiout(iout_process):
    def __init__(self, DefaultSettings, loadDS):
        iout_process.__init__(self, DefaultSettings, loadDS)
        
        self.Publish = DefaultSettings['Publish']
        self.data['DefaultSettings']['Publish'] = self.Publish
    
       
    
    
    def iout_4subplot_multi(self, psi_dic, pol_list, data4_list, color_dic, 
                            A_dic, text_list, graph_dic, topic, series_dat):
        

        self.flux_iout_loader()
    
        self.calc_pol_angle(pol_list = pol_list, plot_angle= False)
        
        
        fig, axs = plt.subplots(4, 1)
        
        
        for ind, dat_name in enumerate(data4_list):
            
            for aa in series_dat:
                
                ang_list = self.data['angle']['angle_list'][aa]
                
                if ind == 3:
                    A_label_bol = False
                else:
                    A_label_bol = True
                
                in_dat = self.data['iout_data'][dat_name][aa]
                
                self.flux_poloidal_sub(itername = aa, pol_list = pol_list, ang_list = ang_list,
                    art_text = text_list[ind], axs = axs[ind], nnp = 1, input_dat = in_dat, input_ls = '-',
                    color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = A_label_bol)
            if topic == 'all_fluxes':
                
                if ind % 2 == 0:
                    
                    axs[ind].axhline(y=0, color = 'black', linestyle = '--', label= '$\Gamma_x$ = 0')
                else:
                    pass
            else:
                pass
            
            axs[ind].legend(loc= 'upper right')
        
        axs[3].set_xlabel('{}'.format(graph_dic['xlabel']))
        axs[0].set_title('{}'.format(graph_dic['title']))
        
        
        plt.subplots_adjust(hspace=.0)
    
    
    
    def iout_2subplot_multi(self, psi_dic, pol_list, data2_list, color_dic, 
                            A_dic, text_list, graph_dic, topic, series_dat):
            
            self.flux_iout_loader()
        
            self.calc_pol_angle(pol_list = pol_list, plot_angle= False)
            
            fig, axs = plt.subplots(2, 1)
            
            for ind, dat_name in enumerate(data2_list):
            
                for aa in series_dat:
                    
                    ang_list = self.data['angle']['angle_list'][aa]
                    
                    in_dat = self.data['iout_data'][dat_name][aa]
                    
                    if ind == 1:
                        A_label_bol = False
                    else:
                        A_label_bol = True
                    
                    self.flux_poloidal_sub(itername = aa, pol_list = pol_list, ang_list = ang_list,
                    nnp = 1, input_dat = in_dat, input_ls = '-', art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = A_label_bol)
                
                if topic == 'fluxes_geo' or topic == 'fluxes_no_geo':
                    
                    axs[ind].legend(loc= 'upper right')
                    
                elif topic == 'annual_coe':
                    
                    axs[ind].legend(loc= 'upper left')
                
                else:
                    pass
                
                
            axs[1].set_xlabel('{}'.format(graph_dic['xlabel']))
            axs[0].set_title('{}'.format(graph_dic['title']))
            
            
            plt.subplots_adjust(hspace=.0)
        
     
    
    def plot_neuden_multi(self, psi_dic, pol_list, color_dic, 
                            A_dic, text_list, graph_dic, solrange, series_dat):
        
        self.calc_pol_angle(pol_list = pol_list, plot_angle= False)
        fig, axs = plt.subplots()
        
        for aa in series_dat:
            
            neuden_data_a = []
            neuden_data_b = []
            
            for kt in pol_list:
                
                neuden_data = self.data['ft44'][aa]['dab2'][int(kt), psi_dic['st']:psi_dic['ed']]
                
                neuden_data_a.append(neuden_data.max())
                neuden_data_b.append(neuden_data.min())
            
            sk = int(pol_list[0])
            sd = int(pol_list[-1]) + 1
            
            ang_list = self.data['angle']['angle_list'][aa]
        
        
            neuden_dat = np.transpose(self.data['ft44'][aa]['dab2'][sk:sd, 
                                                    psi_dic['st']:psi_dic['ed'], 0])
            
            axs.add_artist(text_list[0])
            
            axs.plot(ang_list, neuden_dat[0, :], linestyle='-', color= color_dic[aa], 
                     label = 'A= {}'.format(A_dic[aa]))
            
            
            if solrange == True:
                               
                axs.fill_between(ang_list, neuden_data_a, neuden_data_b, 
                                  color= color_dic[aa], alpha = 0.4)
                
                axs.plot(ang_list, neuden_dat[-1, :], '-', color= color_dic[aa])
            
            else:
                pass
                
            axs.set_xlabel('{}'.format(graph_dic['xlabel']))
            axs.legend(loc= 'upper right')
            axs.set_title('{}'.format(graph_dic['title']))
    
            
    
    def iout_3x2subplot_multi(self, psi_dic, pol_list, data_list_A, data_list_B, color_dic, 
                            A_dic, text_list, graph_dic, topic, series_dat):
        
          
        self.flux_iout_loader()
        
        self.calc_pol_angle(pol_list = pol_list, plot_angle= False)
        
        fig, axs = plt.subplots(3, 2)
        
        
        
        for ind, dat_name in enumerate(data_list_A):
            
            for aa in series_dat:
                
                ang_list = self.data['angle']['angle_list'][aa]
                
                if ind == 2:
                    A_label_bol = False
                else:
                    A_label_bol = True
                
                dat = self.data['iout_data'][dat_name][aa]
                
                self.flux_poloidal_sub(itername = aa, pol_list = pol_list, ang_list = ang_list,
                    input_dat = dat, art_text = text_list[ind], axs = axs[ind, 0], 
                    color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = A_label_bol)
                    
            
            axs[ind, 0].legend(loc= 'lower left')
        
        axs[2, 0].set_xlabel('{}'.format(graph_dic['xlabel']))
        axs[0, 0].set_title('{}'.format(graph_dic['title']))
        
        plt.subplots_adjust(hspace=.0)
        
               
        
        for ind, dat_name in enumerate(mag_list):
            
            for aa in series_dat:
                
                if aa == 'org':
                    
                    pass
                
                else:
                    ang_list = xl.data['angle']['angle_list'][aa]
                    
                    
                    sk = int(pol_list_a[0])
                    sd = int(pol_list_a[-1]) + 1
                        
                    org_case_mag = xl.data['iout_data'][dat_name]['org'][psi_dic['st']:psi_dic['ed'], sk:sd]
                    shift_case_mag = xl.data['iout_data'][dat_name][aa][psi_dic['st']:psi_dic['ed'], sk:sd]
                    
                    
                    mag_ratio = np.divide(shift_case_mag[0, :], org_case_mag[0, :])
                    
                    axs[ind, 1].add_artist(ratiotext_list[ind])
                    
                    axs[ind, 1].plot(ang_list, mag_ratio, linestyle='-', 
                        color= color_dic[aa])
            
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                
                sk = int(pol_list_a[0])
                sd = int(pol_list_a[-1]) + 1
                    
                bx_mag = xl.data['iout_data']['bx'][aa][psi_dic['st']:psi_dic['ed'], sk:sd]
                bz_mag = xl.data['iout_data']['bz'][aa][psi_dic['st']:psi_dic['ed'], sk:sd]
                
                axs[2, 1].add_artist(ratiotext_list[2])
                
                mag_ratio = np.divide(bx_mag[0, :], bz_mag[0, :])
                
                axs[2, 1].plot(ang_list, mag_ratio, linestyle= '-', 
                    color= color_dic[aa])
                
                axs[2, 1].legend(loc= 'lower left') 
                
        axs[2, 1].set_xlabel('poloidal angle')
        axs[0, 1].set_title('magnetic strength ratio at the separatrix') 
            
        plt.subplots_adjust(hspace=.0)
            
        
       

"""

if topic == 'all_fluxes':
    
    
    if xl.withshift == True and xl.withseries == False:
        
        series_dat = xl.data['dircomp']['multi_shift']
        
        quant_list, text_list, graph_dic = graph_item(topic = topic)
      
        xl.iout_4subplot_multi(psi_dic = psi_dic, pol_list = pol_list_a, data4_list = quant_list, color_dic = color_dic, 
                                A_dic = A_dic, text_list = text_list, graph_dic = graph_dic, 
                                topic = topic , series_dat = series_dat)

elif topic == 'fluxes_geo':
    
    
    if xl.withshift == True and xl.withseries == False:
        
        series_dat = xl.data['dircomp']['multi_shift']
        
        quant_list, text_list, graph_dic = graph_item(topic = topic)
      
        xl.iout_2subplot_multi(psi_dic = psi_dic, pol_list = pol_list_a, data2_list = quant_list, color_dic = color_dic, 
                                A_dic = A_dic, text_list = text_list, graph_dic = graph_dic, 
                                topic = topic , series_dat = series_dat)
    
elif topic == 'fluxes_no_geo':
    
    
    if xl.withshift == True and xl.withseries == False:
        
        series_dat = xl.data['dircomp']['multi_shift']
        
        quant_list, text_list, graph_dic = graph_item(topic = topic)
      
        xl.iout_2subplot_multi(psi_dic = psi_dic, pol_list = pol_list_a, data2_list = quant_list, color_dic = color_dic, 
                                A_dic = A_dic, text_list = text_list, graph_dic = graph_dic, 
                                topic = topic , series_dat = series_dat)
    
elif topic == 'neu_den':
    
    if xl.withshift == True and xl.withseries == False:
        
        series_dat = xl.data['dircomp']['multi_shift']
        
        text_list, graph_dic = graph_item(topic = topic)
    
    
        xl.plot_neuden_multi(psi_dic = psi_dic, pol_list = pol_list_a, color_dic = color_dic, 
                                A_dic = A_dic, text_list = text_list, graph_dic =  graph_dic, 
                                solrange = False, series_dat = series_dat)

elif topic == 'coe_check':
    
    
    if xl.withshift == True and xl.withseries == False:
        
        series_dat = xl.data['dircomp']['multi_shift']
        
        quant_list, text_list, graph_dic = graph_item(topic = topic)
      
        xl.iout_4subplot_multi(psi_dic = psi_dic, pol_list = pol_list_a, data4_list = quant_list, color_dic = color_dic, 
                                A_dic = A_dic, text_list = text_list, graph_dic = graph_dic, 
                                topic = topic, series_dat = series_dat)

elif topic == 'annual_coe':
    
    
    if xl.withshift == True and xl.withseries == False:
        
        series_dat = xl.data['dircomp']['multi_shift']
        
        quant_list, text_list, graph_dic = graph_item(topic = topic)
      
        xl.iout_2subplot_multi(psi_dic = psi_dic, pol_list = pol_list_a, data2_list = quant_list, color_dic = color_dic, 
                                A_dic = A_dic, text_list = text_list, graph_dic = graph_dic, 
                                topic = topic , series_dat = series_dat)




----------- data2 old code ---------------

if topic == 'fluxes_geo':
    
    axs[0].axhline(y=0, color = 'black', linestyle = '--', label= '$\Gamma_x$ = 0')
    
elif topic == 'fluxes_no_geo':
    
    axs[0].axhline(y=0, color = 'black', linestyle = '--', 
                   label= '$\Gamma_{\Theta}$ = 0')

else:
    pass


for aa in series_dat:
    
    ang_list = self.data['angle']['angle_list'][aa]
    
    in_dat = self.data['iout_data'][data2_list[1]][aa]
    
    self.flux_poloidal_sub(itername = aa, pol_list = pol_list, ang_list = ang_list,
      nnp = 1, input_dat = in_dat, input_ls = '-', art_text = text_list[1], axs = axs[1], 
    color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = False)



"""


    


        
            
            
    
    