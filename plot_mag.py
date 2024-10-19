# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 13:36:06 2024

@author: ychuang
"""


from R_diff_calc import Diff_R_calc
import load_B2_data_method as lBdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText


class plot_mag(Diff_R_calc):
    def __init__(self, DefaultSettings, loadDS):
        Diff_R_calc.__init__(self, DefaultSettings, loadDS)






    def plot_mag(self, pol_list):
                   
        if self.withshift == True and self.withseries == False:
            
            # pol_list_a = []
            # for i in range(36):
            #     pol_list_a.append('{}'.format(28 + i))
                  
            # xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
            
            
            flux_list = ['B', 'bx', 'bz']
            
            
            fig, axs = plt.subplots(3, 2)
            
            color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                         'dot7': 'blue'}
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8'}
            
            B_text = AnchoredText('{}'.format('B: [T]'), 
                                         loc='upper center')
            
            bx_text = AnchoredText('{}'.format('$B_{pol}$: [T]'), 
                                         loc='upper center')
            
            bz_text = AnchoredText('{}'.format('$B_{tor}$: [T]'), 
                                         loc='lower center')
            
            text_list = [B_text, bx_text, bz_text]
            
            for ind, dat_name in enumerate(flux_list):
                
                for aa in xl.data['dircomp']['multi_shift']:
                    
                    ang_list = xl.data['angle']['angle_list'][aa]
                    
                    if ind == 2:
                        A_label_bol = False
                    else:
                        A_label_bol = True
                    
                    dat = xl.data['iout_data'][dat_name][aa]
                    
                    xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                        input_dat = dat, art_text = text_list[ind], axs = axs[ind, 0], 
                        color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = A_label_bol)
                        
                
                axs[ind, 0].legend(loc= 'lower left')
            
            axs[2, 0].set_xlabel('poloidal angle')
            axs[0, 0].set_title('magnetic strength at the separatrix')
            
            plt.subplots_adjust(hspace=.0)
            
            
            Bxratio_text = AnchoredText('{}'.format('$B_{pol}/B_{pol}^{MAST}$'), 
                                         loc='lower left')
            
            Bzratio_text = AnchoredText('{}'.format('$B_{tor}/B_{tor}^{MAST}$'), 
                                         loc='lower left')
            
            xzratio_text = AnchoredText('{}'.format('$B_{pol}/B_{tor}$'), 
                                         loc='lower center')
            
            ratiotext_list = [Bxratio_text, Bzratio_text, xzratio_text]
            
            
            mag_list = ['bx', 'bz']
            
            
            for ind, dat_name in enumerate(mag_list):
                
                for aa in xl.data['dircomp']['multi_shift']:
                    
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
                        
                    bx_mag = xl.data['b2wdat']['bx'][aa][psi_dic['st']:psi_dic['ed'], sk:sd]
                    bz_mag = xl.data['iout_data']['bz'][aa][psi_dic['st']:psi_dic['ed'], sk:sd]
                    
                    axs[2, 1].add_artist(ratiotext_list[2])
                    
                    mag_ratio = np.divide(bx_mag[0, :], bz_mag[0, :])
                    
                    axs[2, 1].plot(ang_list, mag_ratio, linestyle= '-', 
                        color= color_dic[aa])
                    
                    axs[2, 1].legend(loc= 'lower left') 
                    
            axs[2, 1].set_xlabel('poloidal angle')
            axs[0, 1].set_title('magnetic strength ratio at the separatrix') 
                
            plt.subplots_adjust(hspace=.0)