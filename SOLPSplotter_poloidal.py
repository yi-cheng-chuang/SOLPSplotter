# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 17:25:25 2024

@author: ychuang
"""

from SOLPSplotter_fit import profile_fit
import matplotlib.pyplot as plt
import SOLPS_set as ss
import numpy as np


class poloidal_plot(profile_fit):
    def __init__(self, DefaultSettings, loadDS):
        profile_fit.__init__(self, DefaultSettings, loadDS)
        
        self.Publish = DefaultSettings['Publish']
        self.data['DefaultSettings']['Publish'] = self.Publish
    
    
    def set_plot(self):
        if self.Publish == 'b2plottersetting':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 5, markersize= 9)
            plt.rcParams.update({'font.size': 16})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
  
        else:
            print('Publish setting is incorrect or add another setting')



    
    # def opacity_study_unit(self):
    #     unit = {'efold_length_psiN': 'Neutral penetration length ($\psi_N$)',
    #             'pedestal_width_psiN': 'Pedestal width ($\psi_N$)',
    #               'dimensionless_opaqueness': 'Experimental opaqueness', 
    #               'neutral_density': 'Neutral density ${n_D}$ (m$^{-3}$)', 
    #               'electron_pedestal_density': 'Electron pedestal density: $n_{ped}$ (m$^{-3}$)',
    #               'temperature_pedestal_width': 'Temperature pedestal width: $\Delta T$',
    #               'flux_expansion': 'Flux expansion',
    #               'efold_length': '$\lambda_{n_D}$ [mm]',
    #               'pedestal_width': '$\Delta n_e$ [mm]',
                  
    #               }
    #     return unit
    
    
   
    def opacity_poloidal_plot_method(self, item, pol_angle, result_dic, color_code, 
                                 A_value, unit_dic):
        
        if item == 'electron_pedestal_density':
            plt.scatter(result_dic['electron_pedestal_density'], result_dic['efold_length'], 
                     color= color_code, label= 'aspect ratio = {} {}'.format(A_value, 'efold length'))
            plt.scatter(result_dic['electron_pedestal_density'], result_dic['pedestal_width'], 
                     color= color_code, label= 'aspect ratio = {} {}'.format(A_value, 'pedestal width'))
            plt.xlabel('electron_pedestal_density: $n_{ped}$ (m$^{-3}$)')
            plt.title('{} and {} verses {} from {:.2f} to {:.2f}'.format('efold length', 
                                                      'pedestal width', 'electron pedestal density', max(pol_angle), min(pol_angle)))
            # plt.title('{} and {} verses {}'.format('efold length', 
            #                 'pedestal width', 'electron pedestal density'))
            plt.legend()
            
        elif item == 'pedestal_width_psiN':
            # plt.scatter(pol_loc[aa], np.round_(result_dic[i][aa], 2), 
            #              label= 'modified {} m case'.format(change_ver_dic[aa]))
            # plt.errorbar(pol_loc[aa], np.round_(result_dic[i][aa], 2), yerr= np.std(result_dic[i][aa]), fmt= 'o', 
            #               label= '{}'.format(i))
            plt.plot(pol_angle, np.round_(result_dic[item], 2), '-', color = color_code,
                  label= 'aspect ratio = {}'.format(A_value))
            plt.title('{} verses poloidal angle from {:.2f} to {:.2f}'.format(unit_dic[item], 
                                              max(pol_angle), min(pol_angle)))
            plt.xlabel('poloidal angle')
            plt.legend()

        elif item == 'pedestal_width' or item == 'efold_length':
            new_r = result_dic[item]*1000
            plt.plot(pol_angle, new_r, '-', color= color_code, 
                     label= 'aspect ratio = {}'.format(A_value))
            plt.title('{} verses poloidal angle from {:.2f} to {:.2f}'.format(unit_dic[item], 
                                              max(pol_angle), min(pol_angle)))                   
            plt.xlabel('poloidal angle')
            plt.legend()
            
        else:
            plt.plot(pol_angle, result_dic[item], '-', color= color_code ,
                  label= 'aspect ratio = {}'.format(A_value))
            # plt.plot(pol_loc[aa], result_dic[i][aa], '-', color= color_dic[aa])
            plt.title('{} verses poloidal angle from {:.2f} to {:.2f}'.format(unit_dic[item], 
                                              max(pol_angle), min(pol_angle)))
            # plt.title('{} verses poloidal angle'.format(unit_dic[i]))
            # plt.title('{}'.format(unit_dic[i]))                   
            plt.xlabel('poloidal angle')
            plt.legend()
    
    
    
    def poloidal_label(self, angle_fix, item, xpoint_fix):
        
        if max(angle_fix) > 90 and item != 'electron_pedestal_density' and item != 'width_relation':
            plt.axvline(x= 90, color='black',lw=3, ls='--', label= '90')
        else:
            pass
        if max(angle_fix) > 180 and item != 'electron_pedestal_density':
            plt.axvline(x= 180, color='seagreen',lw=3, ls='--', label= 'inner midplane')
        else:
            pass
        
        if max(angle_fix) > 240 and item != 'electron_pedestal_density':
            plt.axvline(x= 240, color='seagreen',lw=3, ls='--', label= 'poloidal angle 240')
        else:
            pass
        
        
        if min(angle_fix) < 0 and item != 'electron_pedestal_density':
            plt.axvline(x= 0, color='darkorange',lw=3, ls='--', label= 'outer midplane')
        else:
            pass
        
        if min(angle_fix) < -70 and item != 'electron_pedestal_density':
            plt.axvline(x= xpoint_fix, color='darkblue',lw=3, ls='--', label= 'xpoint angle')
            plt.axvline(x= xpoint_fix + 360, color='darkblue',lw=3, ls='--')
        else:
            pass
        # plt.ylabel('{}'.format(unit_dic[i]))
        
        
        
        plt.legend()
        
        
    def neuden_data_check(self, pol_list):
        
        
        if self.withshift == False and self.withseries == False:
            
            list_len = len(pol_list)
            cp_dat = np.zeros([list_len, 3])
            
            
            neuden_dat = self.data['opacity_poloidal']['neutral_density']
                        
            cp_dat[:, 0] = int(pol_list)
            cp_dat[:, 1] = self.data['angle']['angle_list']
            cp_dat[:, 2] = neuden_dat
            
            self.data['neuden_angle'] = cp_dat
        
        elif self.withshift == True and self.withseries == False:
            cp_dat_dic = {}
            list_len = len(pol_list)
            for aa in self.data['dircomp']['multi_shift']:
                
                list_len = len(pol_list)
                cp_dat = np.zeros([list_len, 3])
                
                
                neuden_dat = self.data['opacity_poloidal'][aa]['neutral_density']
                            
                cp_dat[:, 0] = [int(i) for i in pol_list]
                cp_dat[:, 1] = self.data['angle']['angle_list'][aa]
                cp_dat[:, 2] = neuden_dat
                
                cp_dat_dic[aa] = cp_dat
            
            self.data['neuden_angle'] = cp_dat_dic
            
    
    
    
    def opacity_poloidal_plot(self, log_flag, save_pdf):
        
        itemname = self.data['poloidal_itemname']
        # adj_list = list(result_dic.keys())
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        
        # print(adj_list)
        for i in itemname:
            plt.figure(figsize=(7,7))
            if log_flag:
                plt.yscale('log')
            else:
                pass
            
            if self.withshift == False and self.withseries == False:
                
                result = self.data['opacity_poloidal']

                
                unit = self.opacity_study_unit()
                pol_loc = self.data['angle']['angle_list']
                xpoint = self.data['angle']['xpoint_angle']
                a_shift = self.data['dircomp']['a_shift']
                A_val = A_dic[a_shift]
                color = color_dic[a_shift]
                
                self.opacity_poloidal_plot_method(item = i, pol_angle = pol_loc, 
            result_dic = result, color_code = color, A_value = A_val, unit_dic = unit)
                
                self.poloidal_label(angle_fix= pol_loc, item= i, xpoint_fix = xpoint)
                
                
        
            
            elif self.withshift == True and self.withseries == False:
                
                
                for aa in self.data['dircomp']['multi_shift']:
                    
                    result = self.data['opacity_poloidal'][aa]
                    unit = self.opacity_study_unit()
                    pol_loc = self.data['angle']['angle_list'][aa]
                    xpoint = self.data['angle']['xpoint_angle'][aa]
                    A_val = A_dic[aa]
                    color = color_dic[aa]
                    ang_fix = self.data['angle']['angle_list']['org']
                    xp_fix = self.data['angle']['xpoint_angle']['org']
                    
                    self.opacity_poloidal_plot_method(item = i, pol_angle = pol_loc, 
                result_dic = result, color_code = color, A_value = A_val, unit_dic = unit)
                
                
                self.poloidal_label(angle_fix= ang_fix, item= i, xpoint_fix = xp_fix)
                
                if save_pdf:
                    
                    fig_dir  = ss.set_figdir()
                    plt.savefig('{}/{}.pdf'.format(fig_dir, i), format='pdf')
            
            
            elif self.withshift == True and self.withseries == False:
                
                density_dic = {}
                for k in self.data['dircomp']['Attempt'].keys():
                    kk = float(k)*pow(10, 19)
                    density_dic[k] = kk
                    
                
                for aa in list(self.data['dircomp']['Attempt'].keys()):
                    
                    result = self.data['opacity_poloidal'][aa]
                    unit = self.opacity_study_unit()
                    pol_loc = self.data['angle']['angle_list']
                    xpoint = self.data['angle']['xpoint_angle']
                    a_shift = self.data['dircomp']['a_shift']
                    A_val = A_dic[a_shift]
                    color = color_dic[a_shift]
                    
                    self.opacity_poloidal_plot_method(item = i, pol_angle = pol_loc, 
                result_dic = result, color_code = color, A_value = A_val, 
                    unit_dic = unit, xpoint_loc = xpoint, angle_fix= pol_loc)
                
                
                self.poloidal_label(angle_fix= pol_loc, item= i, xpoint_fix = xpoint)
            
            
            elif self.withshift == True and self.withseries == True:
                print('poloidal_plot function is not there yet!')
            
            
            else:
                print('poloidal_plot function has a bug.')
        
        
        def opacity_poloidal_subplot_method(self, itemname):
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            
            if self.withshift == False and self.withseries == False:
                
                result = self.data['opacity_poloidal']

                
                unit = self.opacity_study_unit()
                pol_loc = self.data['angle']['angle_list']
                xpoint = self.data['angle']['xpoint_angle']
                a_shift = self.data['dircomp']['a_shift']
                A_val = A_dic[a_shift]
                color = color_dic[a_shift]
                
                self.opacity_poloidal_plot_method(item = itemname, pol_angle = pol_loc, 
            result_dic = result, color_code = color, A_value = A_val, unit_dic = unit)
                
                self.poloidal_label(angle_fix= pol_loc, item= i, xpoint_fix = xpoint)
                
                
        
            
            elif self.withshift == True and self.withseries == False:
                
                
                for aa in self.data['dircomp']['multi_shift']:
                    
                    result = self.data['opacity_poloidal'][aa]
                    unit = self.opacity_study_unit()
                    pol_loc = self.data['angle']['angle_list'][aa]
                    xpoint = self.data['angle']['xpoint_angle'][aa]
                    A_val = A_dic[aa]
                    color = color_dic[aa]
                    ang_fix = self.data['angle']['angle_list']['org']
                    xp_fix = self.data['angle']['xpoint_angle']['org']
                    
                    self.opacity_poloidal_plot_method(item = itemname, pol_angle = pol_loc, 
                result_dic = result, color_code = color, A_value = A_val, unit_dic = unit)
                
                
                self.poloidal_label(angle_fix= ang_fix, item= i, xpoint_fix = xp_fix)
        
        
        
                
        

        def data_reorder(iter_list, change_var, data_collect, char):
            withshift = char['withshift']
            withseries = char['withseries']
            series_flag = char['series_flag']
            
            if withshift == False and withseries == True:
                plt.figure(figsize=(7,7))
                color_list = ['red', 'salmon', 'orange', 'lime', 'green', 'darkgreen', 
                              'cyan', 'deepskyblue', 'navy', 'purple']
                if series_flag == 'change_den':
                    for p in range(len(change_var)):
                        x_cor = np.ones(len(iter_list))*change_var[p]
                        plt.scatter(x_cor, data_collect[:, p], color= color_list[p], 
                                        label= 'efold length for different core electron density {}'.format(change_var[p]))
                    plt.xlabel('Electron density: ${n_e}$ (m$^{-3}$)')
                    plt.title('dimensionless opaqueness verses shift core electron density')
                    plt.legend()
                    
                elif series_flag == 'eireneN':
                    for p in range(len(change_var)):
                        x_cor = np.ones(len(iter_list))*change_var[p]
                        plt.scatter(x_cor, data_collect[:, p], color= color_list[p], 
                                        label= '{}'.format(change_var[p]))
                    plt.xlabel('Eirene SOL particle number')
                    # plt.yscale('log')
                    plt.xscale('log')
                    # plt.title('dimensionless opaqueness verses different SOL eirene particle number')
                    plt.title('dimensionless opaqueness')
                    # plt.legend()
                
            elif withshift == True and withseries == False:
                plt.figure(figsize=(7,7))
                A_list = [1.4, 2.0, 2.4, 2.8, 3.4]
                color_list = ['red', 'orange', 'green', 'blue', 'purple']
                for p in range(len(change_var)):
                    x_cor = np.ones(len(iter_list))*A_list[p]
                    plt.scatter(x_cor, data_collect[:, p], color= color_list[p], 
                                label= 'aspect ratio= {}'.format(str(A_list[p])))
                plt.axvline(x= 0.7/0.5, color='salmon',lw=3, ls='--', 
                            label= 'MAST aspect ratio')
                plt.axvline(x= 1.67/0.67, color='darkgreen',lw=3, ls='--', label= 'D3D aspect ratio')
                plt.axvline(x= 0.68/0.22, color='cyan',lw=3, ls='--', label= 'C-Mod/ ITER aspect ratio')
                plt.axvline(x= 3.4, color='black',lw=3, ls='--', label= 'JT-60 aspect ratio')
                plt.xlabel('aspect ratio')
                # plt.ylabel('dimensionless opaqueness')
                # plt.title('dimensionless opaqueness verses different modify distance')
                plt.title('dimensionless opaqueness')
                plt.legend()
                
                
                
                
                
                
                
                
            
            
        
        
"""    
    


"""