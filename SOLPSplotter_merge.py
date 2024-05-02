# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 14:25:17 2024

@author: user
"""

from SOLPSplotter_fit import profile_fit
from matplotlib.offsetbox import AnchoredText
import SOLPS_set as ss
import os
import opacity_plot_method as opm
import matplotlib.pyplot as plt
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate
from scipy.optimize import curve_fit
import numpy as np


class paper_poloidal_plot(profile_fit):
    def __init__(self, DefaultSettings, loadDS):
        profile_fit.__init__(self, DefaultSettings, loadDS)
    
    
    def set_plot(self, plot_style):
        if plot_style == 'pol_subplot':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 4, markersize= 7)
            plt.rcParams.update({'font.size': 16})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
            # plt.rcParams['figure.figsize'] = 40, 12
            
  
        else:
            print('Publish setting is incorrect or add another setting')
    
    
     
    def paper_poloidal_method(self, item, pol_angle, result_dic, color_code, 
                                 A_value, unit_dic, ax, plot_order):
        
        
        
        
        if item == 'efold_length_psiN':
            
            anchored_text = AnchoredText('{}{}'.format(plot_order, unit_dic), loc=3)
        else:
            
            anchored_text = AnchoredText('{}{}'.format(plot_order, unit_dic), 
                                         loc=2)
            
        
        
        if item == 'pedestal_width' or item == 'efold_length':
            new_r = result_dic[item]*1000
            ax.plot(pol_angle, new_r, '-', color= color_code, label= 'A= {}'.format(A_value))
            # ax.legend(loc='upper right', bbox_to_anchor=(0.5, 0.5))
                    
        elif item == 'pedestal_width_psiN':
            ax.plot(pol_angle, np.round_(result_dic[item], 2), '-', 
                    color = color_code, label= 'A= {}'.format(A_value))
            # ax.legend(loc='upper right', bbox_to_anchor=(0.5, 0.5))
        
        else:
            ax.plot(pol_angle, result_dic[item], '-', color= color_code,
                    label= 'A= {}'.format(A_value))
            # ax.legend(loc='upper right', bbox_to_anchor=(0.5, 0.5))
        
        if item == 'dimensionless_opaqueness':
            
            ax.set_ylim(ymin = 0, ymax = 2 )
            
        elif item == 'flux_expansion':
            
            ax.set_ylim(ymin = 0, ymax = 10)
        
        elif item == 'pedestal_width':
            
            ax.set_ylim(ymax = 25)
        
        elif item == 'efold_length':
            
            ax.set_ylim(ymax = 35)
            ax.legend(loc= 'upper center', fontsize=10)
        
        
            
        else:
            pass
        
        
        
        
        
        
        ax.add_artist(anchored_text)
                   
    
    
    def paper_singleplot_method(self, item, pol_angle, result_dic, color_code, 
                                 A_value, unit_dic, plot_order, axs):
        
        
        # anchored_text = AnchoredText('{}{}'.format(plot_order, unit_dic), loc=2)
        
        if item == 'pedestal_width' or item == 'efold_length':
            new_r = result_dic[item]*1000
            axs.plot(pol_angle, new_r, '-', color= color_code, label= 'A= {}'.format(A_value))
            axs.set_xlabel('poloidal angle')
            # ax.legend(loc='upper right', bbox_to_anchor=(0.5, 0.5))
                    
        elif item == 'pedestal_width_psiN':
            axs.plot(pol_angle, np.round_(result_dic[item], 2), '-', 
                     color = color_code, label= 'A= {}'.format(A_value))
            axs.set_xlabel('poloidal angle')
            # ax.legend(loc='upper right', bbox_to_anchor=(0.5, 0.5))
        
        else:
            axs.plot(pol_angle, result_dic[item], '-', 
                     color= color_code, label= 'A= {}'.format(A_value))
            axs.set_xlabel('poloidal angle')
            # ax.legend(loc='upper right', bbox_to_anchor=(0.5, 0.5))
        
        # plt.add_artist(anchored_text)    
    
        
    
    def paper_poloidal_label(self, angle_fix, item, xpoint_fix, ax):
        
        # if max(angle_fix) > 90 and item != 'electron_pedestal_density' and item != 'width_relation':
        #     ax.axvline(x= 90, color='black',lw=3, ls='--')
        # else:
        #     pass
        if max(angle_fix) > 180 and item != 'electron_pedestal_density':
            ax.axvline(x= 180, color='gray',lw=3, ls='--', label= 'inner midplane')
        else:
            pass
        if min(angle_fix) < 0 and item != 'electron_pedestal_density':
            ax.axvline(x= 0, color='brown',lw=3, ls='--', label= 'outer midplane')
        else:
            pass
        
        if min(angle_fix) < -70 and item != 'electron_pedestal_density':

            ax.axvline(x= xpoint_fix, color='black',lw=3, ls='--', label= 'xpoint')
            ax.axvline(x= xpoint_fix + 360, color='black',lw=3, ls='--')
        else:
            pass
    
    def nolegend_pol_label(self, angle_fix, item, xpoint_fix, ax):
        

        if max(angle_fix) > 180 and item != 'electron_pedestal_density':
            ax.axvline(x= 180, color='gray',lw=3, ls='--')
        else:
            pass
        if min(angle_fix) < 0 and item != 'electron_pedestal_density':
            ax.axvline(x= 0, color='brown',lw=3, ls='--')
        else:
            pass
        
        if min(angle_fix) < -70 and item != 'electron_pedestal_density':
            
            ax.axvline(x= xpoint_fix, color='black',lw=3, ls='--')
            ax.axvline(x= xpoint_fix + 360, color='black',lw=3, ls='--')
        else:
            pass
        
        
        
        
    
    def neuden_poloidal_label(self, angle_fix, item, xpoint_fix):
        
        # if max(angle_fix) > 90 and item != 'electron_pedestal_density':
        #     plt.axvline(x= 90, color='black',lw=3, ls='--')
        # else:
        #     pass
        if max(angle_fix) > 180 and item != 'electron_pedestal_density':
            plt.axvline(x= 180, color='gray',lw=3, ls='--', label= 'inner midplane')
        else:
            pass
        
        # if max(angle_fix) > 240 and item != 'electron_pedestal_density':
        #     plt.axvline(x= 240, color='gray',lw=3, ls='--', label= '')
        # else:
        #     pass
        
        if min(angle_fix) < 0 and item != 'electron_pedestal_density':
            plt.axvline(x= 0, color= 'brown',lw=3, ls='--', label= 'outer midplane')
        else:
            pass
        
        if min(angle_fix) < -70 and item != 'electron_pedestal_density':
            plt.axvline(x= xpoint_fix, color='black',lw=3, ls='--', label= 'xpoint')
            plt.axvline(x= xpoint_fix + 360, color='black',lw=3, ls='--')
        else:
            pass
        
    
    
        
    def paper_polplot_method(self, log_flag, result, 
                                    i_name, ax, A_dic, color_dic, plot_order):
                    
        if log_flag:
            plt.yscale('log')
        else:
            pass
        
        
        if self.withshift == False and self.withseries == False:
            
            # result = self.data['nxny_sep_data']

            
            unit = opm.opacity_study_unit()
            pol_loc = self.data['angle']['angle_list']
            xpoint = self.data['angle']['xpoint_angle']
            a_shift = self.data['dircomp']['a_shift']
            A_val = A_dic[a_shift]
            color = color_dic[a_shift]
            
            
            self.paper_poloidal_method(item = i_name, pol_angle = pol_loc, 
                    result_dic = result, color_code = color, A_value = A_val, 
                    unit_dic = i_name, ax = ax)
            
            self.subplot_poloidal_label(angle_fix= pol_loc, item= i_name, xpoint_fix = xpoint,
                                ax = ax)
        
        elif self.withshift == True and self.withseries == False:
            
            
            # result = self.data['nxny_sep_data']
            ang_fix = self.data['angle']['angle_list']['org']
            xp_fix = self.data['angle']['xpoint_angle']['org']
            
            
            
            if i_name == 'pedestal_width':
                
                self.paper_poloidal_label(angle_fix= ang_fix, item= i_name, xpoint_fix = xp_fix,
                                    ax = ax)
                
                ax.legend(loc= 'upper center', fontsize=10)
            else:
                self.nolegend_pol_label(angle_fix= ang_fix, item= i_name, xpoint_fix = xp_fix,
                                    ax = ax)
                
                
            
            for aa in self.data['dircomp']['multi_shift']:
                
                dat_set = result[aa]
                
                unit = opm.opacity_study_unit()
                pol_loc = self.data['angle']['angle_list'][aa]
                xpoint = self.data['angle']['xpoint_angle'][aa]

                A_val = A_dic[aa]
                color = color_dic[aa]

                self.paper_poloidal_method(item = i_name, pol_angle = pol_loc, 
                             result_dic = dat_set, color_code = color, 
                             A_value = A_val, unit_dic = unit[i_name], 
                         ax = ax, plot_order = plot_order)
                
            
            
            
        
        else:
            print('sep_poloidal_plot is not there yet!')
    
    
    
    def paper_singlepolplot_method(self, log_flag, result, 
                                    i_name, A_dic, color_dic, plot_order, axs):
                    
        if log_flag:
            plt.yscale('log')
        else:
            pass
        
        
        if self.withshift == False and self.withseries == False:
            
            # result = self.data['nxny_sep_data']

            
            unit = opm.opacity_study_unit()
            pol_loc = self.data['angle']['angle_list']
            xpoint = self.data['angle']['xpoint_angle']
            a_shift = self.data['dircomp']['a_shift']
            A_val = A_dic[a_shift]
            color = color_dic[a_shift]
            
            
            self.paper_singleplot_method(item = i_name, pol_angle = pol_loc, 
                    result_dic = result, color_code = color, A_value = A_val, 
                    unit_dic = i_name, axs= axs)
            
            self.neuden_poloidal_label(angle_fix= pol_loc, item= i_name, xpoint_fix = xpoint)
        
        
        elif self.withshift == True and self.withseries == False:
            
            
            # result = self.data['nxny_sep_data']
            
            for aa in self.data['dircomp']['multi_shift']:
                
                dat_set = result[aa]
                
                unit = opm.opacity_study_unit()
                pol_loc = self.data['angle']['angle_list'][aa]
                xpoint = self.data['angle']['xpoint_angle'][aa]

                A_val = A_dic[aa]
                color = color_dic[aa]
                ang_fix = self.data['angle']['angle_list']['org']
                xp_fix = self.data['angle']['xpoint_angle']['org']
                
                
                self.paper_singleplot_method(item = i_name, pol_angle = pol_loc, 
                             result_dic = dat_set, color_code = color, 
                    A_value = A_val, unit_dic = unit[i_name], plot_order = plot_order, axs = axs)
                
            
            self.neuden_poloidal_label(angle_fix= ang_fix, item= i_name, xpoint_fix = xp_fix)
            
            
        
        else:
            print('sep_poloidal_plot is not there yet!')   
    
    
    
    def paper_poloidal_subplot(self, log_flag):
            
            itemname = ['pedestal_width', 'efold_length', 'dimensionless_opaqueness']
            # adj_list = list(result_dic.keys())
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            alphabat_list = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
            # print(adj_list)
            
            fig_n = 3
            ax_n = 1
            i_n = 0
            
            fig, axs = plt.subplots(fig_n, ax_n, sharex= True)
            for rows in range(fig_n):
                    
                result = self.data['opacity_poloidal']
                
                self.paper_polplot_method(log_flag = log_flag, 
        result = result, i_name = itemname[i_n], ax = axs[rows], 
          A_dic = A_dic, color_dic = color_dic, plot_order = alphabat_list[i_n])
                
                i_n = i_n + 1
            

            axs[fig_n -1].set_xlabel('poloidal angle')
            
            # axs[2, 1].set_xlabel('poloidal angle')
            
            plt.subplots_adjust(hspace=.0)
            # plt.tight_layout()
            fig.savefig('all_pol.pdf')
            

            fig, axs = plt.subplots()

            anchored_text = AnchoredText('{}'.format('neutral density $(n_D)_{sep}$ [$m^{-3}$]'), loc='upper center')
            
            axs.add_artist(anchored_text)

            result = self.data['opacity_poloidal']
            
            self.paper_singlepolplot_method(log_flag = log_flag, 
            result = result, i_name = 'neutral_density', 
     A_dic = A_dic, color_dic = color_dic, plot_order = alphabat_list[i_n], axs = axs)
            
            plt.legend(loc= 'lower left', fontsize=10)
            
            plt.ylim(ymin= 0, ymax = 6*10**16)
            plt.title('fitted neutral density')
            
            plt.savefig('neu_den.pdf')
            
            
            
            
        
        
            
            

            
            