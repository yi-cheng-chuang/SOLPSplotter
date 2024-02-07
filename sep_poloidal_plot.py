# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 16:36:57 2024

@author: ychuang
"""

from SOLPSplotter_sep_data_process import sep_data_process
import os
import opacity_plot_method as opm
import matplotlib.pyplot as plt
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate
from scipy.optimize import curve_fit
import numpy as np


class sep_poloidal_plot(sep_data_process):
    def __init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters):
        sep_data_process.__init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters)
        
   
    
    def set_plot(self, plot_style):
        if plot_style == 'pol_subplot':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 2, markersize= 2)
            plt.rcParams.update({'font.size': 8})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
            # plt.rcParams['figure.figsize'] = 40, 12
            
  
        else:
            print('Publish setting is incorrect or add another setting')
    
    
    
    
    def sep_poloidal_plot_method(self, item, pol_angle, pol_index_list, result_dic, color_code, 
                                 A_value, unit_dic):
        
        st = int(pol_index_list[0]) 
        ed = int(pol_index_list[-1]) + 1
        
        
        plt.plot(pol_angle, result_dic[item][st:ed], '-', color= color_code ,
              label= 'aspect ratio = {}'.format(A_value))
        # plt.plot(pol_loc[aa], result_dic[i][aa], '-', color= color_dic[aa])
        plt.title('{} verses poloidal angle from {:.2f} to {:.2f}'.format(unit_dic, 
                                          max(pol_angle), min(pol_angle)))
        # plt.title('{} verses poloidal angle'.format(unit_dic[i]))
        # plt.title('{}'.format(unit_dic[i]))                   
        plt.xlabel('poloidal angle')
        plt.legend()
    
    
    def sep_poloidal_subplot_method(self, item, pol_angle, pol_index_list, result_dic, color_code, 
                                 A_value, unit_dic, ax):
        
        st = int(pol_index_list[0]) 
        ed = int(pol_index_list[-1]) + 1
        
        
        ax.plot(pol_angle, result_dic[item][st:ed], '-', color= color_code ,
              label= 'aspect ratio = {}'.format(A_value))
        # plt.plot(pol_loc[aa], result_dic[i][aa], '-', color= color_dic[aa])
        ax.set_title('{} verses poloidal angle from {:.2f} to {:.2f}'.format(unit_dic, 
                                          max(pol_angle), min(pol_angle)))
        # plt.title('{} verses poloidal angle'.format(unit_dic[i]))
        # plt.title('{}'.format(unit_dic[i]))                   
        ax.set_xlabel('poloidal angle')
        ax.legend()
    
    
        
    def poloidal_label(self, angle_fix, item, xpoint_fix):
        
        if max(angle_fix) > 90 and item != 'electron_pedestal_density' and item != 'width_relation':
            plt.axvline(x= 90, color='black',lw=3, ls='--', label= '90')
        else:
            pass
        if max(angle_fix) > 180 and item != 'electron_pedestal_density':
            plt.axvline(x= 180, color='seagreen',lw=3, ls='--', label= 'inner midplane')
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
    
    
    
    
    def subplot_poloidal_label(self, angle_fix, item, xpoint_fix, ax):
        
        if max(angle_fix) > 90 and item != 'electron_pedestal_density' and item != 'width_relation':
            ax.axvline(x= 90, color='black',lw=3, ls='--', label= '90')
        else:
            pass
        if max(angle_fix) > 180 and item != 'electron_pedestal_density':
            ax.axvline(x= 180, color='seagreen',lw=3, ls='--', label= 'inner midplane')
        else:
            pass
        if min(angle_fix) < 0 and item != 'electron_pedestal_density':
            ax.axvline(x= 0, color='darkorange',lw=3, ls='--', label= 'outer midplane')
        else:
            pass
        
        if min(angle_fix) < -70 and item != 'electron_pedestal_density':
            ax.axvline(x= xpoint_fix, color='darkblue',lw=3, ls='--', label= 'xpoint angle')
            ax.axvline(x= xpoint_fix + 360, color='darkblue',lw=3, ls='--')
        else:
            pass
        # ax.ylabel('{}'.format(unit_dic[i]))
        ax.legend()
    
    
        
    def iteminput_seppolsubplot_method(self, index_list, log_flag, result, 
                                    i_name, ax, A_dic, color_dic):
        
        # adj_list = list(result_dic.keys())
        # A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
        #           'dot7': '2.8', 'one': '3.4'}
        # color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
        #              'dot7': 'blue', 'one': 'purple'}
        
        # print(adj_list)

            
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
            
            
            self.sep_poloidal_subplot_method(item = i_name, pol_angle = pol_loc, 
    pol_index_list= index_list, result_dic = result, color_code = color, 
                                    A_value = A_val, unit_dic = i_name, ax = ax)
            
            self.subplot_poloidal_label(angle_fix= pol_loc, item= i_name, xpoint_fix = xpoint,
                                ax = ax)
        
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
                
                
                
                
                self.sep_poloidal_subplot_method(item = i_name, pol_angle = pol_loc, 
        pol_index_list= index_list, result_dic = dat_set, color_code = color, 
                                        A_value = A_val, unit_dic = i_name, ax = ax)
                
            
            self.subplot_poloidal_label(angle_fix= ang_fix, item= i_name, xpoint_fix = xp_fix,
                                ax = ax)
        
        else:
            print('sep_poloidal_plot is not there yet!')
    
    
    
    
    def iteminput_seppolplot_method(self, index_list, log_flag, result, 
                                    i_name, A_dic, color_dic):
        
                
        plt.figure(figsize=(7,7))
        if log_flag:
            plt.yscale('log')
        else:
            pass
        
        if self.withshift == False and self.withseries == False:
            
            result = self.data['nxny_sep_data']

            
            unit = opm.opacity_study_unit()
            pol_loc = self.data['angle']['angle_list']
            xpoint = self.data['angle']['xpoint_angle']
            a_shift = self.data['dircomp']['a_shift']
            A_val = A_dic[a_shift]
            color = color_dic[a_shift]
            
            
            self.sep_poloidal_plot_method(item = i_name, pol_angle = pol_loc, 
    pol_index_list= index_list, result_dic = result, color_code = color, 
                                    A_value = A_val, unit_dic = i_name)
            
            self.poloidal_label(angle_fix= pol_loc, item= i_name, xpoint_fix = xpoint)
        
        
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
                
                
                
                
                self.sep_poloidal_plot_method(item = i_name, pol_angle = pol_loc, 
        pol_index_list= index_list, result_dic = dat_set, color_code = color, 
                                        A_value = A_val, unit_dic = i_name)
                
            
            self.poloidal_label(angle_fix= ang_fix, item= i_name, xpoint_fix = xp_fix)
            
        else:
            print('sep_poloidal_plot is not there yet!')
        
        
    
    
    def itemname_method(self, shapename, itemname, it):
        
        if shapename == 'nxny':
            i_name = itemname[it]
        
        elif shapename == 'nxnyns':
            
            plasma_k = itemname[it]['item']
            ns = itemname[it]['ns']
            
            i_name = '{}%ns{}'.format(plasma_k, ns)
        
        elif shapename == 'fluxdim_ns':
            
            nf = itemname[it]['nf']
            ns = itemname[it]['ns']
            plasma_k = itemname[it]['item']
            
            i_name = '{}%nf{}%ns{}'.format(plasma_k, nf, ns)
        
        elif shapename == 'nxnycorner':
            
            i_name = itemname[it]
        
        
        elif shapename == 'nxny_corner_ns':
        
            ns = itemname[it]['ns']
            plasma_k = itemname[it]['item']
            
            i_name = '{}%ns{}'.format(plasma_k, ns)
        
        else:
            print('itemname_method is not there yet!')
        
        
        return i_name
        
        
    def set_figdir(self): #Function to set correct Working Directory Path depending on which machine is in use
        if os.environ['OS'] == 'Windows_NT':
                
            if os.environ['USERNAME'] == 'ychuang':
                
                fig_dir = r"C:\Users\ychuang\Documents\SOLPS_data\simulation_data\mast\027205\dataplot_fig"
        
        return fig_dir
    
    
    
    def sep_poloidal_subplot(self, item_name, result, shapename, index_list, log_flag):
        
        
        if self.withshift == False and self.withseries == False:
            itemname = item_name[shapename]
        
        elif self.withshift == True and self.withseries == False:
            itemname = item_name['org'][shapename]
        
        else:
            print('itername is not there yet!')
            
        
        
        # itemname = self.data['b2fplasmf_key']['nxny']
        # adj_list = list(result_dic.keys())
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        
        # print(adj_list)
        
        it = 0
        
        rows = 2
        cols = 2
        
        sub_num = int(len(itemname)/ (rows * cols))
        
        for fig_num in range(sub_num):
            
            fig, ax = plt.subplots(nrows = rows, ncols = cols)
            
            for rownum in range(rows):
                
                for colnum in range(cols):
                                       
                    i_name = self.itemname_method(shapename = shapename, 
                                                  itemname = itemname, it = it)
                    
                    
                                        
                    self.iteminput_seppolsubplot_method(index_list = index_list, 
                        log_flag = log_flag, result = result, i_name = i_name, 
                        ax = ax[rownum, colnum], A_dic = A_dic, color_dic = color_dic)
                            
                            
                    
                    it = it + 1
                    
                    
            
            plt.tight_layout()
            fig_dir  = self.set_figdir()
            plt.savefig('{}/{}_subplot_{}.png'.format(fig_dir, shapename, fig_num), dpi = 200)
            
            
        
        while it < len(itemname):
            
            
            i_name = self.itemname_method(shapename = shapename, 
                              itemname = itemname, it = it)
            
            
            self.iteminput_seppolplot_method(index_list = index_list, 
                log_flag = log_flag, result = result, i_name = i_name,
                 A_dic = A_dic, color_dic = color_dic)
            
            
            it = it + 1
            
            fig_dir  = self.set_figdir()
            plt.savefig('{}/{}_{}.png'.format(fig_dir, shapename, i_name), dpi = 200)
        
        
        
    def fplasmf_sep_process(self, datashape, pol_loc_list, b2fname):
        
        name_list = self.data['{}_key'.format(b2fname)]
        
        self.sep_data_process(specshape_list = name_list, 
                            shape_spec = datashape, b2f_name = b2fname)

        self.set_plot(plot_style = 'pol_subplot')

        self.calc_pol_angle(pol_list = pol_loc_list, plot_angle= False)

        # name_list = self.data['b2fplasmf_key']
        
        
        
        if '{}_sep_data'.format(datashape) in self.data:
            
            sep_data = self.data['{}_sep_data'.format(datashape)]
                        
            self.sep_poloidal_subplot(index_list = pol_loc_list, item_name = name_list, 
                                 shapename = datashape, result = sep_data, log_flag = False)
            
            
        else:
            
            print('no data to plot!')
            
            
            


        




# ----------------------------------------------------------------------------



"""

def sep_poloidal_subplot(self, item_name, result, shapename, index_list, log_flag):
    
    
    if self.withshift == False and self.withseries == False:
        itemname = item_name[shapename]
    
    elif self.withshift == True and self.withseries == False:
        itemname = item_name['org'][shapename]
    
    else:
        print('itername is not there yet!')
        
    
    
    # itemname = self.data['b2fplasmf_key']['nxny']
    # adj_list = list(result_dic.keys())
    A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
              'dot7': '2.8', 'one': '3.4'}
    color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                 'dot7': 'blue', 'one': 'purple'}
    
    # print(adj_list)
    
    it = 0
    
    rows = 2
    cols = 3
    
    sub_num = int(len(itemname)/ (rows * cols)) -1
    
    for fig in range(sub_num):
        
        fig, ax = plt.subplots(nrows = rows, ncols = cols)
        
        for rownum in range(rows):
            
            for colnum in range(cols):
                
                if shapename == 'nxny':
                    i_name = itemname[it]
                
                elif shapename == 'nxnyns':
                    
                    plasma_k = itemname[it]['item']
                    ns = itemname[it]['ns']
                    i_name = '{}%ns{}'.format(plasma_k, ns)
                
                
                                    
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
                    
                    
                    self.sep_poloidal_subplot_method(item = i_name, pol_angle = pol_loc, 
            pol_index_list= index_list, result_dic = result, color_code = color, 
                                            A_value = A_val, unit_dic = i_name, ax = ax[rownum, colnum])
                    
                    self.subplot_poloidal_label(angle_fix= pol_loc, item= i_name, xpoint_fix = xpoint,
                                        ax = ax[rownum, colnum])
                
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
                        
                        
                        
                        
                        self.sep_poloidal_subplot_method(item = i_name, pol_angle = pol_loc, 
                pol_index_list= index_list, result_dic = dat_set, color_code = color, 
                                                A_value = A_val, unit_dic = i_name, ax = ax[rownum, colnum])
                        
                    
                    self.subplot_poloidal_label(angle_fix= ang_fix, item= i_name, xpoint_fix = xp_fix,
                                        ax = ax[rownum, colnum])
                
                else:
                    print('haha')
                        
                        
                
                it = it + 1
        
        # plt.tight_layout()
        
    
    if it < len(itemname):
        
        
        if shapename == 'nxny':
            i_name = itemname[it]
        
        elif shapename == 'nxnyns':
            
            plasma_k = itemname[it]['item']
            ns = itemname[it]['ns']
            i_name = '{}%ns{}'.format(plasma_k, ns)
        
        
        
        plt.figure(figsize=(7,7))
        if log_flag:
            plt.yscale('log')
        else:
            pass
        
        if self.withshift == False and self.withseries == False:
            
            result = self.data['nxny_sep_data']

            
            unit = opm.opacity_study_unit()
            pol_loc = self.data['angle']['angle_list']
            xpoint = self.data['angle']['xpoint_angle']
            a_shift = self.data['dircomp']['a_shift']
            A_val = A_dic[a_shift]
            color = color_dic[a_shift]
            
            
            self.sep_poloidal_plot_method(item = i_name, pol_angle = pol_loc, 
    pol_index_list= index_list, result_dic = result, color_code = color, 
                                    A_value = A_val, unit_dic = i_name)
            
            self.poloidal_label(angle_fix= pol_loc, item= i_name, xpoint_fix = xpoint)
        
        
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
                
                
                
                
                self.sep_poloidal_plot_method(item = i_name, pol_angle = pol_loc, 
        pol_index_list= index_list, result_dic = dat_set, color_code = color, 
                                        A_value = A_val, unit_dic = i_name)
                
            
            self.poloidal_label(angle_fix= ang_fix, item= i_name, xpoint_fix = xp_fix)
            
        else:
            print('sep_poloidal_plot is not there yet!')
        
        it = it + 1



"""
            
            
        
        
    