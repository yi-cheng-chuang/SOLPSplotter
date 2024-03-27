# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 14:25:17 2024

@author: user
"""


from sep_poloidal_plot import sep_poloidal_plot
from SOLPSplotter_fit import profile_fit
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


class poster_poloidal_plot(profile_fit, sep_poloidal_plot):
    def __init__(self, DefaultSettings, loadDS):
        profile_fit.__init__(self, DefaultSettings, loadDS)
        sep_poloidal_plot.__init__(self, DefaultSettings, loadDS)
    
    
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

    def opacity_poloidal_subplot(self, log_flag, pol_index_list):
            
            itemname = ['efold_length', 'efold_length_psiN', 'pedestal_width', 
    'pedestal_width_psiN', 'dimensionless_opaqueness', 'flux_expansion']
            # adj_list = list(result_dic.keys())
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            # print(adj_list)
            
            fig_n = 3
            ax_n = 2
            i_n = 0
            
            fig, axs = plt.subplots(fig_n, ax_n)
            for rows in range(fig_n):
                
                for cols in range(ax_n):
                    
                    result = self.data['opacity_poloidal']
                    
                    self.iteminput_seppolsubplot_method(index_list = pol_index_list, 
                        log_flag = log_flag, result = result, i_name = itemname[i_n], 
                        ax = axs[rows, cols], A_dic = A_dic, color_dic = color_dic)