# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 23:52:22 2025

@author: ychuang
"""


import os
import numpy as np


class Amal_directory_input:
    
    def __init__(self, DF, data):
        self.DF = DF
        self.data = data
        
    
        
    def set_wdir(self): #Function to set correct Working Directory Path depending on which machine is in use
        if os.environ['OS'] == 'Windows_NT':
            if os.environ['USERNAME'] == 'Yi-Cheng':
                basedrt = r"C:/Users/Yi-Cheng/Documents/SOLPS_Data/Simulation_Data"
                topdrt = r"C:/Users/Yi-Cheng/Documents/SOLPS_Data/Experimental_Data"


            elif os.environ['USERNAME'] == 'user':
                basedrt = r"C:/Users/user/Documents/SOLPS data/simulation data"
                topdrt = r"C:/Users/user/Documents/SOLPS data/experiment data"
                
                
            elif os.environ['USERNAME'] == 'ychuang':
                basedrt = r"C:/Users/ychuang/Documents/SOLPS_data/simulation_data"
                topdrt = r"C:/Users/ychuang/Documents/SOLPS_data/experimental_data"
               

        elif os.environ['OS'] == '5.14.0-362.24.1.el9_3.0.1.x86_64':
            if os.environ['USER'] == 'ychuang':
                basedrt = r"/sciclone/data10/ychuang/solps-iter/runs/mast"
                topdrt = r"/sciclone/data10/ychuang/solps-iter/runs/mast/gnpfiles"
               
        else:
            print('please add new directory in tools')
        
        return basedrt, topdrt



    def mastu_comp_dic(self):
        a_shift = 'org'
        shift = 0

        series = '51704_high_res'
        
        filename = '1_try'
        
        mastu_dircomp_dic = {'Shot': '51704', 'shift': shift, 
                        'series': series, 'filename': filename, 
                        'a_shift': a_shift}
        
        return mastu_dircomp_dic
    
    
    def mast_comp_dic(self):
        a_shift = 'org'
        shift_dic = {'org': 0, 'dot3': 0.3, 'dot5': 0.5, 'dot55': 0.55, 'dot7': 0.7, 'one': 1}
        
        twinscan = False
        
        if twinscan:
            series_name = 'org_cfluxb_std'
            file_name = '77_nf5.52tf4.11_save_a'
        
        else:
            series_name = 'org_cfluxb_std'
            file_name = '80_newt15_a'
            
        
        
        shift_file_dic = {'org': series_name,'dot3': 'dot3','dot5': 'dot5', 'dot55': 'dot55',
                              'dot7': 'dot7'}
        
        
        series_dic = {'org': file_name, 'dot3': '16_n900000_leakbtarnsol_dot3_a', 
            'dot5': '29_n9E5_tallies_dot5_a', 'dot55': '2_dot55_a', 'dot7': '17_n9E5_tallies_dot7_a'}
        
        outputlist = ['Output', 'Output2', 'EirOutput']
        mast_dir_dic = {'Shot': '027205', 'shift_dic': shift_dic, 
                        'shift_file_dic': shift_file_dic, 'series_dic': series_dic, 
                        'a_shift': a_shift, 'Output': outputlist}
        
        return mast_dir_dic