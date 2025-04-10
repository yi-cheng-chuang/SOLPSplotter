# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 20:22:18 2025

@author: ychuang
"""

import os
import numpy as np


class directory_input:
    
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

        series = 'mastu39_CDN'
        
        filename = '13_49404bou1e4_a'
        
        outputlist = ['Output', 'Output2', 'EirOutput']
        mastu_dircomp_dic = {'Shot': '49404', 'shift': shift, 
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


    def mast_twocompare_dir(self):
        a_shift = 'org'
        shift = 0
        series_name = 'org_cfluxb_std'
        fname_list = ['76_nf5.52tf4.11_std_a', '79_nf5.52tf4.11_save_a']
        outputlist = ['Output', 'Output2', 'EirOutput']
        mast_series_dir_dic = {'Shot': '027205', 'series_name': series_name, 'shift_value': shift,
                        'fname_list': fname_list, 'a_shift': a_shift, 'Output': outputlist}
        
        return mast_series_dir_dic




    def mast_comp_dic_withshift(self):
        
        bc = 'fixed'
        
        if bc == 'fixed':
            multi_shift = ['org', 'dot3', 'dot5', 'dot7']
            shift_dic = {'org': 0, 'dot3': 0.3, 'dot5': 0.5, 'dot7': 0.7}
            shift = ['org_cfluxb_std', 'dot3', 'dot5', 'dot7']
            tail = {'org': 'nts_a', 'dot3': 'dot3_a', 'dot5': 'dot5_a', 'dot7': 'dot7_a',
                    }
            
            series = ['76_n900000_leakbsol_nts5_a', '16_n900000_leakbtarnsol_dot3_a', '27_n900000_leakbtarnsol_dot5_a',
                      '15_n900000_leakbtarnsol_dot7_a']
            
        
        elif bc == 'fixed_new':
            multi_shift = ['org', 'dot3', 'dot5', 'dot7']
            shift_dic = {'org': 0, 'dot3': 0.3, 'dot5': 0.5, 'dot7': 0.7}
            shift = ['org_cfluxb_std', 'dot3', 'dot5', 'dot7']
            tail = {'org': 'nts_a', 'dot3': 'dot3_a', 'dot5': 'dot5_a', 'dot7': 'dot7_a',
                    }
            
            
            series = ['78_n9E5_newtallies_a', '19_n9E5_newtallies_dot3_a', '30_n9E5_newtallies_dot5_a', 
                      '18_n9E5_newtallies_dot7_a']
                
        
        elif bc == 'flux':
            multi_shift = ['org', 'dot3', 'dot5', 'dot7']
            shift_dic = {'org': 0, 'dot3': 0.3, 'dot5': 0.5, 'dot7': 0.7}
            shift = ['org_cfluxb_std', 'dot3', 'dot5', 'dot7']
            tail = {'org': 'nts_a', 'dot3': 'dot3_a', 'dot5': 'dot5_a', 'dot7': 'dot7_a',
                    }
            
            series = ['78_n5E4_fluxtallies_a', '19_n5E4_fluxtallies_dot3_a', '30_n5E4_fluxtallies_dot5_a', 
                      '18_n5E4_fluxtallies_dot7_a']
        
        elif bc == 'try':
            
            multi_shift = ['org', 'dot3', 'dot5', 'dot55', 'dot7']
            shift_dic = {'org': 0, 'dot3': 0.3, 'dot5': 0.5, 'dot55': 0.55, 'dot7': 0.7}
            shift = ['org_cfluxb_std', 'dot3', 'dot5', 'dot55', 'dot7']
            tail = {'org': 'nts_a', 'dot3': 'dot3_a', 'dot5': 'dot5_a', 'dot55': 'dot55_a', 'dot7': 'dot7_a',
                    }
            
            series = ['78_n5E4_fluxtallies_a', '19_n5E4_fluxtallies_dot3_a', '30_n5E4_fluxtallies_dot5_a', 
                      '1_dot55_a', '18_n5E4_fluxtallies_dot7_a']

        
        
        outputlist = ['Output', 'Output2', 'EirOutput']
        
        mast_withshift_dic = {'Shot': '027205', 'multi_shift': multi_shift, 'shift_dic': shift_dic, 
                              'shift_filelist': shift, 'tail': tail, 'series': series,
                              'Output': outputlist}
        
        return mast_withshift_dic




    def mastcomp_withshift_compare(self):
        

        multi_shift = ['org', 'dot3', 'dot5', 'dot7']
        shift_dic = {'org': 0, 'dot3': 0.3, 'dot5': 0.5, 'dot7': 0.7}
        shift = ['org_cfluxb_std', 'dot3', 'dot5', 'dot7']
        tail = {'org': 'nts_a', 'dot3': 'dot3_a', 'dot5': 'dot5_a', 'dot7': 'dot7_a'}
        
        series = [('78_n9E5_newtallies_a', '76_n900000_leakbsol_nts5_a'), ('19_n9E5_newtallies_dot3_a', '16_n900000_leakbtarnsol_dot3_a'),
       ('30_n9E5_newtallies_dot5_a', '27_n900000_leakbtarnsol_dot5_a'), ('18_n9E5_newtallies_dot7_a', '15_n900000_leakbtarnsol_dot7_a')]
        
        
        
        outputlist = ['Output', 'Output2', 'EirOutput']
        
        mast_withshift_dic = {'Shot': '027205', 'multi_shift': multi_shift, 'shift_dic': shift_dic, 
                              'shift_filelist': shift, 'tail': tail, 'series': series,
                              'Output': outputlist}
        
        return mast_withshift_dic


    def Ashift_dir_comp(self):
        multi_shift = ['MAST', 'D3D']
        # multi_shift = ['org', 'dot3', 'dot5', 'dot7', 'one']
        shift_dic = {'MAST': 0, 'D3D': 0.55}
        shift = ['org_new_series', 'AD3D']
        tail = {'org': 'nts5_a', 'D3D': 'd3d_a'}
        
        
        series = ['75_n900000_leakbtarnsol_nts5_a']
        
        outputlist = ['Output', 'Output2', 'EirOutput']
        
        mast_withshift_dic = {'Shot': '027205', 'multi_shift': multi_shift, 'shift_dic': shift_dic, 
                              'shift_filelist': shift, 'tail': tail, 'series': series,
                              'Output': outputlist}
        
        return mast_withshift_dic




    def mast_comp_dir_series(self):
        a_shift = 'org'
        shift = 0
        tail = '_leakbsol_nts5_a'
        outputlist = ['Output', 'Output2', 'EirOutput']
        shift_filename = 'org_denscan_fluxb_027205'
        mast_series_dir_dic = {'Shot': '027205', 'shift': shift_filename, 'shift_value': shift,
                        'tail': tail, 'a_shift': a_shift, 'Output': outputlist}
        
        return mast_series_dir_dic


    def mast_comp_dir_tempscan(self):
        a_shift = 'org'
        shift = 0
        tail = '_leakbsol_nts5_a'
        outputlist = ['Output', 'Output2', 'EirOutput']
        shift_filename = 'org_tescan_fluxb_027205'
        mast_series_dir_dic = {'Shot': '027205', 'shift': shift_filename, 'shift_value': shift,
                        'tail': tail, 'a_shift': a_shift, 'Output': outputlist}
        
        return mast_series_dir_dic

        
    def mast_comp_dir_eireneN(self):
        a_shift = 'org'
        shift = 0
        tail = '_nts5_a'
        outputlist = ['Output', 'Output2', 'EirOutput']
        shift_filename = 'org_change_particle_number'
        mast_eireneN_dir_dic = {'Shot': '027205', 'shift': shift_filename, 'shift_value': shift,
                        'tail': tail, 'a_shift': a_shift, 'Output': outputlist}
        
        return mast_eireneN_dir_dic




    def scan_list(self, denscan_dic, tempscan_dic):
        
        ds_start_num = denscan_dic['start']
        ds_stop_num = denscan_dic['stop']
        ds_space_num = denscan_dic['space']

        ds_list = np.linspace(ds_start_num, ds_stop_num, ds_space_num)

        ts_start_num = tempscan_dic['start']
        ts_stop_num = tempscan_dic['stop']
        ts_space_num = tempscan_dic['space']

        ts_list = np.linspace(ts_start_num, ts_stop_num, ts_space_num)
        
        return ds_list, ts_list





    def terminal_series_comp_dir(self):
        a_shift = 'org'
        shift = 0
        outputlist = ['Output', 'Output2', 'EirOutput']
        
        ds_dic = {'start': 5.512, 'stop': 9.512, 'space': 5}
        ts_dic = {'start': 4.115, 'stop': 8.115, 'space': 5}
        
        ds_list, ts_list = self.scan_list(denscan_dic = ds_dic, tempscan_dic = ts_dic)
        
        print('this is from sps')
        print(ds_list)
        
        
        
        mast_eireneN_dir_dic = {'Shot': '027205', 'filename': self.DF.series_filename, 'shift_value': shift,
                        'tail': self.DF.tail, 'a_shift': a_shift, 'Output': outputlist, 
                'denscan_list': ds_list, 'tempscan_list': ts_list}
        
        return mast_eireneN_dir_dic




    def mast_comp_dir_compare(self):
        
        a_shift = 'org'
        shift = 0
        
        multi_shift = ['org', 'dot3', 'dot5', 'dot7', 'one']
        shift_dic = {'org': 0, 'dot3': 0.3, 'dot5': 0.5, 'dot7': 0.7, 'one': 1}
        shift = ['org_new_series', 'dot3', 'dot5', 'dot7', 'one_LS']
        tail = {'org': 'nts_a', 'dot3': 'dot3_a', 'dot5': 'dot5_a', 'dot7': 'dot7_a',
                'one': 'one_a'}
        series = ['72_n100000_n5e3et1e2_nts5_a', '14_n100000_leakagebou_dot3_a', '25_n100000_leakagebou_dot5_a', 
                  '13_n100000_leakagebou_dot7_a', '32_n100000_leakagebou_one_a']
        
        series_2 = ['72_n100000_m12n8e3_nts5_a', '14_n100000_m12_dot3_a', '25_n100000_m12_dot5_a', 
                  '13_n100000_m12_dot7_a', '32_n100000_m12_one_a']
        
        outputlist = ['Output', 'Output2', 'EirOutput']
        shift_filename = 'org_new_series'
        mast_compare_dir_dic = {'Shot': '027205', 'shift': shift, 'shift_value': shift_dic,
                        'tail': tail, 'a_shift': a_shift, 'Output': outputlist}
        
        return mast_compare_dir_dic
    



