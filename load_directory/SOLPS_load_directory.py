# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 11:10:50 2024

@author: ychuang
"""

import numpy as np
from SOLPS_input.directory_input import directory_input
from load_directory.load_dirdata_method import load_dir_method


class load_directory:
    
    
    def __init__(self, DF, data, di: directory_input, ldm: load_dir_method):
        
        
        # DefaultSettings, set_dic = Setting_dic()
        
        self.DF = DF
        self.data = data
        self.di = di
        self.ldm = ldm
        
        
        


        
    def _reset_object(self):
        # self.Shot=None
        # self.Attempts=None
        # self.Parameter=[]
        # self.PARAM={}
        # self.ExpDict={}
        # self.RadCoords={}
        self.data = {}
        
#-------------load-device-simulation-directory---------------------
    
    
    
    
    
    
    def load_cross_machine_dir(self):
        if self.DF.DEV == 'cross_machine' and self.DF.cross_series == False:
            if self.DF.Dnames == 'mast_mastu':
                
                
                
                if self.DF.terminal == False:
                    
                    
                    
                    self.data['mast_dircomp'] = self.di.mastu_comp_dic()
                    mast_basedir, Attempt_dic, shift_value = self.ldm.mastu_base_dir()
                    self.data['mast_dirdata'] = mast_basedir
                    self.data['mast_dircomp']['Attempt'] = Attempt_dic
                    self.data['mast_dircomp']['shift_value'] = shift_value
                    
                    self.data['mastu_dircomp'] = self.di.mast_comp_dic()
                    mast_basedir, Attempt_dic, shift_value = self.ldm.mast_base_dir()
                    self.data['mastu_dirdata'] = mast_basedir
                    self.data['mastu_dircomp']['Attempt'] = Attempt_dic
                    self.data['mastu_dircomp']['shift_value'] = shift_value
                    
                    
                    
                
                elif self.DF.terminal == True:
                    
                    print('load_mastu_dir function is not there yet!')
    
    
    
    
    
   
    
    def load_mastu_dir(self):
        if self.DF.DEV == 'mastu':
            if self.DF.withshift == False and self.DF.withseries == False:
                
                
                if self.DF.terminal == False:
                    self.data['dircomp'] = self.di.mastu_comp_dic()
                    mast_basedir, Attempt_dic, shift_value = self.ldm.mastu_base_dir()
                    self.data['dirdata'] = mast_basedir
                    self.data['dircomp']['Attempt'] = Attempt_dic
                    self.data['dircomp']['shift_value'] = shift_value
                
                elif self.DF.terminal == True:
                    
                    print('load_mastu_dir function is not there yet!')
    
    

    def load_mast_dir(self):
        
        
        terminal = self.DF.terminal
        
        if self.DF.DEV == 'mast':
            
            if self.DF.withshift == False and self.DF.withseries == False:
                
                
                if terminal == False:
                    self.data['dircomp'] = self.di.mast_comp_dic()
                    mast_basedir, Attempt_dic, shift_value = self.ldm.mast_base_dir()
                    self.data['dirdata'] = mast_basedir
                    self.data['dircomp']['Attempt'] = Attempt_dic
                    self.data['dircomp']['shift_value'] = shift_value
                
                elif terminal == True:
                    
                    self.data['dircomp'] = self.di.mast_comp_dic()
                    mast_basedir, Attempt_dic, shift_value = self.ldm.terminal_single_dir()
                    self.data['dirdata'] = mast_basedir
                    self.data['dircomp']['Attempt'] = Attempt_dic
                    self.data['dircomp']['shift_value'] = shift_value
        
        
                
            elif self.DF.withshift == True and self.DF.withseries == False:
                
                # print('series compare is:')
                # print(self.series_compare)
                
                if self.DF.series_compare == True:
                    
                    self.data['dircomp'] = self.di.mastcomp_withshift_compare()
                    shift_dir, att_dic = self.ldm.mast_withshift_cp_dir()
                    self.data['dirdata'] = shift_dir
                    self.data['dircomp']['Attempt'] = att_dic
                
                else:
                    
                    self.data['dircomp'] = self.di.mast_comp_dic_withshift()
                    shift_dir, att_dic = self.ldm.mast_withshift_dir()
                    self.data['dirdata'] = shift_dir
                    self.data['dircomp']['Attempt'] = att_dic
                    
            
            elif self.DF.withshift == False and self.DF.withseries == True:
                
                series_flag = self.DF.series_flag
                
                if series_flag == 'change_den':
                    self.data['dircomp'] = self.di.mast_comp_dir_series()
                    series_dir, att_dic = self.ldm.mast_series_dir()
                    self.data['dirdata'] = series_dir
                    self.data['dircomp']['Attempt'] = att_dic
                    
                elif series_flag == 'eireneN':
                    self.data['dircomp'] = self.di.mast_comp_dir_eireneN()
                    series_dir, att_dic = self.ldm.mast_series_dir()
                    self.data['dirdata'] = series_dir
                    self.data['dircomp']['Attempt'] = att_dic
                    
                elif series_flag == 'change_temp':
                    self.data['dircomp'] = self.di.mast_comp_dir_tempscan()
                    series_dir, att_dic = self.ldm.mast_series_dir()
                    self.data['dirdata'] = series_dir
                    self.data['dircomp']['Attempt'] = att_dic
                
                elif series_flag == 'two_compare':
                    self.data['dircomp'] = self.di.mast_twocompare_dir()
                    series_dir, att_dic = self.ldm.series_twocompare_dir()
                    self.data['dirdata'] = series_dir
                    self.data['dircomp']['Attempt'] = att_dic
                
                
                elif series_flag == 'twin_scan':
                    
                    if terminal == True:
                        
                        # self.data['dircomp'] = sps.terminal_series_comp_dir(tail = self.series_tail, 
                        #                     filename = self.series_filename)
                        
                        dir_comp = self.di.terminal_series_comp_dir()
                        
                        series_dir, att_dic = self.ldm.series_terminal_dir(dir_comp_dic = dir_comp)
                        
                        self.data['dircomp'] = dir_comp
                        self.data['dirdata'] = series_dir
                        self.data['dircomp']['Attempt'] = att_dic
                    
                    elif terminal == False:
                    
                        # self.data['dircomp'] = sps.terminal_series_comp_dir(tail = self.series_tail, 
                        #                     filename = self.series_filename)
                        
                        dir_comp = self.di.terminal_series_comp_dir()
                        
                        print(dir_comp['denscan_list'])
                        
                        series_dir, att_dic = self.ldm.twinscan_local_dir(dir_comp_dic = dir_comp)
                        
                        self.data['dircomp'] = dir_comp
                        self.data['dirdata'] = series_dir
                        self.data['dircomp']['Attempt'] = att_dic
                        
                        
                
            elif self.DF.withshift == True and self.DF.withseries == True:
                print('load_mast_dir is not there yet, to be continue...')      
            else:
                print('There is a bug')

        else:
            
            print('DEV setting is not mast')
        
        
        
        

"""

# self.DEV = DefaultSettings['DEV']
# self.a = DefaultSettings['minor_rad']
# self.withshift = DefaultSettings['withshift']
# self.withseries = DefaultSettings['withseries']
# self.terminal = DefaultSettings['terminal']
# self.series_filename = DefaultSettings['series_filename']
# self.series_tail = DefaultSettings['series_tail']
# self.series_flag = DefaultSettings['series_flag']
# self.series_compare = DefaultSettings['series_compare']
         
    
# "DefaultSettings"    
# if isinstance(DefaultSettings, dict):
#     self.DefaultSettings = DefaultSettings
# else:
#     print('parameter has to be a dictionary')

# if DefaultSettings is None:
#     print('There is no input defaultsettings')
# else:
#     self.DefaultSettings = DefaultSettings
 
# keylist = []
# for key, value in self.DefaultSettings.items():
#     keylist.append(key)

# "Useful data"
# self.data = {'dircomp':{}, 'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
#              'gfile':{}, 'gridsettings': {}, 'psi':{}, 
#              'outputdata':{}, 'iout_data':{}}

"""
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        