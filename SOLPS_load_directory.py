# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 11:10:50 2024

@author: ychuang
"""

import numpy as np
import SOLPS_set as sps
import load_mast_expdata_method as lmem


class load_directory:
    def __init__(self, DefaultSettings):
        
        self.DEV = DefaultSettings['DEV']
        self.a = DefaultSettings['minor_rad']
        self.withshift = DefaultSettings['withshift']
        self.withseries = DefaultSettings['withseries']
        self.terminal = DefaultSettings['terminal']
        self.series_filename = DefaultSettings['series_filename']
        self.series_tail = DefaultSettings['series_tail']
        self.series_flag = DefaultSettings['series_flag']
        self.series_compare = DefaultSettings['series_compare']
                 
            
        "DefaultSettings"    
        if isinstance(DefaultSettings, dict):
            self.DefaultSettings = DefaultSettings
        else:
            print('parameter has to be a dictionary')
        
        if DefaultSettings is None:
            print('There is no input defaultsettings')
        else:
            self.DefaultSettings = DefaultSettings
         
        keylist = []
        for key, value in self.DefaultSettings.items():
            keylist.append(key)
        
        "Useful data"
        self.data = {'defaultkey':keylist,'dircomp':{}, 'DefaultSettings': {},
                     'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
                     'gfile':{}, 'gridsettings': {}, 'psi':{}, 
                     'outputdata':{}, 'iout_data':{}}

        
    def _reset_object(self):
        # self.Shot=None
        # self.Attempts=None
        # self.Parameter=[]
        # self.PARAM={}
        # self.ExpDict={}
        # self.RadCoords={}
        self.data = {}
        
#-------------load-device-simulation-directory---------------------
    


    def load_mast_dir(self):
        if self.DEV == 'mast':
            if self.withshift == False and self.withseries == False:
                
                
                if self.terminal == False:
                    self.data['dircomp'] = sps.mast_comp_dic()
                    mast_basedir, Attempt_dic, shift_value = lmem.mast_base_dir()
                    self.data['dirdata'] = mast_basedir
                    self.data['dircomp']['Attempt'] = Attempt_dic
                    self.data['dircomp']['shift_value'] = shift_value
                
                elif self.terminal == True:
                    
                    self.data['dircomp'] = sps.mast_comp_dic()
                    mast_basedir, Attempt_dic, shift_value = lmem.terminal_single_dir()
                    self.data['dirdata'] = mast_basedir
                    self.data['dircomp']['Attempt'] = Attempt_dic
                    self.data['dircomp']['shift_value'] = shift_value
        
        
                
            elif self.withshift == True and self.withseries == False:
                
                # print('series compare is:')
                # print(self.series_compare)
                
                if self.series_compare == True:
                    
                    self.data['dircomp'] = sps.mastcomp_withshift_compare()
                    shift_dir, att_dic = lmem.mast_withshift_cp_dir()
                    self.data['dirdata'] = shift_dir
                    self.data['dircomp']['Attempt'] = att_dic
                
                else:
                    
                    self.data['dircomp'] = sps.mast_comp_dic_withshift()
                    shift_dir, att_dic = lmem.mast_withshift_dir()
                    self.data['dirdata'] = shift_dir
                    self.data['dircomp']['Attempt'] = att_dic
                    
            
            elif self.withshift == False and self.withseries == True:
                
                
                
                # series_flag = self.DefaultSettings['series_flag']
                
                if self.series_flag == 'change_den':
                    self.data['dircomp'] = sps.mast_comp_dir_series()
                    series_dir, att_dic = lmem.mast_series_dir(series_flag= self.series_flag)
                    self.data['dirdata'] = series_dir
                    self.data['dircomp']['Attempt'] = att_dic
                    
                elif self.series_flag == 'eireneN':
                    self.data['dircomp'] = sps.mast_comp_dir_eireneN()
                    series_dir, att_dic = lmem.mast_series_dir(series_flag= self.series_flag)
                    self.data['dirdata'] = series_dir
                    self.data['dircomp']['Attempt'] = att_dic
                    
                elif self.series_flag == 'change_temp':
                    self.data['dircomp'] = sps.mast_comp_dir_tempscan()
                    series_dir, att_dic = lmem.mast_series_dir(series_flag= self.series_flag)
                    self.data['dirdata'] = series_dir
                    self.data['dircomp']['Attempt'] = att_dic
                
                elif self.series_flag == 'two_compare':
                    self.data['dircomp'] = sps.mast_twocompare_dir()
                    series_dir, att_dic = lmem.series_twocompare_dir(series_flag= self.series_flag)
                    self.data['dirdata'] = series_dir
                    self.data['dircomp']['Attempt'] = att_dic
                
                
                elif self.series_flag == 'twin_scan':
                    
                    if self.terminal == True:
                        
                        self.data['dircomp'] = sps.terminal_series_comp_dir(tail = self.series_tail, 
                                            filename = self.series_filename)
                        
                        dir_comp = sps.terminal_series_comp_dir(tail = self.series_tail, 
                                            filename = self.series_filename)
                        
                        series_dir, att_dic = lmem.series_terminal_dir(series_flag = self.series_flag,
                                                                      dir_comp_dic = dir_comp)
                        self.data['dirdata'] = series_dir
                        self.data['dircomp']['Attempt'] = att_dic
                    
                    elif self.terminal == False:
                    
                        self.data['dircomp'] = sps.terminal_series_comp_dir(tail = self.series_tail, 
                                            filename = self.series_filename)
                        
                        dir_comp = sps.terminal_series_comp_dir(tail = self.series_tail, 
                                            filename = self.series_filename)
                        
                        print(dir_comp['denscan_list'])
                        
                        series_dir, att_dic = lmem.twinscan_local_dir(series_flag = self.series_flag,
                                                                      dir_comp_dic = dir_comp)
                        self.data['dirdata'] = series_dir
                        self.data['dircomp']['Attempt'] = att_dic
                        
                        
                
            elif self.withshift == True and self.withseries == True:
                print('load_mast_dir is not there yet, to be continue...')      
            else:
                print('There is a bug')

        else:
            print('DEV setting is not mast')