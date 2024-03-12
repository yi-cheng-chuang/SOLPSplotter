# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 11:10:50 2024

@author: ychuang
"""

import numpy as np
import SOLPS_set as sps
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate


class load_directory:
    def __init__(self, DefaultSettings):
        
        self.DEV = DefaultSettings['DEV']
        self.withshift = DefaultSettings['withshift']
        self.withseries = DefaultSettings['withseries']
                 
            
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
                self.data['dircomp'] = sps.mast_comp_dic()
                mast_basedir, Attempt_dic, shift_value = lmem.mast_base_dir()
                self.data['dirdata'] = mast_basedir
                self.data['dircomp']['Attempt'] = Attempt_dic
                self.data['dircomp']['shift_value'] = shift_value
                
            elif self.withshift == True and self.withseries == False:
                self.data['dircomp'] = sps.mast_comp_dic_withshift()
                shift_dir, att_dic = lmem.mast_withshift_dir()
                self.data['dirdata'] = shift_dir
                self.data['dircomp']['Attempt'] = att_dic
            
            elif self.withshift == False and self.withseries == True:
                series_flag = self.DefaultSettings['series_flag']
                if series_flag == 'change_den':
                    self.data['dircomp'] = sps.mast_comp_dir_series()
                    series_dir, att_dic = lmem.mast_series_dir(series_flag= series_flag)
                    self.data['dirdata'] = series_dir
                    self.data['dircomp']['Attempt'] = att_dic
                elif series_flag == 'eireneN':
                    self.data['dircomp'] = sps.mast_comp_dir_eireneN()
                    series_dir, att_dic = lmem.mast_series_dir(series_flag= series_flag)
                    self.data['dirdata'] = series_dir
                    self.data['dircomp']['Attempt'] = att_dic
            elif self.withshift == True and self.withseries == True:
                print('load_mast_dir is not there yet, to be continue...')      
            else:
                print('There is a bug')

        else:
            print('DEV setting is not mast')