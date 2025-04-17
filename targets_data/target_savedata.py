# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 19:39:59 2025

@author: ychuang
"""


import numpy as np
from targets_data.target_fileload import target_dataload
from data_organize_tool.loadfiles_tools import series_loadfiles
from load_simulation_data.read_B2simulation_data import load_B2simu_data
from load_directory.load_dirdata_method import load_dir_method



class target_radial:
    
    def __init__(self, DF, data, td: target_dataload, sl: series_loadfiles, lbd: load_B2simu_data, ldm: load_dir_method):
        
        self.DF = DF
        self.data = data
        self.sl = sl
        self.td = td
        self.lbd = lbd
        self.ldm = ldm



    def load_target_profile(self):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        
        if withshift == False and withseries == False:
            
            nx = self.data['b2fgeo']['nx']
            ny = self.data['b2fgeo']['ny']
            dat_struc = {'nx': nx, 'ny': ny}
                
            target_dic = self.td.tarNTdata(itername = None, data_struc = dat_struc)
            
            self.data['target_profile'] = target_dic
        
        elif withshift == True and withseries == False:
            
            
            target_dic = {}
            
            
            for aa in self.data['dircomp']['multi_shift']:
                
                
                nx = self.data['b2fgeo'][aa]['nx']
                ny = self.data['b2fgeo'][aa]['ny']
                dat_struc = {'nx': nx, 'ny': ny}
                
                
                tarprofiles = self.td.tarNTdata(itername = aa, data_struc = dat_struc)
                
                target_dic[aa] = tarprofiles
            
            self.data['target_profile'] = target_dic
        
        elif withshift == False and withseries == True:
            
            target_dic = {}
            
            nx = self.data['b2fgeo']['nx']
            ny = self.data['b2fgeo']['ny']
            dat_struc = {'nx': nx, 'ny': ny}
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            if self.DF.series_flag == 'twin_scan':
                
                ds_key, ts_key = self.lbd.twokeylists(printvalue= False)
                
                mid_dic = self.ldm.two_layer_dic(key_a = ds_key, key_b = ts_key)
                
                target_dic = self.sl.twoDscan_loadfiles(iterlist = scan, iterlist_a = ds_key, iterlist_b = ts_key, 
                                                        dat_dic = mid_dic, dat_struc = dat_struc, data_mod = 'target')
                
            
            else:
                
                mid_dic = {}
                target_dic = self.sl.oneDscan_loadfiles(iterlist = scan, dat_dic = mid_dic, dat_struc = dat_struc, data_mod = 'target')
            
            self.data['target_profile'] = target_dic
        
        elif withshift == True and withseries == True:
            print('calc_midplane_profile is not there yet!')
        
        
        else:
            print('calc_midplane_profile has a bug')