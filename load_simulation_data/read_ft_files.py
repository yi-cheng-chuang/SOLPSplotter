# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 04:04:50 2025

@author: ychuang
"""


from load_simulation_data.load_Eirene_data_method import read_ft46, read_ft44
from load_directory.load_dirdata_method import load_dir_method


class load_ftfiles_data:
    
    def __init__(self, DF, data, ldm: load_dir_method):
        
        self.DF = DF
        self.data = data
        self.ldm = ldm
    
    
    
    """
    there are two option of fort file to load, one is fort.44.i, the other
    one is fort.46.i
          
    """
    
    def one_dim_scan_ft46(self, iterlist):
        
        ft46_dic = {}
        
        for aa in iterlist:
            
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], 'fort.46.i')
            ft46 = read_ft46(fileName = file_loc)
            ft46_dic[aa] = vars(ft46)
            
            
        return ft46_dic
        
    
    
    def two_dim_scan_ft46(self, iterlist, iterlist_a, iterlist_b):
        

        ft46_dic = self.ldm.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
        
        for tp in iterlist:
            aa = tp[0]
            ab = tp[1]
            
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa][ab], 'fort.46.i')
            ft46 = read_ft46(fileName = file_loc)
            ft46_dic[aa][ab] = vars(ft46)
        
        return ft46_dic
    
    
    
    def load_ft46(self):
        
        ftname = 'fort.46.i'
        
        if self.withshift == False and self.withseries == False:
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'], '{}'.format(ftname))
            ft46 = read_ft46(fileName = file_loc)
            ft46_dic = vars(ft46)
            
            self.data['ft46'] = ft46_dic
            # print('the next line is b2fplasmf')
            # print(type(k))
        
        elif self.withshift == True and self.withseries == False:
            ft46_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                
                file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], '{}'.format(ftname))
                ft46 = read_ft46(fileName = file_loc)
                ft46_dic[aa] = vars(ft46)
                
                
            self.data['ft46'] = ft46_dic
        
        elif self.withshift == False and self.withseries == True:
            # ft46_dic = {}
            
            # for aa in list(self.data['dircomp']['Attempt'].keys()):
                
            #     file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], '{}'.format(ftname))
            #     ft46 = read_ft46(fileName = file_loc)
            #     ft46_dic[aa] = vars(ft46)
                
                
            # self.data['ft46'] = ft46_dic
            
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            if self.series_flag == 'twin_scan':
                
                mcds = self.data['dircomp']
                
                ds_key = []
                ts_key = []
                
                print('this is mcds')
                print(type(mcds['denscan_list'][3]))
                
                for x in mcds['denscan_list']:
                    ds_key.append('{:.3f}'.format(x))
                    
                print(ds_key)
                for x in mcds['tempscan_list']:
                    ts_key.append('{:.3f}'.format(x))
                    

                print(ts_key)
                
                ft46_dic = self.two_dim_scan_ft46(iterlist = scan, iterlist_a = ds_key, 
                                       iterlist_b = ts_key)
                
                # state_dic, dim_dic = self.two_dim_scan_b2fstate(iterlist = scan, 
                #                     iterlist_a = ds_key, iterlist_b = ts_key)
                
            else:
                
                ft46_dic = self.one_dim_scan_ft46(iterlist = scan)
                # state_dic, dim_dic = self.one_dim_scan_b2fstate(iterlist = scan)
            
            
            self.data['ft46'] = ft46_dic
            
            
        else:
            print('load_b2fplasmf function is not there yet!')
    
    
    
    def one_dim_scan_ft44(self, iterlist):
        
        
        if self.series_compare == True:
            
            ft44_dic = {}
            
            for aa in iterlist:
                
                ft44_cp = {'fixed': {}, 'flux': {}}
                for kk in ['fixed', 'flux']:
                    file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa][kk], 'fort.44.i')
                    ft44 = read_ft44(fileName = file_loc)
                    ft44_cp[kk] = vars(ft44)
                
                ft44_dic[aa] = ft44_cp
                
        
        else:
            
            ft44_dic = {}
            
            for aa in iterlist:
                
                file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], 'fort.44.i')
                ft44 = read_ft44(fileName = file_loc)
                ft44_dic[aa] = vars(ft44)
        
        
        
        return ft44_dic
    
    
    def two_dim_scan_ft44(self, iterlist, iterlist_a, iterlist_b):
        


        ft44_dic = self.ldm.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
        
        for tp in iterlist:
            aa = tp[0]
            ab = tp[1]
            
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa][ab], 'fort.44.i')
            ft44 = read_ft44(fileName = file_loc)
            ft44_dic[aa][ab] = vars(ft44)
        
        return ft44_dic
     
    
    
    def load_ft44(self):
        
        ftname = 'fort.44.i'
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == False:
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'], '{}'.format(ftname))
            ft44 = read_ft44(fileName = file_loc)
            ft44_dic = vars(ft44)
            
            self.data['ft44'] = ft44_dic
            # print('the next line is b2fplasmf')
            # print(type(k))
        
        elif withshift == True and withseries == False:
            
            
            scan = self.data['dircomp']['multi_shift']
            
            ft44_dic = self.one_dim_scan_ft44(iterlist = scan)
            
            self.data['ft44'] = ft44_dic
        
        elif withshift == False and withseries == True:
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            series_flag = self.DF.series_flag
            
            
            if series_flag == 'twin_scan':
                
                mcds = self.data['dircomp']
                
                ds_key = []
                ts_key = []
                
                print('this is mcds')
                print(type(mcds['denscan_list'][3]))
                
                for x in mcds['denscan_list']:
                    ds_key.append('{:.3f}'.format(x))
                    
                print(ds_key)
                for x in mcds['tempscan_list']:
                    ts_key.append('{:.3f}'.format(x))
                    

                print(ts_key)
                
                ft44_dic = self.two_dim_scan_ft44(iterlist = scan, 
                                    iterlist_a = ds_key, iterlist_b = ts_key)
            
            else:
                ft44_dic = self.one_dim_scan_ft44(iterlist = scan)
                
                            
            self.data['ft44'] = ft44_dic
            
                
        else:
            print('load_b2fplasmf function is not there yet!')