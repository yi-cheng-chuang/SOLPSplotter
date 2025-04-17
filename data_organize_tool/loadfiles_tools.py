# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 17:56:31 2025

@author: ychuang
"""



from midplane_data.midplane_netendS import midplane_radial
from targets_data.target_fileload import target_dataload




class series_loadfiles:
    
    def __init__(self, DF, data, td: target_dataload, mr: midplane_radial):
        
        self.DF = DF
        self.data = data
        self.td = td
        self.mr = mr



    def oneDscan_loadfiles(self, iterlist, dat_dic, dat_struc, data_mod):
        
        
        series_compare = self.DF.series_compare
        
        if series_compare == True:
            
            print('we will improve this in the future!')
        
        else:
            
            
            for aa in iterlist:
                
                
                if data_mod == 'midplane':
                    
                    midprofiles = self.mr.calc_midplane_profile_method(itername = aa, data_struc = dat_struc)
                    
                    dat_dic[aa] = midprofiles
                
                elif data_mod == 'target':
                    
                    target_dic = self.td.tarNTdata(itername = aa, data_struc = dat_struc)
                    
                    dat_dic[aa] = target_dic
                
                
                
                
                
        return dat_dic
    
        
    def twoDscan_loadfiles(self, iterlist, iterlist_a, iterlist_b, dat_dic, dat_struc, data_mod):
        

            
        for tp in iterlist:
            aa = tp[0]
            ab = tp[1]
            
            
            if data_mod == 'midplane':
                
                midprofiles = self.mr.calc_midplane_profile_method(itername = tp, data_struc = dat_struc)
                
                dat_dic[aa][ab] = midprofiles
            
            elif data_mod == 'target':
                
                target_dic = self.td.tarNTdata(itername = tp, data_struc = dat_struc)
                
                dat_dic[aa][ab] = target_dic
            
  
        return dat_dic