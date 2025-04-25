# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 17:37:59 2025

@author: ychuang
"""

from SOLPS_input.input_setting import Setting_dic, loadDS_dic


class test_module_load:
    
       
    def __init__(self):
        # self.publish = 'b2plottersetting'
        
        DefaultSettings = Setting_dic()
        
        self.DF = DefaultSettings
        
        
        "Useful data"
        self.data = {'dircomp':{}, 'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
                     'gfile':{}, 'gridsettings': {}, 'psi':{}, 'DefaultSettings': {}, 
                     'outputdata':{}, 'iout_data':{}}
    
    
        
        
    
    
    

if __name__ == "__main__":
    dpl = contour_datapipline()
    
    poloidal_loc_list = []
    for i in range(40):
        poloidal_loc_list.append('{}'.format(28 + i))
    
    self.DF.add_class_attribute(cls, name, value)
    
    dpl.run_contour_plot()