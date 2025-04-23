# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 15:06:54 2025

@author: ychuang
"""

from SOLPS_input.input_setting import Setting_dic, loadDS_dic
from module_load_for_plot.mastu_module_load import *





class mastu_datapipline:
    
       
    def __init__(self):
        # self.publish = 'b2plottersetting'
        
        DefaultSettings = Setting_dic()
        
        self.DF = DefaultSettings
        
        
        "Useful data"
        self.data = {'dircomp':{}, 'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
                     'gfile':{}, 'gridsettings': {}, 'psi':{}, 'DefaultSettings': {}, 
                     'outputdata':{}, 'iout_data':{}}



    def run_mastu_plot(self):
        
        xml = mastu_load_prepare_module(DF = self.DF, data = self.data)
        
        
        
        if self.DF.DEV == 'mastu':
            
            
            xml.mastu_prepare_for_plots()
            
            
            
            
            
            


if __name__ == "__main__":
    dpl = mastu_datapipline()
    dpl.run_mastu_plot()