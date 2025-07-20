# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 23:22:15 2025

@author: ychuang
"""

from SOLPS_input.input_setting import Setting_dic, loadDS_dic
from module_load_for_plot.Amal_mastu_module_load import *
from vesselgeo_contour_plot.SOLPSplotter_plotgeo import plot_geo
from vesselgeo_contour_plot.plan_step import step_planning
from vesselgeo_contour_plot.B2_boundary_plot import B2_boundary_contour



class Amal_mastu_datapipline:
    
       
    def __init__(self):
        # self.publish = 'b2plottersetting'
        
        DefaultSettings = Setting_dic()
        
        self.DF = DefaultSettings
        
        
        "Useful data"
        self.data = {'dircomp':{}, 'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
                     'gfile':{}, 'gridsettings': {}, 'psi':{}, 'DefaultSettings': {}, 
                     'outputdata':{}, 'iout_data':{}}



    def Amal_mastuTS_plot(self):
        
        xml = mastu_load_prepare_module(DF = self.DF, data = self.data)


        
        
        
        if self.DF.DEV == 'mastu':
            
            
            xml.mastu_prepare_for_plots()
            
            
            
                
            


if __name__ == "__main__":
    dpl = Amal_mastu_datapipline()
    dpl.Amal_mastuTS_plot()