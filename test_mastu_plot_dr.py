# -*- coding: utf-8 -*-
"""
Created on Tue May  6 19:31:53 2025

@author: ychuang
"""

from SOLPS_input.input_setting import Setting_dic, loadDS_dic
from module_load_for_plot.mastu_module_load import *
from vesselgeo_contour_plot.SOLPSplotter_plotgeo import plot_geo
from vesselgeo_contour_plot.plan_step import step_planning
from vesselgeo_contour_plot.B2_boundary_plot import B2_boundary_contour



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
        xpg = plot_geo(DF = self.DF, data = self.data)
        xsp = step_planning(DF = self.DF, data = self.data)
        xbbc = B2_boundary_contour(DF = self.DF, data = self.data)

        
        
        
        if self.DF.DEV == 'mastu':
            
            
            xml.mastu_prepare_for_plots()
            
            
            plot_theme = 'step_plan'
            
            if plot_theme == 'geo':
                
                             
                xpg.plot_g()
                xpg.plot_sec()
            
            elif plot_theme == 'step_plan':
                
                xsp.mastu_Ra_measure()
                xsp.vessel_change(plan_major_radius = 3.61, plan_minor_radius = 2)
            
            elif plot_theme == 'B2_boundary':
                
                
                xbbc.plot_B2boundary()
                
            


if __name__ == "__main__":
    dpl = mastu_datapipline()
    dpl.run_mastu_plot()