# -*- coding: utf-8 -*-
"""
Created on Mon May  5 14:00:23 2025

@author: ychuang
"""



from SOLPS_input.input_setting import Setting_dic, loadDS_dic
from module_load_for_plot.module_load import *



class contour_datapipline:
    
       
    def __init__(self):
        # self.publish = 'b2plottersetting'
        
        DefaultSettings = Setting_dic()
        
        self.DF = DefaultSettings
        
        
        "Useful data"
        self.data = {'dircomp':{}, 'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
                     'gfile':{}, 'gridsettings': {}, 'psi':{}, 'DefaultSettings': {}, 
                     'outputdata':{}, 'iout_data':{}}

    def run_contour_plot(self):
        
        
        xta = twscan_assist(DF = self.DF, data = self.data)
        xlp = load_prepare_module(DF = self.DF, data = self.data)
        xps = plot_setting(DF = self.DF, data = self.data)
        xfm = fit_method_collection(DF = self.DF, data = self.data)
        xlc = load_coordgeo_method(DF = self.DF, data = self.data)
        xrp = RP_mapping(DF = self.DF, data = self.data, lcm = xlc, fmc= xfm)       
        
        
        if self.DF.DEV == 'mast':
            
            
            poloidal_loc_list = xlp.prepare_for_plots(plot_type = 'contour')
            
            
            
            """
            geo_contour_theme = ['fluxexpansion', 'shift_vessel']
            
            
            """
            
            # elif geo_contour_theme == 'shift_vessel':
                
            #     xsc.shift_vessel_in_one()
            
            
            # else:
            #     pass
                       
            
            withshift = self.DF.withshift
            withseries = self.DF.withseries
            
            
            if withshift == False and withseries == True: 
                
                contour_theme = 'twscan_rectangular_contour'

        

            elif withshift == True and withseries == False:
                
                
                geo_contour_theme = 'None'
                
                if geo_contour_theme == 'shift_vessel':
                    
                    xsc.shift_vessel_in_one()
                
                
                else:
                    pass
                


if __name__ == "__main__":
    dpl = contour_datapipline()
    dpl.run_contour_plot()