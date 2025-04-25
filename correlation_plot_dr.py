# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 16:17:02 2025

@author: ychuang
"""



from SOLPS_input.input_setting import Setting_dic, loadDS_dic
from module_load_for_plot.module_load import *
from twscan_module.twinscan_prepare import twscan_assist
from correlation_plot.Ashift_correlation_plot import aspect_ratio_correlate
from correlation_plot.twscan_correlate import twcorrelate
from correlation_plot.twscan_bin_plot import bin_plot
from correlation_plot.twscan_core_edge import core_edge

"""
This code plot all the radial neutral density and source for all 25 case.

"""



class correlation_datapipline:
    
       
    def __init__(self):
        # self.publish = 'b2plottersetting'
        
        DefaultSettings = Setting_dic()
        
        self.DF = DefaultSettings
        
        
        "Useful data"
        self.data = {'dircomp':{}, 'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
                     'gfile':{}, 'gridsettings': {}, 'psi':{}, 'DefaultSettings': {}, 
                     'outputdata':{}, 'iout_data':{}}

    def run_correlation_plot(self):
        
        
        xta = twscan_assist(DF = self.DF, data = self.data)
        xps = plot_setting(DF = self.DF, data = self.data)
        xlp = load_prepare_module(DF = self.DF, data = self.data)
        xarc = aspect_ratio_correlate(DF = self.DF, data = self.data)
        xbp = bin_plot(DF = self.DF, data = self.data)
        xtc = twcorrelate(DF = self.DF, data = self.data, twa= xta, bp = xbp)
        xce = core_edge(DF = self.DF, data = self.data, twa= xta, bp = xbp)
        
        
        
        
        
        if self.DF.DEV == 'mast':
            
            
            poloidal_loc_list = xlp.prepare_for_plots(plot_type = 'poloidal')
            
            
            withshift = self.DF.withshift
            withseries = self.DF.withseries
            
            
            if withshift == False and withseries == True:
                
                
                corrplot_theme = 'binplot'
                
                
                if corrplot_theme == 'twscan_corr':
                    
                    xtc.twscan_corr(format_option = 'Te', pol_list = poloidal_loc_list, 
                                    plot_case = 'single', scan_style = 'density')


                
                elif corrplot_theme == 'binplot':
                    
                    xbp.bin_check(pol_list = poloidal_loc_list)

                    xbp.pol_bin_plot(pol_list = poloidal_loc_list)
                
                
                elif corrplot_theme == 'core-edge':
                    
                    xbp.bin_with_angle()

                    xce.twscan_CE(scan_var = 'ne_sep', format_option = 'Te', 
                                 pol_list = poloidal_loc_list, plot_case = 'single')
                    

                    
                    
                    
                
                
                
                
                
            
            
            elif withshift == True and withseries == False:
                
                
                xarc.plot_data_reorder(pol_list = poloidal_loc_list)
            
            
            
            
                
                
            
            


if __name__ == "__main__":
    dpl = correlation_datapipline()
    dpl.run_correlation_plot()