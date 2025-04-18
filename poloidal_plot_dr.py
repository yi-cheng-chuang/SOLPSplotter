# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 17:43:09 2025

@author: ychuang
"""


from SOLPS_input.input_setting import Setting_dic, loadDS_dic
from module_load_for_plot.module_load import *
from twscan_module.twinscan_prepare import twscan_assist
from poloidal_plot.twscan_polfluxndS_poloidal import twscan_polfluxndS_polplot
from poloidal_plot.twscan_opacity_poloidal import twscan_opacity_polplot




class poloidal_datapipline:
    
       
    def __init__(self):
        # self.publish = 'b2plottersetting'
        
        DefaultSettings = Setting_dic()
        
        self.DF = DefaultSettings
        
        
        "Useful data"
        self.data = {'dircomp':{}, 'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
                     'gfile':{}, 'gridsettings': {}, 'psi':{}, 'DefaultSettings': {}, 
                     'outputdata':{}, 'iout_data':{}}



    def run_poloidal_plot(self):
        
        xlp = load_prepare_module(DF = self.DF, data = self.data)
        xta = twscan_assist(DF = self.DF, data = self.data)
        xtp = twscan_polfluxndS_polplot(DF = self.DF, data = self.data, twa = xta)
        xto = twscan_opacity_polplot(DF = self.DF, data = self.data, twa = xta)
        
        
        
        if self.DF.DEV == 'mast':
            
            
            poloidal_loc_list = xlp.prepare_for_plots(plot_type = 'poloidal')
            
            
            """
            polplot_theme_list = ['opacity_study', 'polfluxndS']
            
            
            """
            
            withshift = self.DF.withshift
            withseries = self.DF.withseries
            
            
            if withshift == False and withseries == True:
                
                
                polplot_theme = 'polfluxndS'
                
                if polplot_theme == 'polfluxndS':
                    
                    xtp.twpolfluxndS_plot(scan_style = 'denscan', pol_list = poloidal_loc_list, 
                                          log_flag = True, rad_loc= 18)
                
                elif polplot_theme == 'opacity_study':
                    
                    xto.twscan_opacity_polplot(scan_style = 'tempscan', plot_option = 'opacity study poloidal plot', 
                                               format_option = '3x1')
                
                
            
            
            elif withshift == True and withseries == False:
                
                pass
            
            
            
            


            # xl.neuden_data_check(pol_list= poloidal_index_list)

            # xl.opacity_poloidal_plot(log_flag = False, save_pdf = False)
            
            
            


if __name__ == "__main__":
    dpl = poloidal_datapipline()
    dpl.run_poloidal_plot()