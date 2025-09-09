# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 06:50:03 2025

@author: ychuang
"""

from SOLPS_input.input_setting import Setting_dic, loadDS_dic
from module_load_for_plot.module_load import *
from twscan_module.twinscan_prepare import twscan_assist
from radial_plot.SOLPSplotter_NTplot import NT_plot
from radial_plot.SOLPSplotter_ndscan import neuden_scan
from radial_plot.plot_ndmid import midnd_plot
from radial_plot.SOLPSplotter_radial import radial_plot_opacity
from radial_plot.twscan_rad import twscan_radial
from radial_plot.twscan_target import twinscan_target_radial





class module_load_tester_datapipline:
    
       
    def __init__(self):
        # self.publish = 'b2plottersetting'
        
        DefaultSettings = Setting_dic()
        
        self.DF = DefaultSettings
        
        
        "Useful data"
        self.data = {'dircomp':{}, 'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
                     'gfile':{}, 'gridsettings': {}, 'psi':{}, 'DefaultSettings': {}, 
                     'outputdata':{}, 'iout_data':{}}

    def module_tester(self):
        
        
        xta = twscan_assist(DF = self.DF, data = self.data)
        xnt = NT_plot(DF = self.DF, data = self.data, twa = xta)
        xmp = midnd_plot(DF = self.DF, data = self.data, twa = xta)
        xns = neuden_scan(DF = self.DF, data = self.data, twa = xta)
        xlp = load_prepare_module(DF = self.DF, data = self.data)
        xrpo = radial_plot_opacity(DF = self.DF, data = self.data)
        xps = plot_setting(DF = self.DF, data = self.data)
        xtr = twscan_radial(DF = self.DF, data = self.data, twa = xta)
        xttr = twinscan_target_radial(DF = self.DF, data = self.data, twa = xta)
        
        
        if self.DF.DEV == 'mast':
            
            
            poloidal_index_list = xlp.prepare_for_plots(plot_type = 'radial')
            
            
            withshift = self.DF.withshift
            withseries = self.DF.withseries
            
            
            if withshift == False and withseries == False:
                
                print('we come back to single case!')
            
            
            
            elif withshift == False and withseries == True:
                
                pass
                
                
                # dat_type = 'radial'
                
                # if dat_type == 'midplane scans':
                    
                    
                #     """
                #     radial plot for quantities at the midplane
                    
                #     """
                    
                    
                    
                    
                #     midprof = 'norm_nd'
                    
                #     if midprof == 'NT':
                        
                #         xnt.neteTS_plot(scan_style = 'tempscan', xcoord_type = 'psi')
                    
                #     elif midprof == 'nd':
                        
                #         xmp.midnd_plot(scan_style = 'denscan', xcoord_type = 'psi')
                    
                #     elif midprof == 'norm_nd':
                        
                #         xns.twinscan_ndrad_plot(scan_style = 'tempscan', log_flag = True, 
                #                                 format_option= 'neuden', xcoord_type = 'psi', withcut = False)
                    
                #     elif midprof == 'norm_S':
                        
                #         xns.twinscan_ndrad_plot(scan_style = 'tempscan', log_flag = True, 
                #                                 format_option= 'ionize', xcoord_type = 'psi', withcut = True)
                
                # elif dat_type == 'radial':
                    
                #     """
                #     radial plot
                    
                #     """
                    
                    
                #     xtr.scan_rad(format_option = 'source_peak', pol_loc = poloidal_index_list[0], 
                #                  plot_case = 'fivescan', scan_style = 'density', 
                #                  log_scale= True)
                
                
                # elif dat_type == 'targets':
                    
                    
                #     """
                #     plot ne,te,neuden,source at the target
                                      
                #     """
                    
                    
                    
                #     xttr.twinscan_targetNT(scan_style = 'denscan', log_flag = True, match= True)
                    
       
            

            
            
            
            
                
                
            
            


if __name__ == "__main__":
    dpl = module_load_tester_datapipline()
    dpl.module_tester()