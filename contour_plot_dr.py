# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 18:27:05 2025

@author: ychuang
"""



from SOLPS_input.input_setting import Setting_dic, loadDS_dic
from module_load_for_plot.module_load import *
from twscan_module.twinscan_prepare import twscan_assist
from contour_plot.twscan_eirene_contour import Eirene_contour
from contour_plot.geo_contour_plot import fluxexpansion_contour
from contour_plot.contourplot_toolbox import contour_plot_method_collect
from contour_plot.target_contour import target_contour_plot





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
        xec = Eirene_contour(DF = self.DF, data = self.data, twa = xta)
        xcpm = contour_plot_method_collect(DF = self.DF, data = self.data)
        xfc = fluxexpansion_contour(DF = self.DF, data = self.data, fmc = xfm, cpmc = xcpm, rp = xrp)
        xtc = target_contour_plot(DF = self.DF, data = self.data)
        
        
        if self.DF.DEV == 'mast':
            
            
            poloidal_index_list = xlp.prepare_for_plots(plot_type = 'contour')
            
            
            plotfluxexpansion = False
            
            if plotfluxexpansion:
                
                xfc.flux_expansion_contour_plot()
            
            else:
                pass
                       
            
            withshift = self.DF.withshift
            withseries = self.DF.withseries
            
            
            if withshift == False and withseries == True:
                
                
                contour_theme = 'twcontour'
                
                if contour_theme == 'twcontour':
                    
                    xec.twscan_eirene_contourplot(scan_style = 'tempscan', plot_option = 'Neuden contour', 
                                                  format_option = '1x1')
                
                else:
                    pass
                    
                    
       
            elif withshift == True and withseries == False:
                
                
                # xl.shift_vessel_in_one()


                dat_list = ['neutral density', 'Poloidal flux', 'Source', 'hx']

                xtc.iout_paper_plot(plotstyle = 'paper', dataname = 'Source', sideswitch = 'both')
                
                # xec.eirene_contour_plot()
                
                
            
            
                
                
            
            


if __name__ == "__main__":
    dpl = contour_datapipline()
    dpl.run_contour_plot()