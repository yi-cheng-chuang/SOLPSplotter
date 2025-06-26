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
from poloidal_plot.twscan_boundarynendS_poloidal import twscan_boundary_nendS_polplot
from poloidal_plot.SOLPSplotter_paperpolplot import paper_poloidal_plot
from poloidal_plot.SOLPSplotter_poloidal import poloidal_plot
from poloidal_plot.grant_proposal_plot import result_explain
from poloidal_plot.integrate_flux_plot import integrate_flux
from poloidal_plot.show_flux import show_flow







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
        xtbn = twscan_boundary_nendS_polplot(DF = self.DF, data = self.data, twa = xta)
        xppp = paper_poloidal_plot(DF = self.DF, data = self.data)
        xsf = show_flow(DF = self.DF, data = self.data)
        xpp = poloidal_plot(DF = self.DF, data = self.data)
        xre = result_explain(DF = self.DF, data = self.data)
        # xif = integrate_flux(DF = self.DF, data = self.data, sc = xsc)
        
        
        
        if self.DF.DEV == 'mast':
            
            
            poloidal_loc_list = xlp.prepare_for_plots(plot_type = 'poloidal')
            
            
            """
            polplot_theme_list = ['opacity_study', 'polfluxndS']
            
            
            """
            
            withshift = self.DF.withshift
            withseries = self.DF.withseries
            
            
            if withshift == False and withseries == True:
                
                
                polplot_theme = 'opacity_study'
                
                if polplot_theme == 'polfluxndS':
                    
                    xtp.twpolfluxndS_plot(scan_style = 'denscan', pol_list = poloidal_loc_list, 
                                          log_flag = True, rad_loc= 18)
                
                elif polplot_theme == 'opacity_study':
                    
                    xto.twscan_opacity_polplot(scan_style = 'denscan', plot_option = 'opacity study poloidal plot', 
                                               format_option = '3x1', space_option = 'psiN')
                
                elif polplot_theme == 'nendS':
                    
                    xtbn.twpolfluxndS_plot(scan_style = 'denscan', pol_list = poloidal_loc_list, 
                                          log_flag = True, rad_loc= -1)
                
                
            
            
            elif withshift == True and withseries == False:
                
                
                polplot_theme = 'PSIpaper'
                
                if polplot_theme == 'PSIpaper':
                    
                    xppp.paper_poloidal_subplot(log_flag = False)
                
                
                elif polplot_theme == 'poloidal':
                    
                    # xl.neuden_data_check(pol_list= poloidal_index_list)

                    xpp.opacity_poloidal_plot(log_flag = False, save_pdf = False)
                
                
                elif polplot_theme == 'integrate_flux':

                    xif.allcover_three()

                    # xl.allcover_twosum()
                
                elif polplot_theme == 'grant':

                    xre.plot_addition_ndop()
                    
                    xre.plot_addition_ndflux()
                    
                    xre.shape_plot_addition(pol_list = poloidal_loc_list)
                
                
                elif polplot_theme == 'showflux':
                    
                                       
                    xsf.totpolflux()

                    # xl.totpolfluxNR()


                    # xl.totsource()
                    # xl.totsourceNR()


                    # xl.totnd()

                    # xl.totndNR()


                    xsf.totfluxes(pol_list = poloidal_loc_list)

                    # xl.totfluxesNR(pol_list = poloidal_index_list)

                    # xl.totfluxesNG(pol_list = poloidal_index_list)

                    # xl.rebu_fluxesNG(pol_list = poloidal_index_list)


                    # xl.avgfluxesNG(pol_list = poloidal_index_list)

                    # xl.triNR_eq(side = 'HFS')polneteNR_HFS
                    # xl.triNR_eq(side = 'LFS')

                    # xl.polneteNG(side = 'HFS')
                    # xl.polneteNG(side = 'LFS')

                    # xl.triNG(side = 'HFS')
                    # xl.triNG(side = 'LFS')

                    # xl.allcover_fhix()

                    # xl.spac_fhix()

                    # xl.hflux_bar()

                    # xl.hfluxe_bar()

                    # xl.hflux_tot_bar()

                    # xl.hflux_value_bar()

                    # xl.pflux_tot_bar()
                

                    
                    
                
            
            
            
            


            
            
            
            


if __name__ == "__main__":
    dpl = poloidal_datapipline()
    dpl.run_poloidal_plot()