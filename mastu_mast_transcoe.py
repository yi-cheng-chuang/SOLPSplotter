# -*- coding: utf-8 -*-
"""
Created on Sun Jul  6 15:04:50 2025

@author: ychuang
"""

from SOLPS_input.input_setting import Setting_dic
from SOLPS_input.directory_input import directory_input
from load_directory.SOLPS_load_directory import load_directory
from load_directory.load_dirdata_method import load_dir_method
from load_directory.grab_attempt_number import grab_aptn_method
from load_coordinate.load_coord_method import load_coordgeo_method
from load_coordinate.SOLPSplotter_geo import load_geometry
from load_coordinate.CM_geo import CrossMachine_load_geo
from set_plot.plot_format import plot_setting
from modify_transport_coefficient.transport_coefficient_adjust_method import transcoe_method
from modify_transport_coefficient.Transcoe_adj import transport_coefficient_adjustment
from modify_transport_coefficient.transcoe_align import transport_coefficient_alignment



class transcoe_mod_datapipline:
    
    def __init__(self):
        # self.publish = 'b2plottersetting'
        
        DefaultSettings = Setting_dic()
        
        self.DF = DefaultSettings
        
        
        "Useful data"
        
        if self.DF.DEV == 'cross_machine' and self.DF.Dnames == 'mastu_mast':
            
           
            self.data = {'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
                         'gridsettings': {}, 'psi':{}, 'iout_data':{}}
        
        else:
           
            self.data = {'dircomp':{}, 'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
                         'gfile':{}, 'gridsettings': {}, 'psi':{}, 'iout_data':{}}
            


    def transportcoe_modifier(self):
        
        xdi = directory_input(DF = self.DF, data = self.data)
        xga = grab_aptn_method(DF = self.DF, data = self.data)
        xldm = load_dir_method(DF = self.DF, data = self.data, di= xdi, gam = xga)
        xld = load_directory(DF = self.DF, data = self.data, di= xdi, ldm= xldm)
        xlc = load_coordgeo_method(DF = self.DF, data = self.data)
        xlg = load_geometry(DF = self.DF, data = self.data, lcm = xlc)
        xcmlg = CrossMachine_load_geo(DF = self.DF, data = self.data, lcm = xlc, lg = xlg)
        xps = plot_setting(DF = self.DF, data = self.data)
        xtm = transcoe_method(DF = self.DF, data = self.data, lcm = xlc)
        xtca = transport_coefficient_adjustment(DF = self.DF, data = self.data, gam = xga, tm= xtm)
        xtcal = transport_coefficient_alignment(DF = self.DF, data = self.data, gam = xga, geo = xlg, tm= xtm)
        
        
        if self.DF.DEV == 'cross_machine':
            
            if self.DF.Dnames == 'mastu_mast':
                
                xld.load_cross_machine_dir()
                xcmlg.CM_load_solpsgeo()
                xcmlg.CM_calcpsi_avcr()
                # xps.set_plot()
                
                # xtcal.transport_coe_align_plot(plot_transcoe = True, paper_transcoe = True, save_eps = False)
                xtcal.align_transco(plot_align = True, log_flag = False)

        
        # mod_detail = False
        
        # if mod_detail:
            
        #     ped_adj = [False, False, False]
        #     de_cut = [7, 32]
        #     ki_cut = [13, 18]
        #     ke_cut = [13, 21]
            
            
        #     xtca.transco_mod_detail(withmod = True, de_cut = de_cut, ki_cut = ki_cut, 
        #                 ke_cut = ke_cut, log_flag = True, modnew = True, ped_adj = ped_adj)
            
        # else:
            
        #     xtca.mod_transco(withmod = True, de_SOL = 18, ki_SOL = 18, 
        #                      ke_SOL = 21, log_flag = True, modnew = False)
        
       

        # xtcal.transport_coe_align_plot(plot_transcoe = True, paper_transcoe = True, save_eps = False)
        # xtcal.align_transco(plot_align = True, log_flag = False)

            
        


if __name__ == "__main__":
    pipeline = transcoe_mod_datapipline()
    pipeline.transportcoe_modifier()