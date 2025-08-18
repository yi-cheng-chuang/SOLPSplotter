# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 23:21:05 2025

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
from load_simulation_data.read_B2simulation_data import load_B2simu_data
from calc_core_boundary_method import core_boundary_calc




class mastu_mast_set_flux_boundary_datapipline:
    
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
            


    def flux_boundary_converter(self):
        
        xdi = directory_input(DF = self.DF, data = self.data)
        xga = grab_aptn_method(DF = self.DF, data = self.data)
        xldm = load_dir_method(DF = self.DF, data = self.data, di= xdi, gam = xga)
        xld = load_directory(DF = self.DF, data = self.data, di= xdi, ldm= xldm)
        xlc = load_coordgeo_method(DF = self.DF, data = self.data)
        xlg = load_geometry(DF = self.DF, data = self.data, lcm = xlc)
        xcmlg = CrossMachine_load_geo(DF = self.DF, data = self.data, lcm = xlc, lg = xlg)
        xps = plot_setting(DF = self.DF, data = self.data)
        xlbd = load_B2simu_data(DF = self.DF, data = self.data, ldm = xldm)
        xcbc = core_boundary_calc(DF = self.DF, data = self.data)
        
        
        if self.DF.DEV == 'cross_machine':
            
            if self.DF.Dnames == 'mastu_mast':
                
                xld.load_cross_machine_dir()
                xcmlg.CM_load_solpsgeo()
                xcmlg.CM_calcpsi_avcr()
                xlbd.load_cross_machine_b2fstate()
                xlbd.load_cross_mechine_b2wdat()
                # xcbc.calc_core_boundary_length(plot = False)
                xcbc.calc_core_boundary_flux()




            
        


if __name__ == "__main__":
    pipeline = mastu_mast_set_flux_boundary_datapipline()
    pipeline.flux_boundary_converter()