# -*- coding: utf-8 -*-
"""
Created on Thu May 29 14:58:58 2025

@author: ychuang
"""

from SOLPS_input.input_setting import Setting_dic
from SOLPS_input.directory_input import directory_input
from load_directory.SOLPS_load_directory import load_directory
from load_directory.load_dirdata_method import load_dir_method
from load_directory.grab_attempt_number import grab_aptn_method
from load_coordinate.load_coord_method import load_coordgeo_method
from load_coordinate.SOLPSplotter_geo import load_geometry
from set_plot.plot_format import plot_setting
from fit_data.fitting_method import fit_method_collection
from midplane_data.SOLPSplotter_PRmap import RP_mapping
from separatrix_data.pol_sep_calculation import poloidal_separatrix_calc
from load_experimental_data.load_expdata_method import read_expdata_method
from load_experimental_data.SOLPSplotter_load_expdata import load_expdata
from load_experimental_data.load_piecewise_expdata import load_piecewise_expdata




class Setfixedboundary_datapipline:
    
    def __init__(self):
        # self.publish = 'b2plottersetting'
        
        DefaultSettings = Setting_dic()
        
        self.DF = DefaultSettings
        
        
        "Useful data"
        self.data = {'dircomp':{}, 'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
                     'gfile':{}, 'gridsettings': {}, 'psi':{}, 'DefaultSettings': {}, 
                     'outputdata':{}, 'iout_data':{}}



    def Set_fixedboundary(self):
        
        xdi = directory_input(DF = self.DF, data = self.data)
        xga = grab_aptn_method(DF = self.DF, data = self.data)
        xldm = load_dir_method(DF = self.DF, data = self.data, di= xdi, gam = xga)
        xld = load_directory(DF = self.DF, data = self.data, di= xdi, ldm= xldm)
        xlc = load_coordgeo_method(DF = self.DF, data = self.data)
        xlg = load_geometry(DF = self.DF, data = self.data, lcm = xlc)
        xps = plot_setting(DF = self.DF, data = self.data)
        xfm = fit_method_collection(DF = self.DF, data = self.data)
        xrp = RP_mapping(DF = self.DF, data = self.data, lcm = xlc, fmc= xfm)
        xpsc = poloidal_separatrix_calc(DF = self.DF, data = self.data, lcm = xlc, fmc= xfm)
        xre = read_expdata_method(DF = self.DF, data = self.data)
        xle = load_expdata(DF = self.DF, data = self.data, fmc= xfm, rem = xre, lg= xlg)
        xlpe = load_piecewise_expdata(DF = self.DF, data = self.data, fmc= xfm, rem = xre, lg= xlg)

        
        
        if self.DF.DEV == 'mast':
            
            xld.load_mast_dir()
            xlg.load_solpsgeo()
        
        else:
            print('we left for the future contributer to apply different profiles')
        
        xlg.calcpsi_avcr()
        xps.set_plot() 
        xps.opacity_study_unit()
        
        
        xrp.calc_RRsep(plotRR= True, plot_psi_dsa_align= False, midplane_loc = 'maxis')
        fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                                  'plot_exp_and_fit': True, 'plot_shift_compare': False,
                                  'data_print': True, 'piecewise': True}
        xlpe.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
        xrp.calc_sep_dsa()
        

        poloidal_index_list = ['56']
        xrp.calc_dsa(pol_loc= poloidal_index_list[0])
        

            
        


if __name__ == "__main__":
    pipeline = Setfixedboundary_datapipline()
    pipeline.Set_fixedboundary()