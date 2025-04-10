# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 01:14:04 2024

@author: ychuang
"""


from SOLPS_input.input_setting import Setting_dic, loadDS_dic
from SOLPS_input.directory_input import directory_input
from load_directory.SOLPS_load_directory import load_directory
from load_directory.load_dirdata_method import load_dir_method
from load_directory.grab_attempt_number import grab_aptn_method
from load_coordinate.load_coord_method import load_coordgeo_method
from load_coordinate.SOLPSplotter_geo import load_geometry
from set_plot.plot_format import plot_setting
from fit_data.fitting_method import fit_method_collection
from PRmap.SOLPSplotter_PRmap import RP_mapping
from load_experimental_data.load_expdata_method import read_expdata_method
from load_experimental_data.SOLPSplotter_load_expdata import load_expdata


"""
This code plot all the radial neutral density and source for all 25 case.

"""



class twndscan_datapipline:
    
    def __init__(self):
        # self.publish = 'b2plottersetting'
        
        DefaultSettings = Setting_dic()
        
        self.DF = DefaultSettings
        
        
        "Useful data"
        self.data = {'dircomp':{}, 'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
                     'gfile':{}, 'gridsettings': {}, 'psi':{}, 'DefaultSettings': {}, 
                     'outputdata':{}, 'iout_data':{}}

    def run_twndscan(self):
        
        xdi = directory_input(DF = self.DF, data = self.data)
        xga = grab_aptn_method(DF = self.DF, data = self.data)
        xld = load_dir_method(DF = self.DF, data = self.data, di= xdi, gam = xga)
        xld = load_directory(DF = self.DF, data = self.data, di= xdi, ldm= xld)
        xlc = load_coordgeo_method(DF = self.DF, data = self.data)
        xlg = load_geometry(DF = self.DF, data = self.data, lcm = xlc)
        xps = plot_setting(DF = self.DF, data = self.data)
        xfm = fit_method_collection(DF = self.DF, data = self.data)
        xrp = RP_mapping(DF = self.DF, data = self.data, lcm = xlc, fmc= xfm)
        xre = read_expdata_method(DF = self.DF, data = self.data)
        xle = load_expdata(DF = self.DF, data = self.data, fmc= xfm, rem = xre, lg= xlg)
        
        
        
        if self.DF.DEV == 'mast':
            
            xld.load_mast_dir()
            xlg.load_solpsgeo()
            xlg.calcpsi_avcr()
            xps.set_plot()
            
            
            xrp.calc_RRsep(plotRR= False, plot_psi_dsa_align= False)
            fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                                      'plot_exp_and_fit': True, 'plot_shift_compare': False,
                                      'data_print': True}
            xle.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
            xl.load_b2fstate()
            # xl.load_ft44()
            # xl.calc_sep_dsa()
            # xl.set_plot()

            # poloidal_index_list = ['40']
            # xl.calc_dsa(pol_loc= poloidal_index_list[0])


            # xl.opacity_data_fit(pol_list = poloidal_index_list, dat_size = 'small', check_ne = False)
            # xl.radial_data_fit(pol_loc = poloidal_index_list[0], dat_size = 'small', check_ne = False)

            # xl.load_fluxes_iout()


            # xl.twinscan_ndrad_plot(scan_style = 'denscan', dat_size = 'small', log_flag = False, 
            #                        format_option= 'neuden', xcoord_type = 'psi')

            # xl.twinscan_ndrad_plot(scan_style = 'denscan', dat_size = 'small', log_flag = False, 
            #                        format_option= 'ionize', xcoord_type = 'psi')
        


        
            
        


if __name__ == "__main__":
    dpl = twndscan_datapipline()
    dpl.run_twndscan()
















