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
from load_simulation_data.read_B2simulation_data import load_B2simu_data
from load_simulation_data.read_ft_files import load_ftfiles_data
from fit_data.SOLPSplotter_fit import profile_fit
from PRmap.midplane_profile_calculation import midplane_radial
from radial_plot.SOLPSplotter_NTplot import NT_plot
from radial_plot.plot_ndmid import midnd_plot
from twscan_module.twinscan_prepare import twscan_assist

"""
This code plot all the radial neutral density and source for all 25 case.

"""



class twscan_radial_datapipline:
    
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
        xldm = load_dir_method(DF = self.DF, data = self.data, di= xdi, gam = xga)
        xld = load_directory(DF = self.DF, data = self.data, di= xdi, ldm= xldm)
        xlc = load_coordgeo_method(DF = self.DF, data = self.data)
        xlg = load_geometry(DF = self.DF, data = self.data, lcm = xlc)
        xps = plot_setting(DF = self.DF, data = self.data)
        xfm = fit_method_collection(DF = self.DF, data = self.data)
        xrp = RP_mapping(DF = self.DF, data = self.data, lcm = xlc, fmc= xfm)
        xre = read_expdata_method(DF = self.DF, data = self.data)
        xle = load_expdata(DF = self.DF, data = self.data, fmc= xfm, rem = xre, lg= xlg)
        xlb = load_B2simu_data(DF = self.DF, data = self.data, ldm= xldm)
        xlf = load_ftfiles_data(DF = self.DF, data = self.data, ldm= xldm)
        xpf = profile_fit(DF = self.DF, data = self.data, fmc = xfm, rp= xrp)
        xta = twscan_assist(DF = self.DF, data = self.data)
        xmr = midplane_radial(DF = self.DF, data = self.data, lbd = xlb, ldm = xldm)
        xnt = NT_plot(DF = self.DF, data = self.data, twa = xta)
        xmp = midnd_plot(DF = self.DF, data = self.data, twa = xta)
        
        
        
        
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
            xlb.load_b2fstate()
            xlf.load_ft44()
            xrp.calc_sep_dsa()

            poloidal_index_list = ['40']
            xrp.calc_dsa(pol_loc= poloidal_index_list[0])


            xpf.opacity_data_fit(pol_list = poloidal_index_list, check_ne = False)
            xpf.radial_data_fit(pol_loc = poloidal_index_list[0], check_ne = False)
            
            xmr.calc_midplane_profile()
            xmr.calc_ndmid_cut()
            
            
            "midplane scans"
            
            midprof = 'NT'
            
            if midprof == 'NT':
                
                xnt.neteTS_plot(scan_style = 'denscan', xcoord_type = 'psi')
            
            elif midprof == 'nd':
                
                xmp.midnd_plot(scan_style = 'denscan', xcoord_type = 'psi')
            
            
            # xl.twinscan_ndrad_plot(scan_style = 'denscan', dat_size = 'small', log_flag = False, 
            #                        format_option= 'neuden', xcoord_type = 'psi')

            # xl.twinscan_ndrad_plot(scan_style = 'denscan', dat_size = 'small', log_flag = False, 
            #                        format_option= 'ionize', xcoord_type = 'psi')
                
                
                
                
            
            
            
            

            # xl.load_fluxes_iout()


            


if __name__ == "__main__":
    dpl = twscan_radial_datapipline()
    dpl.run_twndscan()
















