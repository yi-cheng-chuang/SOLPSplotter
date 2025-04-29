# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 17:02:15 2025

@author: ychuang
"""

# from SOLPS_input.input_setting import Setting_dic, loadDS_dic
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
from load_simulation_data.read_B2simulation_data import load_B2simu_data
from load_simulation_data.read_ft_files import load_ftfiles_data
from fit_data.SOLPSplotter_fit import profile_fit
from midplane_data.midplane_netendS import midplane_radial
from midplane_data.midplane_ndS_cut import midplane_ndsource_withcut
from targets_data.target_fileload import target_dataload
from data_organize_tool.loadfiles_tools import series_loadfiles
from targets_data.target_savedata import target_radial



class load_prepare_module:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
        

    def prepare_for_plots(self, plot_type):
        
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
        xlb = load_B2simu_data(DF = self.DF, data = self.data, ldm= xldm)
        xlf = load_ftfiles_data(DF = self.DF, data = self.data, ldm= xldm)
        xpf = profile_fit(DF = self.DF, data = self.data, fmc = xfm, rp= xrp)
        xmr = midplane_radial(DF = self.DF, data = self.data, lbd = xlb, ldm = xldm)
        xmc = midplane_ndsource_withcut(DF = self.DF, data = self.data, lbd = xlb, ldm = xldm)
        xtd = target_dataload(DF = self.DF, data = self.data)
        xsl = series_loadfiles(DF = self.DF, data = self.data, td = xtd, mr = xmr)
        xtr = target_radial(DF = self.DF, data = self.data, td = xtd, sl = xsl, lbd = xlb, ldm = xldm)
        
        
        
        
        
        if self.DF.DEV == 'mast':
            
            xld.load_mast_dir()
            xlg.load_solpsgeo()
            xlg.calcpsi_avcr()
            xlg.load_vessel()
            xps.set_plot()
            
            
            xrp.calc_RRsep(plotRR= False, plot_psi_dsa_align= False, midplane_loc = 'maxis')
            fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                                      'plot_exp_and_fit': True, 'plot_shift_compare': False,
                                      'data_print': True}
            xle.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
            xlb.load_b2fstate()
            xlb.load_b2wdat()
            xlf.load_ft44()
            xrp.calc_sep_dsa()

            poloidal_index_list = ['59']
            xrp.calc_dsa(pol_loc= poloidal_index_list[0])
            
            
            if plot_type == 'radial':
                
                xpf.radial_data_fit(pol_loc = poloidal_index_list[0], check_ne = False)
                xpf.opacity_data_fit(pol_list = poloidal_index_list, check_ne = False)
                
                xmr.calc_midplane_profile()
                xmc.calc_ndSmid_cut()
                xtr.load_target_profile()
                
                
                
                return poloidal_index_list
            
            
            elif plot_type == 'poloidal':
                
                poloidal_loc_list = []
                for i in range(40):
                    poloidal_loc_list.append('{}'.format(28 + i))
                
                xpf.opacity_data_fit(pol_list = poloidal_loc_list, check_ne = False)
                xpsc.calc_pol_angle(pol_list = poloidal_loc_list, plot_angle= False)
                xtr.load_target_profile()
                

            
                return poloidal_loc_list
            
            
            elif plot_type == 'contour':
                
                xlf.load_ft46()
                
                poloidal_loc_list = []
                for i in range(40):
                    poloidal_loc_list.append('{}'.format(28 + i))
                
                xpsc.calc_pol_angle(pol_list = poloidal_loc_list, plot_angle= False)
                
                
                return poloidal_loc_list
                
                

            
            
            
            
            
            
            
            
            
            
            