# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 06:38:50 2025

@author: ychuang
"""

import SOLPS_set as sps
import twscan_rad as tr
import SOLPS_transcoe_adj as sta



"""
This code plot all the radial neutral density and source for all 25 case.

"""


d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = tr.twsc_rad(DefaultSettings = d, loadDS = lex)

xl.load_mast_dir()
xl.load_solpsgeo()
# xl.calcpsi()
xl.calcpsi_avcr()
xl.calc_RRsep(plotRR= False, plot_psi_dsa_align= False)
fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                          'plot_exp_and_fit': True, 'plot_shift_compare': False,
                          'data_print': True}
xl.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
xl.load_b2fstate()
xl.load_ft44()
# xl.load_b2fplasmf()
# xl.b2fplasmf_filter()
# xl.load_output_data(param= 'Te')
xl.calc_sep_dsa()
xl.set_plot()

poloidal_index_list = ['43']
xl.calc_dsa(pol_loc= poloidal_index_list[0])


xl.opacity_data_fit(pol_list = poloidal_index_list, dat_size = 'small', check_ne = False)
xl.radial_data_fit(pol_loc = poloidal_index_list[0], dat_size = 'small', check_ne = False)

xl.load_b2wdat()

xl.scan_rad(dat_size = 'small', scan_var = 'ne_sep', format_option = 'source_peak', 
             pol_loc = poloidal_index_list[0], plot_case = 'fivescan', scan_style = 'density', 
             log_scale= True)


