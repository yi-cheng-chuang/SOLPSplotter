# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 15:09:47 2024

@author: ychuang
"""

import SOLPS_set as sps
import SOLPSplotter_NTplot as spn

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spn.NT_plot(DefaultSettings = d, loadDS = lex)

xl.load_mast_dir()
xl.load_solpsgeo()
xl.calcpsi_avcr()
xl.calc_RRsep(plotRR= False, plot_psi_dsa_align= False)
fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                          'plot_exp_and_fit': True, 'plot_shift_compare': False,
                          'data_print': True}
xl.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
xl.load_b2fstate()
xl.load_ft44()
xl.calc_sep_dsa()
xl.set_plot()

poloidal_index_list = ['59']
xl.calc_dsa(pol_loc= poloidal_index_list[0])



xl.opacity_data_fit(pol_list = poloidal_index_list, dat_size = 'small', check_ne = False)
xl.radial_data_fit(pol_loc = poloidal_index_list[0], dat_size = 'small', check_ne = False)


xl.neteTS_plot(scan_style = 'denscan', dat_size = 'small', xcoord_type = 'psi')


# xl.AtomNT_plot(AM_flag = 'atom', scan_style = 'denscan')