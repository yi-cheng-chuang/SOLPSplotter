# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 17:10:22 2024

@author: ychuang
"""

import SOLPS_set as sps
import SOLPSplotter_radial as spr

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spr.radial_plot(DefaultSettings = d, loadDS = lex)

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
xl.calc_sep_dsa()
xl.set_plot()

poloidal_index_list = ['59']
xl.calc_dsa(pol_loc= poloidal_index_list[0])



xl.opacity_data_fit(pol_list = poloidal_index_list)
xl.radial_data_fit(pol_loc = poloidal_index_list[0])


xl.ne_te_TS_plot()

