# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:51:44 2024

@author: ychuang
"""

import SOLPS_set as sps
import SOLPSplotter_load_simulation_data as slsd
import SOLPS_transcoe_adj as sta

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = slsd.load_simu_data(DefaultSettings = d, loadDS = lex)

xl.load_mast_dir()
xl.load_solpsgeo()
# xl.calcpsi()
xl.calcpsi_avcr()
# xl.calc_RRsep(plotRR= False, plot_psi_dsa_align= False)
fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                          'plot_exp_and_fit': False, 'plot_shift_compare': False,
                          'data_print': True}
xl.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
xl.load_iout(filename = 'b2npc11_fnax001.dat')
