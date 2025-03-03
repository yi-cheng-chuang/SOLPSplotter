# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 15:52:51 2025

@author: ychuang
"""

import SOLPS_set as sps
import SOLPSplotter_mastuTS as spm

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spm.preprocess_mastuTS(DefaultSettings = d, loadDS = lex)

xl.load_mastu_dir()
xl.load_mastusolpsgeo()

xl.calcpsi_avcr()
xl.calc_RRsep(plotRR= True, plot_psi_dsa_align= False)
xl.plot_g()
xl.plot_sec()
xl.load_mastu_TS(plot_OD = False, plot_P = True, writefile = True)
# fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
#                           'plot_exp_and_fit': True, 'plot_shift_compare': False,
#                           'data_print': True}
# xl.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
