# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 23:11:00 2024

@author: user
"""

import SOLPS_set as sps
import SOLPSplotter_sep_data_process as spsdp
import SOLPS_transcoe_adj as sta

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spsdp.sep_data_process(DEV = d['DEV'], withshift= d['withshift'], withseries= d['withseries'],
            DefaultSettings = d['DefaultSettings'], loadDS = lex, 
            Parameters= d['Parameters'])

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
xl.load_b2fplasmf()
# xl.load_ft46()
xl.b2fplasmf_filter()