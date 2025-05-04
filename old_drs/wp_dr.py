# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 19:29:40 2024

@author: ychuang
"""


import SOLPS_set as sps
import matplotlib.pyplot as plt
import withparticles as wp
import numpy as np

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = wp.withparticles(DefaultSettings = d, loadDS = lex)

xl.load_mast_dir()
xl.load_solpsgeo()
# xl.calcpsi()
xl.calcpsi_avcr()
xl.calc_RRsep(plotRR= False, plot_psi_dsa_align= False)
fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                          'plot_exp_and_fit': False, 'plot_shift_compare': False,
                          'data_print': True}
xl.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
xl.load_vessel()
xl.load_ft44()
xl.load_b2fstate()
xl.load_b2wdat()
xl.set_plot()


poloidal_index_list = []
for i in range(44):
    poloidal_index_list.append('{}'.format(25 + i))
xl.calc_pol_angle(pol_list = poloidal_index_list, plot_angle= False)
xl.calc_polsep(pol_list = poloidal_index_list)


xl.sep_length()

xl.allcover_three()

# xl.allcover_twosum()