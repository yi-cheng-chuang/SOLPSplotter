# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 10:17:17 2024

@author: ychuang
"""

"""
This code plot all the particle and heat flux scan for all 25 case across the core and separatrix

"""




import SOLPS_set as sps
import matplotlib.pyplot as plt
import target_contour as tc
import numpy as np

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = tc.target_contour(DefaultSettings = d, loadDS = lex)

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


xl.load_vessel()
xl.shift_vessel_in_one()


dat_list = ['neutral density', 'Poloidal flux', 'Source', 'hx']

xl.iout_paper_plot(plotstyle = 'paper', dataname = 'Source', sideswitch = 'both')




