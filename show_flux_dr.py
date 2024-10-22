# -*- coding: utf-8 -*-
"""
Created on Sun Oct 20 02:06:18 2024

@author: ychuang
"""

import SOLPS_set as sps
import matplotlib.pyplot as plt
import show_flux as sf
import numpy as np

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = sf.show_flow(DefaultSettings = d, loadDS = lex)

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


# xl.totpolflux()

# xl.totpolfluxNR()


# xl.totsource()
# xl.totsourceNR()


# xl.totnd()

# xl.totndNR()


# xl.totfluxes(pol_list = poloidal_index_list)

# xl.totfluxesNR(pol_list = poloidal_index_list)

xl.totfluxesNG(pol_list = poloidal_index_list)


xl.avgfluxesNG(pol_list = poloidal_index_list)

# xl.triNR_eq(side = 'HFS')polneteNR_HFS
# xl.triNR_eq(side = 'LFS')

xl.polneteNG(side = 'HFS')
xl.polneteNG(side = 'LFS')

xl.triNG(side = 'HFS')
xl.triNG(side = 'LFS')




