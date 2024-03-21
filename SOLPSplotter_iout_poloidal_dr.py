# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 11:08:33 2024

@author: ychuang
"""

import SOLPS_set as sps
import SOLPSplotter_poloidal as spp
import SOLPS_transcoe_adj as sta

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spp.poloidal_plot(DefaultSettings = d, loadDS = lex)

xl.load_mast_dir()
xl.load_solpsgeo()
# xl.calcpsi()
xl.calcpsi_avcr()
xl.calc_RRsep(plotRR= False, plot_psi_dsa_align= False)
fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                          'plot_exp_and_fit': False, 'plot_shift_compare': False,
                          'data_print': True}
xl.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
xl.load_b2fstate()
xl.load_b2fplasmf()
xl.load_output_data(param= 'Te')
xl.calc_sep_dsa()


poloidal_index_list = ['59']
xl.calc_dsa(pol_loc= poloidal_index_list[0])


xl.set_plot()
poloidal_index_list = []
for i in range(47):
    poloidal_index_list.append('{}'.format(25 + i))
    

xl.calc_pol_angle(pol_list = poloidal_index_list, plot_angle= False)

