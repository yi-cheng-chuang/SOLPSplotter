# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 17:34:38 2024

@author: ychuang
"""

import SOLPS_set as sps
import SOLPSplotter_poloidal as spp
import SOLPS_transcoe_adj as sta

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spp.poloidal_plot(DEV = d['DEV'], withshift= d['withshift'], withseries= d['withseries'],
            DefaultSettings = d['DefaultSettings'], loadDS = lex, 
            Parameters= d['Parameters'], Publish= d['Publish'])

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
# xl.b2fplasmf_filter()
xl.load_output_data(param= 'Te')
xl.calc_sep_dsa()


poloidal_index_list = ['59']
xl.calc_dsa(pol_loc= poloidal_index_list[0])



xl.set_plot()
poloidal_index_list = []
for i in range(47):
    poloidal_index_list.append('{}'.format(25 + i))
    
xl.opacity_data_fit(pol_list = poloidal_index_list)
xl.calc_pol_angle(pol_list = poloidal_index_list, plot_angle= False)

xl.opacity_poloidal_plot(log_flag = False)







