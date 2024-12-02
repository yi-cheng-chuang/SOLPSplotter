# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 18:55:14 2024

@author: ychuang
"""

import SOLPS_set as sps
import SOLPSplotter_eirene as spe
import SOLPS_transcoe_adj as sta

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spe.Eirene_contour(DefaultSettings = d, loadDS = lex)

xl.load_mast_dir()
xl.load_solpsgeo()

xl.calcpsi_avcr()
xl.calc_RRsep(plotRR= False, plot_psi_dsa_align= False)
fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                          'plot_exp_and_fit': False, 'plot_shift_compare': False,
                          'data_print': True}
xl.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
xl.load_ft46()
xl.calc_sep_dsa()


poloidal_index_list = ['59']
xl.calc_dsa(pol_loc= poloidal_index_list[0])



# xl.set_plot()

xl.eirene_contour_plot()
# xl.twcontourplot(scan_style = 'tempscan', plot_option = 'Neuden contour', format_option = '1x1')

# xl.calcpsi()
# poloidal_index_list = []
# for i in range(47):
#     poloidal_index_list.append('{}'.format(25 + i))
    
# xl.opacity_data_fit(pol_list = poloidal_index_list)
# xl.calc_pol_angle(pol_list = poloidal_index_list, plot_angle= False)

# xl.opacity_poloidal_plot(log_flag = False, save_pdf = True)