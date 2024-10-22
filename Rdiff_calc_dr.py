# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 15:38:37 2024

@author: ychuang
"""

import SOLPS_set as sps
import matplotlib.pyplot as plt
import R_diff_calc as Rdc
import numpy as np

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = Rdc.Diff_R_calc(DefaultSettings = d, loadDS = lex)

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



# xl.heat_fluxes(pol_list = poloidal_index_list)

# xl.pol_psiN(pol_list = poloidal_index_list)

xl.paper_fluxes(pol_list = poloidal_index_list)

# xl.poloidal_transcoe(pol_list = poloidal_index_list)

# xl.polNT(pol_list = poloidal_index_list)

# xl.polsource(pol_list = poloidal_index_list)

# xl.sqrtg_change(pol_list = poloidal_index_list)

# xl.R_change(pol_list = poloidal_index_list)

# xl.hxhy(pol_list = poloidal_index_list)


# xl.parallel_flow(pol_list = poloidal_index_list, side = 'outer target')


xl.pllflow_evolve(pol_list = poloidal_index_list, side = 'HFS')





# xl.term_three(pol_list = poloidal_index_list)

# xl.term_one(pol_list = poloidal_index_list)

# xl.test_case(pol_list = poloidal_index_list)


# xl.continuity_terms(pol_list = poloidal_index_list)


xl.check_targetNT(pol_list = poloidal_index_list)

# xl.check_source(pol_list = poloidal_index_list)

# xl.heatS_no_geo(pol_list = poloidal_index_list)



