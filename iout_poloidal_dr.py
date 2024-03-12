# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 20:52:15 2024

@author: ychuang
"""

import SOLPS_set as sps
import sep_poloidal_plot as spp
import SOLPS_transcoe_adj as sta

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spp.sep_poloidal_plot(DefaultSettings = d, loadDS = lex)

xl.load_mast_dir()
xl.load_solpsgeo()
# xl.calcpsi()
xl.calcpsi_avcr()
xl.calc_RRsep(plotRR= False, plot_psi_dsa_align= False)
xl.calc_sep_dsa()
fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                          'plot_exp_and_fit': False, 'plot_shift_compare': False,
                          'data_print': True}
xl.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)

iout_quant = [('b2npc11_fnax001.dat', False), ('hz.dat', True)]


for item in iout_quant:
    
    qu = xl.load_iout(filename = item[0], simple_quant = item[1])
    print(qu)


poloidal_index_list = []
for i in range(47):
    poloidal_index_list.append('{}'.format(25 + i))





# b2fplasmf = xl.data[b2file_name]

# xl.b2f_file_filter(b2f_file= b2fplasmf, b2f_name = b2file_name)

# xl.fplasmf_sep_process(datashape = dat_shape, 
#             pol_loc_list = poloidal_index_list, b2fname = b2file_name)