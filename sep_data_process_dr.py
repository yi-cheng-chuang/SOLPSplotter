# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 23:11:00 2024

@author: user
"""

import SOLPS_set as sps
import sep_poloidal_plot as spp
import SOLPS_transcoe_adj as sta

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spp.sep_poloidal_plot(DEV = d['DEV'], withshift= d['withshift'], withseries= d['withseries'],
            DefaultSettings = d['DefaultSettings'], loadDS = lex, 
            Parameters= d['Parameters'])

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
xl.load_b2fstate()
xl.load_b2fplasmf()
# xl.load_ft46()

b2file_name = 'b2fstate'
dat_shape = 'nxnyns'

['b2fstate', 'b2fplasmf']
['nxny', 'nxnyns', 'nxnycorner', 'fluxdim_ns', 'nxny_corner_ns']


poloidal_index_list = []
for i in range(47):
    poloidal_index_list.append('{}'.format(25 + i))

b2fplasmf = xl.data[b2file_name]

xl.b2f_file_filter(b2f_file= b2fplasmf, b2f_name = b2file_name)

xl.fplasmf_sep_process(datashape = dat_shape, 
            pol_loc_list = poloidal_index_list, b2fname = b2file_name)



"""

backup:
    
    

if b2file_name == 'b2fplasmf':
    
    b2plamsf = xl.data['b2fplasmf']

    xl.b2f_file_filter(b2f_file= b2plamsf)
    
    if dat_shape == 'nxny':
        
        xl.nxny_sep_data_process()

        xl.set_plot(plot_style = 'pol_subplot')
        
            
            
        xl.calc_pol_angle(pol_list = poloidal_index_list, plot_angle= False)

        name_list = xl.data['b2fplasmf_key']
        sep_data = xl.data['nxny_sep_data']


        xl.sep_poloidal_subplot(index_list = poloidal_index_list, item_name = name_list, 
                             shapename = dat_shape, result = sep_data, log_flag = False)
        
    elif dat_shape == 'nxnyns':
        
        xl.nxny_sep_data_process()

        xl.set_plot(plot_style = 'pol_subplot')
        poloidal_index_list = []
        for i in range(47):
            poloidal_index_list.append('{}'.format(25 + i))
            
            
        xl.calc_pol_angle(pol_list = poloidal_index_list, plot_angle= False)

        name_list = xl.data['b2fplasmf_key']
        sep_data = xl.data['nxny_sep_data']


        xl.sep_poloidal_subplot(index_list = poloidal_index_list, item_name = name_list, 
                             shapename = dat_shape, result = sep_data, log_flag = False)
    
    # elif dat_shape == 'nxnyns':
    

"""
        
        
        
    
    










