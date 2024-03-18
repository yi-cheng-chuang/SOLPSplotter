# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 19:02:41 2024

@author: ychuang
"""

import SOLPS_set as sps
import SOLPSplotter_radial as spr
import SOLPS_transcoe_adj as sta

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spr.radial_plot(DefaultSettings = d, loadDS = lex)

xl.load_mast_dir()
xl.load_solpsgeo()
# xl.calcpsi()
xl.calcpsi_avcr()
xl.calc_RRsep(plotRR= False, plot_psi_dsa_align= False)
fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                          'plot_exp_and_fit': True, 'plot_shift_compare': False,
                          'data_print': True}
xl.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
xl.calc_sep_dsa()

poloidal_index_list = ['59']
xl.calc_dsa(pol_loc= poloidal_index_list[0])


f_tuple5 = [('hx.dat', 'hx'),('vol.dat', 'sqrt_g')]
f_tuple6 = [('hy.dat', 'hy'), ('vol.dat', 'sqrt_g')]

flux_tuple = [('b2npc11_fnax001.dat', 'flux_x_0'), ('b2npc11_fnay001.dat', 'flux_y_0')]

# flux_tuple = [('b2npc11_fnax001.dat', 'flux_x_0'), ('b2npc11_fnay001.dat', 'flux_y_0')]


f_list2 = [f_tuple5, f_tuple6]
fcoe_list = []


for tuples in f_list2:
    
    qu_dic = {}
    for aa in xl.data['dircomp']['multi_shift']:
        data, qu = xl.load_iout_ratio(file_tuple = tuples, itername = aa)
        qu_dic[aa] = data
    
    xl.data['iout_data'][qu] = qu_dic
        
    print(qu)
    fcoe_list.append(qu)
    
    
print(fcoe_list)

flux_qu_list = []

for flux_tpl in flux_tuple:
    
    flux_qu = xl.load_iout(filename = flux_tpl[0], simple_quant = flux_tpl[1])
    print(flux_qu)
    # print(flux_qu.split('_'))
    flux_qu_list.append(flux_qu)

res_qu_list = []

for item in flux_qu_list:
    
    res_qu_dic = {}
    for aa in xl.data['dircomp']['multi_shift']:
        if item.split('_')[1] == 'x':
            data, res_qu = xl.load_iout_multi(name1 = fcoe_list[0], name2 = item, 
                                input_name= 'poloidal_flux', itername = aa)
            
            res_qu_dic[aa] = data
        
        elif item.split('_')[1] == 'y':
            data, res_qu = xl.load_iout_multi(name1 = fcoe_list[1], name2 = item,
                                input_name= 'radial_flux', itername = aa)
            
            res_qu_dic[aa] = data
  
        else:
            print('a bug! flag is {}'.format(item.split('_')[1]))
    
    xl.data['iout_data'][res_qu] = res_qu_dic
    

    print(res_qu)
    res_qu_list.append(res_qu)


print(res_qu_list)

for rqu in res_qu_list:
    
    percent = False
    xpoint = False
    
    if percent:
        xl.plot_iout_radial_percent(quant = rqu, log_scale= False)
    
    elif percent == False and xpoint == False:      
        xl.plot_iout_radial_divertor(quant = rqu, log_scale= False)
    
    elif percent  == False and xpoint == True:
        
        xl.plot_iout_radial_xpoint(quant = rqu, log_scale= False)


    
xl.plot_radial_xpoint()     
        
    
    
    
    


