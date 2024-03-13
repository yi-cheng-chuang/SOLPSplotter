# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:51:44 2024

@author: ychuang
"""

import SOLPS_set as sps
import SOLPSplotter_contour as spc
import SOLPS_transcoe_adj as sta

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spc.PlotContour(DefaultSettings = d, loadDS = lex)

xl.load_mast_dir()
xl.load_solpsgeo()
# xl.calcpsi()
xl.calcpsi_avcr()
# xl.calc_RRsep(plotRR= False, plot_psi_dsa_align= False)
fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                          'plot_exp_and_fit': False, 'plot_shift_compare': False,
                          'data_print': True}
xl.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
xl.load_vessel()

topic = 'Q1'


"""

('b2tfnb_bxuanax001.dat', False), ('hx.dat', True), ('hy.dat', True), 
('b2npc11_fnax001.dat', False), ('hz.dat', True)

"""

if topic == 'Q1':
    test_equ = 0
    if test_equ == 0:
        
        f_tuple1 = [('b2npco_sna000.dat', False), ('vol.dat', True)]
        f_tuple2 = [('b2npc11_sna001.dat', False), ('vol.dat', True)]
        f_tuple3 = [('b2npco_sna000.dat', False), ('vol.dat', True)]
        f_tuple4 = [('b2npco_sna000.dat', False), ('vol.dat', True)]
        
        f_list = [f_tuple1, f_tuple2, f_tuple3, f_tuple4]
        
        qu = xl.load_iout_ratio(file_tuple = f_list[0])
        print(qu)
        
        xl.iout_contour_plot(quant = qu)
        
        iout_quant = [ ('b2npco_sna000.dat', False), ('b2npc11_sna001.dat', False), 
                      ('b2stel_sna_ion000.dat', False), ('b2stbr_sna_eir001.dat', False), 
                      ('b2stcx_sna_001.dat', False), ('b2npc11_dnadt001.dat', False), ('b2stel_sna_rec000.dat', False)]

elif topic == 'test':
    
    iout_quant = [ ('b2npco_sna000.dat', False), ('b2npc11_sna001.dat', False), 
                  ('b2stel_sna_ion000.dat', False), ('b2stbr_sna_eir001.dat', False), 
                  ('b2stcx_sna_001.dat', False), ('b2npc11_dnadt001.dat', False), ('b2stel_sna_rec000.dat', False)]



    for item in iout_quant:
        
        qu = xl.load_iout(filename = item[0], simple_quant = item[1])
        print(qu)

        xl.iout_contour_plot(quant = qu)
    


    


