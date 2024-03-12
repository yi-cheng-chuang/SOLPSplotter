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
# xl.load_iout(filename = 'b2npc11_fnax001.dat')

iout_quant = [('b2npc11_fnax001.dat', False), ('hz.dat', True), ('hx.dat', True), ('hy.dat', True)]

#('b2tfnb_bxuanax001.dat', False)

for item in iout_quant:
    
    qu = xl.load_iout(filename = item[0], simple_quant = item[1])
    print(qu)

    xl.iout_contour_plot(quant = qu)
    


