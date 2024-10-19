# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 15:51:28 2024

@author: ychuang
"""

import SOLPS_set as sps
import SOLPSplotter_twscan_iout as spt



d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spt.twinscan_ioutflux(DefaultSettings = d, loadDS = lex)

xl.load_mast_dir()
xl.load_solpsgeo()
# xl.calcpsi()
xl.calcpsi_avcr()
xl.calc_RRsep(plotRR= False, plot_psi_dsa_align= False)
fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                          'plot_exp_and_fit': False, 'plot_shift_compare': False,
                          'data_print': True}
xl.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
xl.load_output_data(param = 'IonFlx')
xl.load_vessel()
xl.load_ft44()
xl.load_b2fstate()
# xl.load_b2fplasmf()

xl.load_b2fstate_fna()
# xl.set_plot()
xl.load_fluxes_iout()
xl.twinscan_ioutplot(scan_style = 'tempscan', plot_option = 'radial flux poloidal plot', format_option = '2x2')
# a = ('b2tfnb_fnbx001.dat', 'pol_flux_no_psch')
# xl.load_iout(filename = a[0], simple_quant = a[1])

