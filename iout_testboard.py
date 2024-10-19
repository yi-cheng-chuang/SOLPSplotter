# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 19:38:20 2024

@author: ychuang
"""

import SOLPS_set as sps
import SOLPSplotter_iout_plotpresent as sip



d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = sip.Present_ioutplot(DefaultSettings = d, loadDS = lex)

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
xl.set_plot()
xl.load_b2fstate()

xl.calc_sep_dsa()


psi_st = 19
psi_ed = 38

psi_dic = {'st': psi_st, 'ed': psi_ed}

pol_list_a = []
for i in range(36):
    pol_list_a.append('{}'.format(28 + i))

xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)


color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
             'dot7': 'blue', 'one': 'purple'}
A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
          'dot7': '2.8', 'one': '3.4'}



['all_fluxes', 'fluxes_geo', 'fluxes_no_geo', 'neu_den', 'coe_check', 'annual_coe',
 'mag_contour', 'fluxes_no_psch', 'mag_pol', 'annual_review_mag', 'varify_vpara']


topic = 'fluxes_no_geo'


topic_label_dic = {'all_fluxes': '4', 'fluxes_geo': '2', 'fluxes_no_geo': '2', 
                   'neu_den': 'neuden_sp', 'coe_check': '4', 'annual_coe': '2',
                   }



xl.iout_plotcomb(topic_label_dic = topic_label_dic ,topic = topic, psi_dic = psi_dic, 
                  pol_list = pol_list_a, color_dic = color_dic, A_dic = A_dic)
        
 
        
            
        
        
               



