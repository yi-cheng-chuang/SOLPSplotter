# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 20:11:57 2024

@author: ychuang
"""

import SOLPS_set as sps
import SOLPSplotter_ioutflux as spi



d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spi.iout_flux(DefaultSettings = d, loadDS = lex)

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
xl.load_b2fplasmf()

xl.load_b2fstate_fna()

xl.load_fluxes_iout()

function_flag = 'balance'

if xl.withseries == True:
    
    if function_flag == 'two_compare':
        
        xl.plot_flux_compare()
    
    else:
        print('ioutflux_test, please check function flag')
    
    

elif xl.withseries == False:
    
    if function_flag == 'plot_flux':
        
        xl.plot_radial_particleflux()
        
        xl.radial_cell_plot()
    
    elif function_flag == 'balance':
        
        # print('balance flag is not there yet!')
        
        xl.calculate_balance()
    
    else:
        
        print('ioutflux_test, please check function flag')

else:
    
    print('ioutflux_test, there is a bug!')
        
        
        
        
    
    
        















