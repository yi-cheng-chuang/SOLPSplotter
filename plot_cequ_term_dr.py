# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 01:08:41 2024

@author: ychuang
"""

import SOLPS_set as sps
import matplotlib.pyplot as plt
import load_iout_data as lid
import load_B2_data_method as lBdm
import numpy as np
from matplotlib.offsetbox import AnchoredText
import scipy.stats as stats
from matplotlib.colors import LogNorm
from matplotlib import cm
from numpy import ma
from scipy.optimize import curve_fit



d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = lid.Plotiout(DefaultSettings = d, loadDS = lex)

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
poloidal_index_list = []
for i in range(10):
    poloidal_index_list.append('{}'.format(28 + i))

topic = 'differantiate result'

b2wdat_dic = {}

for aa in xl.data['dircomp']['multi_shift']:
    
    file_loc = '{}/'.format(xl.data['dirdata']['simudir'][aa])
    na_dat = xl.data['b2fstate']['org']['na']
    
    
    b2wdat = lBdm.read_b2wdat(b2wdatLoc = file_loc, 
                              nSpec = np.shape(na_dat)[2])
    b2wdat_dic[aa] = vars(b2wdat)   

xl.data['b2wdat'] = b2wdat_dic
    

def diff_quant_y(iout_dat):


    # del_x_Line = 
    del_y = np.zeros(np.shape(iout_dat))
    
    # print(np.shape(iout_dat))
    
    for ix in range(np.shape(iout_dat)[1]):
        for iy in range(np.shape(iout_dat)[0]):
            if (iy +1 == np.shape(iout_dat)[0]):
                del_y[iy][ix] = iout_dat[iy][ix]
            else:
                del_y[iy][ix] = iout_dat[iy][ix] - iout_dat[iy+1][ix]
    
    return del_y


def diff_quant_x(iout_dat):


    # del_x_Line = 
    del_x = np.zeros(np.shape(iout_dat))
    
    # print(np.shape(iout_dat))
    
    for ix in range(np.shape(iout_dat)[1]):
        for iy in range(np.shape(iout_dat)[0]):
            if (ix + 1 == np.shape(iout_dat)[1]):
                del_x[iy][ix] = iout_dat[iy][ix]
            else:
                del_x[iy][ix] = iout_dat[iy][ix] - iout_dat[iy][ix+1]
    
    return del_x

def sum_quant_x(iout_dat):


    # del_x_Line = 
    sum_x = np.zeros(np.shape(iout_dat))
    
    # print(np.shape(iout_dat))
    
    for ix in range(np.shape(iout_dat)[1]):
        for iy in range(np.shape(iout_dat)[0]):
            if (ix +1 == np.shape(iout_dat)[1]):
                sum_x[iy][ix] = iout_dat[iy][ix]
            else:
                sum_x[iy][ix] = iout_dat[iy][ix] + iout_dat[iy][ix+1]
    
    return sum_x


xl.flux_iout_loader()


fig, axs = plt.subplots()

color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
             'dot7': 'blue', 'one': 'purple'}
A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
          'dot7': '2.8', 'one': '3.4'}

for aa in xl.data['dircomp']['multi_shift']:
    
    pol_list_a = []
    for i in range(36):
        pol_list_a.append('{}'.format(26 + i))
    
    sk = int(pol_list_a[0])
    sd = int(pol_list_a[-1]) + 1
    
    choose_rad = 20
    
    xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
    
    ang_list = xl.data['angle']['angle_list'][aa]
    
    pnpx_dat = xl.data['iout_data']['pnpx'][aa]
    
    pnpy_dat = xl.data['iout_data']['pnpy'][aa]
    
    tcoe_dat = xl.data['iout_data']['tcoe_x'][aa]
    check_dat = np.multiply(pnpx_dat, tcoe_dat)
    
    check_list = check_dat[choose_rad, sk:sd]
    
    pnpx_list = pnpx_dat[choose_rad, sk:sd]
    
    pnpy_list = pnpy_dat[choose_rad, sk:sd]
    
    
    gradnx_dat = xl.data['iout_data']['gradn_x'][aa]
    gradnx_list = gradnx_dat[choose_rad, sk:sd]
    
    
    polgradnx_dat = xl.data['iout_data']['poloidal_gradn'][aa]
    polgradnx_list = polgradnx_dat[choose_rad, sk:sd]
    
    
    
    dat = xl.data['iout_data']['ni'][aa]
    
    hx_dat = xl.data['iout_data']['hx'][aa]
    hy_dat = xl.data['iout_data']['hy'][aa]
    
    diff_dat = diff_quant_y(iout_dat = dat)
    # try_dat = np.divide(diff_dat, hy_dat)
    
    diffx_dat = diff_quant_x(iout_dat = dat)
    
    xl.data['check'] = diff_dat
    
    diff_list = diff_dat[choose_rad -1, sk:sd]
    
    diffx_list = diffx_dat[choose_rad, sk-1:sd-1]
    
    
    sum_dat = sum_quant_x(iout_dat = dat)
    sumx_list = sum_dat[choose_rad, sk-1:sd-1]
    
    hy_list = hy_dat[choose_rad, sk-1:sd-1]
    # tryx_dat = np.divide(diffx_dat, hy_dat)
    tryx_list = np.divide(diffx_list, hy_list)
    
    pol_loc = list(map(int, pol_list_a))
    # axs.set_yscale('log')
    
    axs.plot(pol_loc, pnpx_list, linestyle = '--', 
            color= color_dic[aa])
    
    # axs.plot(pol_loc, pnpy_list, linestyle = '--', 
    #         color= color_dic[aa])
    
    # axs.plot(pol_loc, polgradnx_list, linestyle = '-', 
    #         color= color_dic[aa])
    
    # axs.plot(pol_loc, tryx_list, linestyle = '-', 
    #         color= color_dic[aa])
    
    # axs.plot(pol_loc, gradnx_list, linestyle = '--', 
    #         color= color_dic[aa])
    
    # axs.plot(pol_loc, diff_list, linestyle = '-', 
    #         color= color_dic[aa])
    
    # axs.plot(pol_loc, 250*diffx_list, linestyle = '-', 
    #         color= color_dic[aa])
    
    axs.plot(pol_loc, diffx_list, linestyle = '-', 
            color= color_dic[aa])
    
    # axs.plot(pol_loc, sumx_list, linestyle = '-', 
    #         color= color_dic[aa])
    
    
    # ratio_list = np.divide(check_list, diff_list)
    
    # axs.plot(ang_list, ratio_list, linestyle = '-', 
    #         color= color_dic[aa])
    


# b2wdats = []
# for i in range(len(caseDirs)):
    
#     na_dat = xl.data['b2fstate']['org']['na']
    
#     nSpec = np.shape(na_dat)[2]
#     b2wdats.append(read_b2wdat(rundir + caseDirs[i],nSpec))


if topic == 'differantiate result':
    
    
    if xl.withshift == True and xl.withseries == False:
        
        xl.flux_iout_loader()
    
        
        pol_list_a = []
        for i in range(32):
            pol_list_a.append('{}'.format(29 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        psi_st = 17
        psi_ed = 38
        
        psi_dic = {'st': psi_st, 'ed': psi_ed}
        
        flux_list = xl.data['iout_load_quant']['flux']
        
        quant_list = [flux_list[0], 'ni']
        
        
        fig, axs = plt.subplots(2, 1)
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        vpara_text = AnchoredText('{}'.format('$n_i$: [$m^{-3}$]'), 
                                     loc='upper center')
        
        ngvpara_text = AnchoredText('{}'.format('$\Gamma_x$: [$m^{-2}s^{-1}$]'), 
                                     loc='upper center')
        
        text_list = [ngvpara_text, vpara_text]
        
        for ind, dat_name in enumerate(quant_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                if ind == 1:
                    A_label_bol = False
                    
                    pnpx_dat = xl.data['iout_data']['pnpx'][aa]
                    
                    xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, 
                        ang_list = ang_list, psi_dic = psi_dic,
                        input_dat = pnpx_dat, art_text = text_list[ind], axs = axs[ind], 
                        color_dic = color_dic, A_dic = A_dic, nnp = 1,
                        no_A_label = A_label_bol, input_ls = '--')
                    
                else:
                    A_label_bol = True
                

                dat = xl.data['iout_data'][dat_name][aa]
                
                print(np.shape(dat.T))
                
                diff_dat = diff_quant_y(iout_dat = dat.T)
                
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, 
                    ang_list = ang_list, psi_dic = psi_dic,
                    input_dat = diff_dat.T, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, nnp = 1,
                    no_A_label = A_label_bol, input_ls = '-')
                    
            axs[ind].legend(loc= 'upper right')
        
        
        axs[1].set_xlabel('poloidal angle')
        axs[0].set_title('v parallel flux at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)








