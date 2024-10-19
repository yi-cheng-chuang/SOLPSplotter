# -*- coding: utf-8 -*-
"""
Created on Thu May 23 12:45:06 2024

@author: ychuang
"""

import SOLPS_set as sps
import matplotlib.pyplot as plt
import SOLPSplotter_contour as spc
import SOLPS_transcoe_adj as sta
import numpy as np
import load_iout_data as lid
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


xl.calc_pol_angle(pol_list = poloidal_index_list, plot_angle= False)

# xl.opacity_poloidal_plot(log_flag = False, save_pdf = False)

# xl.neuden_percent()


['all_fluxes', 'fluxes_no_geo', 'fluxes_geo', 'neu_den', 'annual_coe',
 'coe_check', 'mag_contour', 'fluxes_no_psch', 'mag_pol', 'annual_review_mag', 
 'varify_vpara']


topic = 'all_fluxes'





"""
Retire functions:
                    
if topic == 'all_fluxes':
    
    if xl.withshift == True and xl.withseries == False:
        
        xl.flux_iout_loader()
    
        psi_st = 19
        psi_ed = 38
        
        psi_dic = {'st': psi_st, 'ed': psi_ed}
        
        
        pol_list_a = []
        for i in range(36):
            pol_list_a.append('{}'.format(28 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        
        flux_list = ['poloidal_flux','radial_flux', 'flux_x_0', 'flux_y_0']
        
        
        fig, axs = plt.subplots(4, 1)
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        pol_text = AnchoredText('{}'.format('Poloidal flux $\Gamma_x$ [$m^{-1} s^{-1}$]'), 
                                     loc='upper center')
        
        rad_text = AnchoredText('{}'.format('Radial flux $\Gamma_y$ [$m^{-1} s^{-1}$]'), 
                                     loc='upper center')
        
        x0_text = AnchoredText('{}'.format('Poloidal flux $\sqrt{g}/h_x*\Gamma_x$ [$m*s^{-1}$]'), 
                                     loc='lower center')
        
        y0_text = AnchoredText('{}'.format('radial flux $\sqrt{g}/h_y*\Gamma_y$ [$m*s^{-1}$]'), 
                                     loc='upper center')
        
        text_list = [pol_text, rad_text, x0_text, y0_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                if ind == 3:
                    A_label_bol = False
                else:
                    A_label_bol = True
                
                in_dat = xl.data['iout_data'][dat_name][aa]
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    art_text = text_list[ind], axs = axs[ind], nnp = 1, input_dat = in_dat, input_ls = '-',
                    color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = A_label_bol)
                
            if ind % 2 == 0:
                
                axs[ind].axhline(y=0, color = 'black', linestyle = '--', label= '$\Gamma_x$ = 0')
            
            else:
                pass
            
            axs[ind].legend(loc= 'upper right')
        
        axs[3].set_xlabel('poloidal angle')
        axs[0].set_title('Particle flux at separatrix')
        
        
        plt.subplots_adjust(hspace=.0)



if topic == 'fluxes_geo':
    
    if xl.withshift == True and xl.withseries == False:
        
        xl.flux_iout_loader()
    
        psi_st = 19
        psi_ed = 38
        
        psi_dic = {'st': psi_st, 'ed': psi_ed}
        
        
        pol_list_a = []
        for i in range(36):
            pol_list_a.append('{}'.format(28 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        

        fig, axs = plt.subplots(2, 1)
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        pol_text = AnchoredText('{}'.format('Poloidal flux $\Gamma_x$ [$m^{-1} s^{-1}$]'), 
                                     loc='upper center')
        
        neu_text = AnchoredText('{}'.format('Neutral density [$m^{-3}$]'), 
                                     loc='upper center')
        
        rad_text = AnchoredText('{}'.format('Radial flux $\Gamma_y$ [$m^{-1} s^{-1}$]'), 
                                     loc='upper center')
        

        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
            i_name = 'flux_x_0', art_text = pol_text, axs = axs[0], 
            color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = True)
            
        axs[0].axhline(y=0, color = 'black', linestyle = '--', label= '$\Gamma_x$ = 0')
        
        
        
        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    i_name = 'flux_y_0', art_text = rad_text, axs = axs[1], 
            color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = False)
            
            
        axs[1].set_xlabel('poloidal angle')
        axs[1].legend(loc= 'upper right')
        axs[0].legend(loc= 'upper right')
        axs[0].set_title('Particle flux at separatrix')
        
        
        plt.subplots_adjust(hspace=.0)



if topic == 'fluxes_no_geo':
    
    if xl.withshift == True and xl.withseries == False:
        
        xl.flux_iout_loader()
    
        psi_st = 19
        psi_ed = 38
        
        psi_dic = {'st': psi_st, 'ed': psi_ed}
        
        
        pol_list_a = []
        for i in range(36):
            pol_list_a.append('{}'.format(28 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        

        fig, axs = plt.subplots(2, 1)
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        pol_text = AnchoredText('{}'.format('(a) Poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                     loc='upper center')
        
        neu_text = AnchoredText('{}'.format('Neutral density [$m^{-3}$]'), 
                                     loc='upper center')
        
        rad_text = AnchoredText('{}'.format('(b) Radial flux $\Gamma_r$ [$m^{-2} s^{-1}$]'), 
                                     loc='upper center')
        

        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            pflux_dat = xl.data['iout_data']['poloidal_flux'][aa]
            
            xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
            input_dat = pflux_dat, art_text = pol_text, axs = axs[0], 
            color_dic = color_dic, A_dic = A_dic, no_A_label = True, input_ls = '-')
            
        axs[0].axhline(y=0, color = 'black', linestyle = '--', label= '$\Gamma_{\Theta}$ = 0')
        
        
        
        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            rflux_dat = xl.data['iout_data']['radial_flux'][aa]
            
            xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = rflux_dat, art_text = rad_text, axs = axs[1], 
            color_dic = color_dic, A_dic = A_dic, no_A_label = False, input_ls = '-')
            
            
        axs[1].set_xlabel('poloidal angle')
        axs[1].legend(loc= 'upper right')
        axs[0].legend(loc= 'upper right')
        axs[0].set_title('Particle flux at separatrix')
        
        
        plt.subplots_adjust(hspace=.0)





if topic == 'neu_den':
    
    if xl.withshift == True and xl.withseries == False:
        
        pol_list_a = []
        for i in range(36):
            pol_list_a.append('{}'.format(26 + i))
        
        neu_text = AnchoredText('{}'.format('Neutral density [$m^{-3}$]'), 
                                     loc='upper center')


        for aa in xl.data['dircomp']['multi_shift']:
            
            neuden_data_a = []
            neuden_data_b = []
            
            for kt in pol_list_a:
                
                neuden_data = xl.data['ft44'][aa]['dab2'][int(kt), psi_st:psi_ed]
                
                neuden_data_a.append(neuden_data.max())
                neuden_data_b.append(neuden_data.min())
            
            sk = int(pol_list_a[0])
            sd = int(pol_list_a[-1]) + 1
            
            ang_list = xl.data['angle']['angle_list'][aa]
        
        
            neuden_dat = np.transpose(xl.data['ft44'][aa]['dab2'][sk:sd, psi_st:psi_ed, 0])
            
            axs[1].add_artist(neu_text)
            
            # axs[1].fill_between(ang_list, neuden_data_a, neuden_data_b, 
            #                  color= color_dic[aa], alpha = 0.4)
            
            axs[1].plot(ang_list, neuden_dat[0, :], linestyle='-', color= color_dic[aa])
            
            # axs[1].plot(ang_list, neuden_dat[-1, :], '-', color= color_dic[aa])



if topic == 'coe_check':
    
    if xl.withshift == True and xl.withseries == False:
        
        xl.flux_iout_loader()
    
        psi_st = 19
        psi_ed = 38
        
        psi_dic = {'st': psi_st, 'ed': psi_ed}
        
        
        pol_list_a = []
        for i in range(36):
            pol_list_a.append('{}'.format(28 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        
        flux_list = ['hx','hx_divide_sqrt_g', 'hy','hy_divide_sqrt_g']
        
        
        fig, axs = plt.subplots(4, 1)
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        hx_text = AnchoredText('{}'.format('$h_x$: [m]'), 
                                     loc='upper center')
        
        hxg_text = AnchoredText('{}'.format('$h_x/ \sqrt{g}$: [$m^{-2}$]'), 
                                     loc='upper center')
        
        hy_text = AnchoredText('{}'.format('$h_y$: [m]'), 
                                     loc='upper center')
        
        hyg_text = AnchoredText('{}'.format('$h_y/ \sqrt{g}$: [$m^{-2}$]'), 
                                     loc='upper center')
        
        text_list = [hx_text, hxg_text, hy_text, hyg_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                if ind == 3:
                    A_label_bol = False
                else:
                    A_label_bol = True
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    i_name = dat_name, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = A_label_bol)
                    
            
            axs[ind].legend(loc= 'upper right')
        
        axs[3].set_xlabel('poloidal angle')
        axs[0].set_title('coefficients at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)



if topic == 'annual_coe':
    
    if xl.withshift == True and xl.withseries == False:
        
        xl.flux_iout_loader()
    
        psi_st = 19
        psi_ed = 38
        
        psi_dic = {'st': psi_st, 'ed': psi_ed}
        
        
        pol_list_a = []
        for i in range(36):
            pol_list_a.append('{}'.format(28 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        
        xl.flux_list = ['hx', 'hy']
        
        
        fig, axs = plt.subplots(2, 1)
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        hx_text = AnchoredText('{}'.format('$h_x$'), 
                                     loc='upper right')
        
        hxg_text = AnchoredText('{}'.format('$h_x/ \sqrt{g}$: [$m^{-2}$]'), 
                                     loc='upper center')
        
        hy_text = AnchoredText('{}'.format('$h_y$'), 
                                     loc='upper center')
        
        hyg_text = AnchoredText('{}'.format('$h_y/ \sqrt{g}$: [$m^{-2}$]'), 
                                     loc='upper center')
        
        text_list = [hx_text, hy_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                if ind == 1:
                    A_label_bol = False
                else:
                    A_label_bol = True
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    i_name = dat_name, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = A_label_bol)
                    
            
            axs[ind].legend(loc= 'upper left')
        
        axs[1].set_xlabel('poloidal angle')
        axs[0].set_title('geometry coefficients at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)




if topic == 'fluxes_no_psch':
    
    if xl.withshift == True and xl.withseries == False:
        
        xl.flux_iout_loader()
    
        psi_st = 19
        psi_ed = 38
        
        psi_dic = {'st': psi_st, 'ed': psi_ed}
        
        
        pol_list_a = []
        for i in range(36):
            pol_list_a.append('{}'.format(28 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        
        flux_list = ['pol_flux_no_psch', 'rad_flux_no_psch', 'bx', 'bz']
        
        
        fig, axs = plt.subplots(4, 1)
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        pol_text = AnchoredText('{}'.format('Poloidal flux with no PSch $\Gamma_x$: [$m^{-1} s^{-1}$]'), 
                                     loc='upper center')
        
        rad_text = AnchoredText('{}'.format('Radial flux with no PSch$\Gamma_y$: [$m^{-1} s^{-1}$]'), 
                                     loc='upper center')
        
        bx_text = AnchoredText('{}'.format('$b_x$: [T]'), 
                                     loc='upper center')
        
        bz_text = AnchoredText('{}'.format('$b_z$: [T]'), 
                                     loc='upper center')
        
        text_list = [pol_text, rad_text, bx_text, bz_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                # if ind % 2 == 0:
                #     text = pol_text
                    
                # else:
                #     text = rad_text
                
                if ind == 3:
                    A_label_bol = False
                else:
                    A_label_bol = True
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    i_name = dat_name, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = A_label_bol)
                    
            
            axs[ind].legend(loc= 'upper right')
        
        axs[3].set_xlabel('poloidal angle')
        axs[0].set_title('fluxes and magnetic strength at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)

"""





if topic == 'annual_review_mag':
    
    if xl.withshift == True and xl.withseries == False:
        
        xl.flux_iout_loader()
    
        # psi_st = 19
        # psi_ed = 38
        
        # psi_dic = {'st': psi_st, 'ed': psi_ed}
        
        
        pol_list_a = []
        for i in range(36):
            pol_list_a.append('{}'.format(28 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        
        flux_list = ['bx', 'bz']
        
        
        fig, axs = plt.subplots(3, 1)
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        pol_text = AnchoredText('{}'.format('Poloidal flux with no PSch $\Gamma_x$: [$m^{-1} s^{-1}$]'), 
                                     loc='upper center')
        
        bx_text = AnchoredText('{}'.format('(a) $B_{pol}$: [T]'), 
                                     loc='upper center')
        
        bz_text = AnchoredText('{}'.format('(b) $B_{tor}$: [T]'), 
                                     loc='lower center')
        
        xzratio_text = AnchoredText('{}'.format('(c) $B_{pol}/B_{tor}$'), 
                                     loc='lower center')
        
        text_list = [bx_text, bz_text, xzratio_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                if ind == 1:
                    A_label_bol = False
                else:
                    A_label_bol = True
                
                dat = xl.data['iout_data'][dat_name][aa]
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = dat, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, no_A_label = A_label_bol)
                
                axs[ind].legend(loc= 'lower left') 
                
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                
                sk = int(pol_list_a[0])
                sd = int(pol_list_a[-1]) + 1
                    
                bx_mag = xl.data['iout_data']['bx'][aa][psi_dic['st']:psi_dic['ed'], sk:sd]
                bz_mag = xl.data['iout_data']['bz'][aa][psi_dic['st']:psi_dic['ed'], sk:sd]
                
                axs[2].add_artist(text_list[2])
                
                mag_ratio = np.divide(bx_mag[0, :], bz_mag[0, :])
                
                axs[2].plot(ang_list, mag_ratio, linestyle= '-', 
                    color= color_dic[aa], label= 'A = {}'.format(A_dic[aa]))
                
                 
                
        
        axs[2].set_xlabel('poloidal angle')
        axs[0].set_title('magnetic strength and ratio at the separatrix')
        
        plt.subplots_adjust(hspace=.0)
        
        

if topic == 'mag_contour':
    
    if xl.withshift == True and xl.withseries == False:
        
        xl.flux_iout_loader()
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        
        for aa in xl.data['dircomp']['multi_shift']:
            
            data = xl.data['iout_data']['bz'][aa]
            xl.plot_change_data(data = data, log_bar = False, bounds = None,
                    itername = aa, quant = 'magnetic field z', ma100 = False,
                    color_dic = color_dic, A_dic = A_dic)

        

if topic == 'varify_vpara':
    
    if xl.withshift == True and xl.withseries == False:
        
        xl.flux_iout_loader()
    
        
        pol_list_a = []
        for i in range(36):
            pol_list_a.append('{}'.format(28 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        vpara = xl.data['iout_load_quant']['vpara']
        drvpara = xl.data['iout_load_quant']['derive vpara']
        
        flux_list = [vpara[0], drvpara[0]]
        
        
        fig, axs = plt.subplots(2, 1)
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        vpara_text = AnchoredText('{}'.format('$\sqrt{g}/h_x b_x v_{\parallel} n_i$: [$m s^{-1}$]'), 
                                     loc='upper center')
        
        ngvpara_text = AnchoredText('{}'.format('$b_x v_{\parallel} n_i$: [$m^{-2}s^{-1}$]'), 
                                     loc='upper center')
        
        text_list = [vpara_text, ngvpara_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                if ind == 1:
                    A_label_bol = False
                else:
                    A_label_bol = True
                
                dat = xl.data['iout_data'][dat_name][aa]
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = dat, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
                    
            # axs[ind].set_yscale("log")
            axs[ind].legend(loc= 'upper right')
        
        axs[1].axhline(y=0, color = 'black', linestyle = '--', label= '$b_x v_{\parallel} n_i$ = 0')
        axs[1].set_xlabel('poloidal angle')
        axs[0].set_title('v parallel flux at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)
        
        
        
        # vpara = xl.data['iout_load_quant']['vpara']
        # drvpara = xl.data['iout_load_quant']['derive vpara']
        
        mag = xl.data['iout_load_quant']['mag']
        
        fig, axs = plt.subplots()
        
        # color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
        #              'dot7': 'blue', 'one': 'purple'}
        # A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
        #           'dot7': '2.8', 'one': '3.4'}
        vpara_text = AnchoredText('{}'.format('v parallel flux ratio'), 
                                     loc='upper center')
        
        text_list = [vpara_text]
        
        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            
            A_label_bol = True
            
            if aa == 'org':
                pass
            else:
                
            
                vpara_derive = xl.data['iout_data'][drvpara[0]][aa]
                drvpara_org = xl.data['iout_data'][drvpara[0]]['org']
                
                drvpara_ratio = np.divide(vpara_derive, drvpara_org)
                
                dat = drvpara_ratio
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = dat, art_text = text_list[0], axs = axs, 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
                
                """
                
                bx = xl.data['iout_data'][mag[0]][aa]
                bx_org = xl.data['iout_data'][mag[0]]['org']
                
                
                bz = xl.data['iout_data'][mag[1]][aa]
                bz_org = xl.data['iout_data'][mag[1]]['org']
                
                bxbz = np.multiply(bx, bz)
                bxbz_org = np.multiply(bx_org, bz_org)
                
                bb = xl.data['iout_data'][mag[2]][aa]
                bb_org = xl.data['iout_data'][mag[2]]['org']
                
                
                squbb = np.multiply(bb, bb)
                bmag = np.divide(bxbz, squbb)
                
                squbb_org = np.multiply(bb_org, bb_org)
                bmag_org = np.divide(bxbz_org, bb_org)
                
                m_dat = np.divide(bmag, bmag_org)
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = m_dat, art_text = text_list[0], axs = axs, 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '--')
                
                """
                
                """
                
                sqg = xl.data['iout_data']['sqrt_g'][aa]
                sqg_org = xl.data['iout_data']['sqrt_g']['org']
                
                g_dat = np.divide(sqg, sqg_org)
                
                           
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = g_dat, art_text = text_list[0], axs = axs, 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '--')
                
                """
                
                hz = xl.data['iout_data']['hz'][aa]
                hz_org = xl.data['iout_data']['hz']['org']
                
                hz_dat = np.divide(hz, hz_org)
                
                           
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = hz_dat, art_text = text_list[0], axs = axs, 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '--')
                
                
        
        axs.legend(loc= 'upper right')
        
        axs.set_xlabel('poloidal angle')
        axs.set_title('v parallel flux ratio at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)
        
        
        gradn = xl.data['iout_load_quant']['gradn']
        drgradn = xl.data['iout_load_quant']['derive gradn']
        
        flux_list = [gradn[0], drgradn[0]]
        
        
        fig, axs = plt.subplots(2, 1)
        
        polgradn_text = AnchoredText('{}'.format('$\sqrt{g}/h_x^2 D_n \partial n_i/ \partial x$: [$s^{-1}$]'), 
                                     loc='upper center')
        
        ngpolgradn_text = AnchoredText('{}'.format('$1/h_x D_n \partial n_i/ \partial x$: [$m^{-2}s^{-1}$]'), 
                                     loc='upper center')
        
        radgradn_text = AnchoredText('{}'.format('$\sqrt{g}/h_y^2 D_n \partial n_i/ \partial x$: [$s^{-1}$]'), 
                                     loc='lower center')
        
        ngradgradn_text = AnchoredText('{}'.format('$1/h_y D_n \partial n_i/ \partial x$: [$m^{-2}s^{-1}$]'), 
                                     loc='lower center')
        
        text_list = [polgradn_text, ngpolgradn_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                if ind == 1:
                    A_label_bol = False
                else:
                    A_label_bol = True
                
                dat = xl.data['iout_data'][dat_name][aa]
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = dat, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
                    
            # axs[ind].set_yscale("log")
            axs[ind].legend(loc= 'upper left')
        
        axs[1].set_xlabel('poloidal angle')
        axs[0].set_title('gradn poloidal flux at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)
        
        flux_list = [gradn[1], drgradn[1]]
        
        
        fig, axs = plt.subplots(2, 1)
        radgradn_text = AnchoredText('{}'.format('$\sqrt{g}/h_y^2 D_n \partial n_i/ \partial x$: [$s^{-1}$]'), 
                                     loc='lower center')
        
        ngradgradn_text = AnchoredText('{}'.format('$1/h_y D_n \partial n_i/ \partial x$: [$m^{-2}s^{-1}$]'), 
                                     loc='lower center')
        
        text_list = [radgradn_text, ngradgradn_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                if ind == 1:
                    A_label_bol = False
                else:
                    A_label_bol = True
                    
                dat = xl.data['iout_data'][dat_name][aa]
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = dat, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
                    
            
            axs[ind].legend(loc= 'upper right')
        
        axs[1].set_xlabel('poloidal angle')
        axs[0].set_title('gradn radial flux at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)
        
        
        
        ni = xl.data['iout_load_quant']['ni']
        
        flux_list = ni
        
        
        fig, axs = plt.subplots(2, 1)
        
        ni_text = AnchoredText('{}'.format('ion density: [$m^{-3}$]'), 
                                     loc='upper center')
        
        ne_text = AnchoredText('{}'.format('electron density: [$m^{-3}$]'), 
                                     loc='upper center')
        
        text_list = [ni_text, ne_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                ni = xl.data['iout_data'][dat_name][aa]
                
                A_label_bol = False
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = ni, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
                
                    
            
        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            A_label_bol = True
            
            ne = xl.data['b2fstate'][aa]['ne']
            
            xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                input_dat = ne, art_text = text_list[1], axs = axs[1], 
                color_dic = color_dic, A_dic = A_dic, 
                no_A_label = A_label_bol, input_ls = '-')
            
            
          
        
        axs[1].legend(loc= 'upper left')
        
        axs[1].set_xlabel('poloidal angle')
        axs[0].set_title('ion and electron density at the separatrix')
        
        plt.subplots_adjust(hspace=.0)
        
        test = xl.data['b2fstate']['org']['ne']
        
        print(test.shape)


        
        vpp = xl.data['iout_load_quant']['vp and p']
        
        flux_list = vpp
        
        
        fig, axs = plt.subplots(2, 1)
        
        vp_text = AnchoredText('{}'.format('parallel velocity: [$m/s$]'), 
                                     loc='upper center')
        
        p_text = AnchoredText('{}'.format('electric potential: [V]'), 
                                     loc='upper center')
        
        text_list = [vp_text, p_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                vpp_dat = xl.data['iout_data'][dat_name][aa]
                
                if ind == 1:
                    A_label_bol = False
                else:
                    A_label_bol = True
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = vpp_dat, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
            
            
            axs[ind].legend(loc= 'upper right')
        
        
        
        axs[1].set_xlabel('poloidal angle')
        axs[0].set_title('parallel velocity and potential at the separatrix')
        
        plt.subplots_adjust(hspace=.0)
        
        
        fig, axs = plt.subplots()
        
        vpd_text = AnchoredText('{}'.format('parallel velocity difference: [$m/s$]'), 
                                     loc='upper center')
        
        text_list = [vpd_text]
        
            
        for aa in xl.data['dircomp']['multi_shift']:
            
            if aa == 'org':
                pass
            
            else:
                
            
                ang_list = xl.data['angle']['angle_list'][aa]
                
                vpp_dat = xl.data['iout_data'][vpp[0]][aa]
                vpp_dat_org = xl.data['iout_data'][vpp[0]]['org']
                
                sub_vpp_dat = vpp_dat_org - vpp_dat
                
                A_label_bol = True
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = sub_vpp_dat, art_text = text_list[0], axs = axs, 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
            
            
        axs.legend(loc= 'upper right')
        axs.set_xlabel('poloidal angle')
        axs.set_title('parallel velocity and potential at the separatrix')
        
        plt.subplots_adjust(hspace=.0)
        
        
        fig, axs = plt.subplots()
        
        R_text = AnchoredText('{}'.format('R: [m]'), 
                                     loc='upper center')
        
        text_list = [R_text]
        
            
        for aa in xl.data['dircomp']['multi_shift']:
            
                
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            R_dat = xl.data['grid']['RadLoc'][aa]
            
            
            A_label_bol = True
            
            xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                input_dat = R_dat, art_text = text_list[0], axs = axs, 
                color_dic = color_dic, A_dic = A_dic, 
                no_A_label = A_label_bol, input_ls = '-')
            
            
        axs.legend(loc= 'upper right')
        axs.set_xlabel('poloidal angle')
        axs.set_title('R coordinate at the separatrix')
        
        plt.subplots_adjust(hspace=.0)
        
        
        def Rfit(x, r, d, z, b): 
            return r*np.cos(x + d*np.sin(x) - z*np.sin(2*x) + b*np.sin(3*x))
        
        ang_list = xl.data['angle']['angle_list']['org']
        
        
        sk = int(pol_list_a[0])
        sd = int(pol_list_a[-1]) + 1
            
        psi_st = 17
        psi_ed = 38
        
        psi_dic = {'st': psi_st, 'ed': psi_ed}
        nnp = 3
        
        R_dat = xl.data['grid']['RadLoc']['org'][psi_dic['st']:psi_dic['ed'], sk:sd]
        avg_Rdat = []
        
        print('len pol_list is :{}, len psi ref is: {}'.format(len(pol_list_a), len(R_dat[0, :])))
        
        for xa in range(len(pol_list_a)):
            avg_Rdat.append(0.5*(float(R_dat[nnp + 1, xa]) + float(R_dat[nnp + 2, xa])))
        
        
        def Fit_R(ang, R_dat, rcenter):
            
            ang_rad = np.zeros([len(ang)])
            
            for k, at in enumerate(ang):
                ang_rad[k] = at* np.pi / 180
            
            # ar = np.asarray(ang_rad)
            print(type(ang_rad))
            # print(ang_rad)
                
            R_var = np.zeros([len(R_dat)])
            for k, avr in enumerate(R_dat):
                
                R_var[k] = avr - rcenter
                
            
            # Rvr = np.asarray(R_var)
            print(type(R_var))
            # print(R_var)
            
            
            p0 = [0.5, -0.01, 0.01, -0.01]
            popt_R, pcov_R = curve_fit(Rfit, ang_rad, R_var, p0)
            
            Rcoord_fit = rcenter*np.ones([len(R_dat)]) + Rfit(ang_rad, popt_R[0], 
                                    popt_R[1], popt_R[2], popt_R[3])
            
            
            
              
            fit_R_dic = {'Rcoord_fit': Rcoord_fit, 'popt_R': popt_R}
               
            return fit_R_dic
        
        
        rct = xl.data['gfile']['g']['rcentr']
        R_fitdic = Fit_R(ang = ang_list, R_dat = avg_Rdat, rcenter = rct)
        print(R_fitdic['popt_R'])
        
        fig, axs = plt.subplots()
        
        Rfit_text = AnchoredText('{}'.format('R and R fit: [m]'), 
                                     loc='upper center')
        
           
        text_list = [Rfit_text]
        shift_dic = {'org': 0, 'dot3': 0.3, 'dot5': 0.5, 'dot7': 0.7, 'one': 1}
        
        
        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            R_dat = xl.data['grid']['RadLoc'][aa]
            

            A_label_bol = True
            
            xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                input_dat = R_dat, art_text = text_list[0], axs = axs, 
                color_dic = color_dic, A_dic = A_dic, 
                no_A_label = A_label_bol, input_ls = '-')
            
            
            CoordR_fit = R_fitdic['Rcoord_fit'] + shift_dic[aa]
            
            axs.plot(ang_list, CoordR_fit, linestyle = '--', 
                color= color_dic[aa])
            
            
        
        
        axs.legend(loc= 'upper right')
    
        axs.set_xlabel('poloidal angle')
        axs.set_title('R and R fit')
        
        plt.subplots_adjust(hspace=.0)
        
        
        
        
        
        
        drtcoe = xl.data['iout_load_quant']['derive tcoe']
        sing_tcoe = xl.data['iout_load_quant']['single tcoe']
        
        
        flux_list = [sing_tcoe[0]]
        
        
        fig, axs = plt.subplots()
        
        drtcoex_text = AnchoredText('{}'.format('transport coefficient x: [$m^2/s^{-1}$]'), 
                                     loc='upper center')
        
        drtcoey_text = AnchoredText('{}'.format('transport coefficient y: [$m^2/s^{-1}$]'), 
                                     loc='upper center')
        
        singletcoe_text = AnchoredText('{}'.format('transport coefficient : [$m^2/s^{-1}$]'), 
                                     loc='upper center')
        
        
        
        text_list = [drtcoex_text, drtcoey_text, singletcoe_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                tcoe_dat = xl.data['iout_data'][dat_name][aa]
                

                A_label_bol = True
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = tcoe_dat, art_text = text_list[ind], axs = axs, 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
            
            
            axs.legend(loc= 'upper right')
        
        
        
        axs.set_xlabel('poloidal angle')
        axs.set_title('transport coefficient at the separatrix')
        
        plt.subplots_adjust(hspace=.0)
        
        
        corrdpc = xl.data['iout_load_quant']['corrdpc']
        drcorrdpc = xl.data['iout_load_quant']['derive corrdpc']
        
        
        flux_list = [corrdpc[0], drcorrdpc[0]]
        
        
        fig, axs = plt.subplots(2, 1)
        
        corrdpc_text = AnchoredText('{}'.format('corrdpc flux: [$1/s$]'), 
                                     loc='upper center')
        
        drcorrdpc_text = AnchoredText('{}'.format('derive corrdpc flux: [$m^{-2} s^{-1}$]'), 
                                     loc='upper center')
        
        text_list = [corrdpc_text, drcorrdpc_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                dpc_dat = xl.data['iout_data'][dat_name][aa]
                
                if ind == 1:
                    A_label_bol = False
                else:
                    A_label_bol = True
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = dpc_dat, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
            
            
            axs[ind].legend(loc= 'upper right')
        
        
        
        axs[1].set_xlabel('poloidal angle')
        axs[0].set_title('parallel velocity and potential at the separatrix')
        
        plt.subplots_adjust(hspace=.0)



if topic == 'geo_coe_test':
    
    if xl.withshift == True and xl.withseries == False:
        
        pol_list_st = 29
        pol_list_range = 32
        
        
        
        pol_list_a = []
        for i in range(pol_list_range):
            pol_list_a.append('{}'.format(pol_list_st + i))
            
        
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        pol_list_b = []
        for i in range(pol_list_range + 2):
            pol_list_b.append('{}'.format(pol_list_st -1 + i))
        
        
        fig, axs = plt.subplots()
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        
        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            R_dat = xl.data['grid']['RadLoc'][aa]
            Z_dat = xl.data['grid']['VertLoc'][aa]
            
            R_avg = xl.weight_generater(pol_list = pol_list_b, input_dat = R_dat, 
                                     itername = aa, nnp= 1)
            Z_avg = xl.weight_generater(pol_list = pol_list_b, input_dat = Z_dat, 
                                     itername = aa, nnp = 1)
            
            R_mid = []
            Z_mid = []
            
            for k in range(len(R_avg)-1):
                
                R_mid.append(0.5*(R_avg[k]+ R_avg[k + 1]))
                Z_mid.append(0.5*(Z_avg[k]+ Z_avg[k + 1]))
            
            # print(len(R_avg))
            # print(len(R_mid))
            
            R_diff = []
            Z_diff = []
            arc = []
            
            for k in range(len(R_mid)-1):
                
                length_squre = (R_mid[k] - R_mid[k + 1])**2 + (Z_mid[k] - Z_mid[k + 1])**2
                
                arc.append(np.sqrt(length_squre))
                
                R_length = (R_mid[k] - R_mid[k + 1])**2
                
                R_diff.append(np.sqrt(R_length))
                
                Z_length = (Z_mid[k] - Z_mid[k + 1])**2
                
                Z_diff.append(np.sqrt(Z_length))
            
            # print('the length of R_diff is {}'.format(len(R_diff)))
        
            
            axs.plot(ang_list, Z_diff, linestyle = '-', 
                color= color_dic[aa])

        
        
        
        fig, axs = plt.subplots()
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        
        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            R_dat = xl.data['grid']['RadLoc'][aa]
            Z_dat = xl.data['grid']['VertLoc'][aa]
            
            R_avg_a = xl.weight_generater(pol_list = pol_list_a, input_dat = R_dat, 
                                     itername = aa, nnp= 1)
            R_avg_b = xl.weight_generater(pol_list = pol_list_a, input_dat = R_dat, 
                                     itername = aa, nnp= 0)
            
            Z_avg_a = xl.weight_generater(pol_list = pol_list_a, input_dat = Z_dat, 
                                     itername = aa, nnp = 1)
            Z_avg_b = xl.weight_generater(pol_list = pol_list_a, input_dat = Z_dat, 
                                     itername = aa, nnp = 0)
            
            
            R_diff = []
            Z_diff = []
            arc = []
            
            for k in range(len(pol_list_a)):
                
                length_squre = (R_avg_a[k] - R_avg_b[k])**2 + (Z_avg_a[k] - Z_avg_b[k])**2
                
                arc.append(np.sqrt(length_squre))
                
                R_length = (R_avg_a[k] - R_avg_b[k])**2
                
                R_diff.append(np.sqrt(R_length))
                
                Z_length = (Z_avg_a[k] - Z_avg_b[k])**2
                
                Z_diff.append(np.sqrt(Z_length))
            
            # print('the length of R_diff is {}'.format(len(R_diff)))
        
            
            axs.plot(ang_list, arc, linestyle = '-', 
                color= color_dic[aa])


                    
if topic == 'geo_coe_fit':
    
    
    if xl.withshift == True and xl.withseries == False:
        
        xl.flux_iout_loader()
    
        
        pol_list_a = []
        for i in range(32):
            pol_list_a.append('{}'.format(29 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        
        
        xl.flux_list = ['hx', 'hy']
        
        
        fig, axs = plt.subplots(2, 1)
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        vpara_text = AnchoredText('{}'.format('$h_x$: [m]'), 
                                     loc='upper center')
        
        ngvpara_text = AnchoredText('{}'.format('$h_y$: [m]'), 
                                     loc='upper center')
        
        text_list = [vpara_text, ngvpara_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                if ind == 1:
                    A_label_bol = False
                else:
                    A_label_bol = True
                
                dat = xl.data['iout_data'][dat_name][aa]
                
                xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = dat, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, nnp = 1,
                    no_A_label = A_label_bol, input_ls = '-')
                    
            axs[ind].legend(loc= 'upper right')
        
        
        axs[1].set_xlabel('poloidal angle')
        axs[0].set_title('v parallel flux at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)
    

if topic == 'geo_coe_ymatch':
    
    
    if xl.withshift == True and xl.withseries == False:
        
        xl.flux_iout_loader()
    
        
        pol_list_a = []
        for i in range(30):
            pol_list_a.append('{}'.format(29 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        
        
        flux_list = ['hy']
        
        
        fig, axs = plt.subplots()
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        vpara_text = AnchoredText('{}'.format('$h_x$: [m]'), 
                                     loc='upper center')
        
        ngvpara_text = AnchoredText('{}'.format('$h_y$: [m]'), 
                                     loc='upper center')
        
        text_list = [vpara_text, ngvpara_text]
        
            
        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            
            A_label_bol = False

            
            dat = xl.data['iout_data']['hy'][aa]
            
            xl.flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                input_dat = dat, art_text = text_list[0], axs = axs, 
                color_dic = color_dic, A_dic = A_dic, nnp = 1,
                no_A_label = A_label_bol, input_ls = '-')
            
            
            
            R_dat = xl.data['grid']['RadLoc'][aa]
            Z_dat = xl.data['grid']['VertLoc'][aa]
            
            R_avg_a = xl.weight_generater(pol_list = pol_list_a, input_dat = R_dat, 
                                     itername = aa, nnp= 1)
            
            R_avg_b = xl.weight_generater(pol_list = pol_list_a, input_dat = R_dat, 
                                     itername = aa, nnp= 0)
            
            Z_avg_a = xl.weight_generater(pol_list = pol_list_a, input_dat = Z_dat, 
                                     itername = aa, nnp = 1)
            
            Z_avg_b = xl.weight_generater(pol_list = pol_list_a, input_dat = Z_dat, 
                                     itername = aa, nnp = 0)
            
            
            R_diff = []
            Z_diff = []
            arc = []
            
            for k in range(len(pol_list_a)):
                
                length_squre = (R_avg_a[k] - R_avg_b[k])**2 + (Z_avg_a[k] - Z_avg_b[k])**2
                
                arc.append(np.sqrt(length_squre))
            
            # print('the length of R_diff is {}'.format(len(R_diff)))
        
            
            axs.plot(ang_list, arc, linestyle = '--', 
                color= color_dic[aa])
            
            
                
        axs.legend(loc= 'upper right')
        
        
        axs.set_xlabel('poloidal angle')
        axs.set_title('v parallel flux at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)




"""

for item in flux_qu_list:
    
    res_qu_dic = {}
    for aa in xl.data['dircomp']['multi_shift']:
        if item.split('_')[1] == 'x':
            data, res_qu = xl.load_iout_multi(name1 = fcoe_list[0], name2 = item, 
                                input_name= 'poloidal_flux', itername = aa)
            print('flux check: {}'.format(fcoe_list[0]))
            print('flux check: {}'.format(item))
            res_qu_dic[aa] = data
        
        elif item.split('_')[1] == 'y':
            data, res_qu = xl.load_iout_multi(name1 = fcoe_list[1], name2 = item,
                                input_name= 'radial_flux', itername = aa)
            
            print('flux check: {}'.format(fcoe_list[1]))
            print('flux check: {}'.format(item))
            
            res_qu_dic[aa] = data
  
        else:
            print('a bug! flag is {}'.format(item.split('_')[1]))
    
    xl.data['iout_data'][res_qu] = res_qu_dic
    

    print(res_qu)
    res_qu_list.append(res_qu)

"""

"""

Retire the supporting function:
    

def flux_poloidal_sub(input_dat, itername, pol_list, ang_list, art_text,
                      axs, color_dic, A_dic, no_A_label, input_ls):
        
        sk = int(pol_list[0])
        sd = int(pol_list[-1]) + 1
            
        psi_st = 17
        psi_ed = 38
        
        psi_dic = {'st': psi_st, 'ed': psi_ed}
        np = 1
        
        psi_ref = xl.data['psi']['psival']['org'][psi_dic['st']:psi_dic['ed'], sk:sd]
        avg_psi = []
        
        print('len pol_list is :{}, len psi ref is: {}'.format(len(pol_list), len(psi_ref[0, :])))
        
        for xa in range(len(pol_list)):
            avg_psi.append(0.5*(float(psi_ref[np + 1, xa]) + float(psi_ref[np + 2, xa])))
        # avg_psi = 0.5*(psi_ref[:, np + 1] + psi_ref[:, np + 2])
        
        psi_w_dic = {}
        
        for aa in xl.data['dircomp']['multi_shift']:
            
            psi_dat = xl.data['psi']['psival'][aa][psi_dic['st']:psi_dic['ed'], sk:sd]
            
            if aa == 'org':
                
                # ln = int(len(pol_list))
                # print('look at type: {}'.format(type(ln)))
                w_psi = []
                for xa in range(len(pol_list)):
                    w_psi.append(0.5)
                
                
            else:
                
                # w_psi = np.zeros(len(pol_list))
                # print('avg shape: {}'.format(avg_psi.shape))
                w_psi = []
                for xa in range(len(pol_list)):
                    w_psi.append((avg_psi[xa] - psi_dat[np + 1, xa])/(psi_dat[np + 2, xa] - psi_dat[np + 1, xa]))
                
            psi_w_dic[aa] = w_psi 
        
        axs.add_artist(art_text)
        
        if input_dat.shape == (98, 38):
            
            plot_dat = input_dat[sk:sd, psi_dic['st']:psi_dic['ed']]
            
            wei = psi_w_dic[itername]
            avg_dat = []
            
            for wa in range(len(pol_list)):
                avg_dat.append((1- wei[wa])*plot_dat[wa, np + 1] + wei[wa]*plot_dat[wa, np + 2])
                
            
        
        elif input_dat.shape == (36, 96):
            
            wei = psi_w_dic[itername]
            avg_dat = []
            
            
            plot_dat = input_dat[psi_dic['st']:psi_dic['ed'], sk:sd]
            
            
            for wa in range(len(pol_list)):
                avg_dat.append((1- wei[wa])*plot_dat[np + 0, wa] + wei[wa]*plot_dat[np + 1, wa])
        
        
        elif input_dat.shape == (38, 98):
            
            wei = psi_w_dic[itername]
            avg_dat = []
            
            
            plot_dat = input_dat[psi_dic['st']:psi_dic['ed'], sk:sd]
            
            
            for wa in range(len(pol_list)):
                avg_dat.append((1- wei[wa])*plot_dat[np + 1, wa] + wei[wa]*plot_dat[np + 2, wa])
    
        
        
        if no_A_label:
            
            axs.plot(ang_list, avg_dat, linestyle = input_ls, 
                color= color_dic[itername])
        else:
            
            axs.plot(ang_list, avg_dat, linestyle = input_ls, 
                color= color_dic[itername], label= 'A = {}'.format(A_dic[itername]))
            

def general_iout_loader(tpl_list):
    
    qu_list = []
    
    for tpl_qu in tpl_list:
        
        quant = xl.load_iout(filename = tpl_qu[0], simple_quant = tpl_qu[1])
        print(quant)
        # print(flux_qu.split('_'))
        qu_list.append(quant)
    
    return qu_list


def derive_no_geo_quant(qu_list, gcoe_list, pol_name, rad_name):
    
    res_qu_list = []
    
    for item in qu_list:
        
        res_qu_dic = {}
        for aa in xl.data['dircomp']['multi_shift']:
            if item.split('_')[1] == 'x':
                data, res_qu = xl.load_iout_multi(name1 = gcoe_list[0], name2 = item, 
                                    input_name= pol_name, itername = aa)
                print('flux check: {}'.format(gcoe_list[0]))
                print('flux check: {}'.format(item))
                res_qu_dic[aa] = data
            
            elif item.split('_')[1] == 'y':
                data, res_qu = xl.load_iout_multi(name1 = gcoe_list[1], name2 = item,
                                    input_name= rad_name, itername = aa)
                
                print('flux check: {}'.format(gcoe_list[1]))
                print('flux check: {}'.format(item))
                
                res_qu_dic[aa] = data
      
            else:
                print('a bug! flag is {}'.format(item.split('_')[1]))
        
        xl.data['iout_data'][res_qu] = res_qu_dic
        
    
        print(res_qu)
        res_qu_list.append(res_qu)
        
        
    
    return res_qu_list


def flux_iout_loader():
    
    
    f_tuple5 = [('hx.dat', 'hx'),('vol.dat', 'sqrt_g')]
    f_tuple6 = [('hy.dat', 'hy'), ('vol.dat', 'sqrt_g')]
    
    flux_tuple = [('b2npc11_fnax001.dat', 'flux_x_0'), ('b2npc11_fnay001.dat', 'flux_y_0')]
        
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
    
    flux_qu_list = general_iout_loader(tpl_list = flux_tuple)
    
    res_qu_list = derive_no_geo_quant(qu_list = flux_qu_list, gcoe_list = fcoe_list, 
                                pol_name = 'poloidal_flux', rad_name = 'radial_flux')
    
    "=============== geometric data loader ==========="
    
    hz_tuple = [('hz.dat', 'hz')]
    
    hz_list = general_iout_loader(tpl_list = hz_tuple)
    
    
    
    "=============== magnetic data loader ==========="
    
    mag_tuple = [('bbx.dat', 'bx'), ('bbz.dat', 'bz'), ('bb.dat', 'B')]
    
    mag_list = general_iout_loader(tpl_list = mag_tuple)
    
    
    "=============== flux no psch data loader ==========="
    
    
    flux_no_psch_tuple = [('b2tfnb_fnbx001.dat', 'pol_flux_no_psch'), 
                          ('b2tfnb_fnby001.dat', 'rad_flux_no_psch')] 
    
    flux_no_psch_list = general_iout_loader(tpl_list = flux_no_psch_tuple)    
    
    "=============== psch flux data loader ==========="
    
    
    psch_flux_tuple = [('b2tfnb_fnbPSchx001.dat', 'psch_x'), 
                       ('b2tfnb_fnbPSchx001.dat','psch_y')]
    
    psch_list = general_iout_loader(tpl_list = psch_flux_tuple)
    
    derive_psch_list = derive_no_geo_quant(qu_list = psch_list, gcoe_list = fcoe_list, 
                                pol_name = 'poloidal_psch_flux', 
                                rad_name = 'radial_psch_flux')
    
    "=============== v_parallel flux data loader ==========="
    
    vpara_tuple = [('b2tfnb_bxuanax001.dat', 'vpara_x')]
    
    vpara_list = general_iout_loader(tpl_list = vpara_tuple)
    
    derive_vpara_list = derive_no_geo_quant(qu_list = vpara_list, gcoe_list = fcoe_list, 
                                pol_name = 'poloidal_vpara_flux', 
                                rad_name = 'radial_vpara_flux')
    
    
    "=============== grad den flux data loader ==========="
    
    
    gradn_tuple = [('b2tfnb_dPat_mdf_gradnax001.dat', 'gradn_x'), 
                      ('b2tfnb_dPat_mdf_gradnay001.dat', 'gradn_y')]
    
    gradn_list = general_iout_loader(tpl_list = gradn_tuple)
    
    derive_gradn_list = derive_no_geo_quant(qu_list = gradn_list, gcoe_list = fcoe_list, 
                                pol_name = 'poloidal_gradn_flux', 
                                rad_name = 'radial_gradn_flux')
    
    "=============== ion density data loader ==========="
    
    ni_tuple = [('b2npc11_na001.dat', 'ni')]
    
    ni_list = general_iout_loader(tpl_list = ni_tuple)
    
    
    "=============== parallel velocity data loader ==========="
    
    vp_and_p_tuple = [('b2npmo_ua001.dat', 'vp'), ('b2nppo_po.dat', 'e_potential')]
    
    vp_and_p_list = general_iout_loader(tpl_list = vp_and_p_tuple)
    
    
    "=============== transcoe data loader ==========="
    
    tcoe_tuple = [('b2trno_cdnax001.dat', 'tcoe_x'), ('b2trno_cdnay001.dat', 'tcoe_x')]
    
    tcoe_list = general_iout_loader(tpl_list = tcoe_tuple)
    
    derive_tcoe_list = derive_no_geo_quant(qu_list = tcoe_list, gcoe_list = fcoe_list, 
                                pol_name = 'poloidal_tcoe', 
                                rad_name = 'radial_tcoe')
    
    single_tcoe_tuple = [('b2tqna_dna0001.dat', 'single_tcoe')]
    
    single_tcoe_list = general_iout_loader(tpl_list = single_tcoe_tuple)
    
    
    "=============== v_corr_dpc data loader ==========="
    
    corrdpc_tuple = [('b2tfnb_dpccornax001.dat', 'corrdpc_x')]
    
    corrdpc_list = general_iout_loader(tpl_list = corrdpc_tuple)
    
    derive_corrdpc_list = derive_no_geo_quant(qu_list = corrdpc_list, gcoe_list = fcoe_list, 
                                pol_name = 'poloidal_corrdpc', 
                                rad_name = 'radial_corrdpc')
    
    
    qu_list_dic = {'geo_coe': fcoe_list, 'mag': mag_list, 'flux': flux_qu_list,
                   'psch': psch_list, 'derive psch': derive_psch_list, 
                   'vpara': vpara_list, 'derive vpara': derive_vpara_list,
                   'gradn': gradn_list, 'derive gradn': derive_gradn_list,
                   'ni': ni_list, 'hz': hz_list, 'vp and p': vp_and_p_list,
                   'tcoe': tcoe_list, 'derive tcoe': derive_tcoe_list,
                   'single tcoe': single_tcoe_list, 'corrdpc': corrdpc_list, 
                   'derive corrdpc': derive_corrdpc_list,
                   'derived flux': res_qu_list, 'flux no psch': flux_no_psch_list}
    
    
    xl.data['iout_load_quant'] = qu_list_dic


"""



