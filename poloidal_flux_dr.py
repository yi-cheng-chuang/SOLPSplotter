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
from matplotlib.offsetbox import AnchoredText
import scipy.stats as stats
from matplotlib.colors import LogNorm
from matplotlib import cm
from numpy import ma

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



xl = spc.PlotContour(DefaultSettings = d, loadDS = lex)

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
    
xl.opacity_data_fit(pol_list = poloidal_index_list)
xl.calc_pol_angle(pol_list = poloidal_index_list, plot_angle= False)

# xl.opacity_poloidal_plot(log_flag = False, save_pdf = False)

xl.neuden_percent()


['all_fluxes', 'fluxes_no_geo', 'fluxes_geo', 'neu_den', 
 'coe_check', 'mag_contour', 'fluxes_no_psch', 'mag_pol']

topic = 'mag_pol'


"""

('b2tfnb_bxuanax001.dat', False), ('hx.dat', True), ('hy.dat', True), 
('b2npc11_fnax001.dat', False), ('hz.dat', True)

"""


def flux_poloidal_sub(itername, pol_list, ang_list, i_name, art_text,
                      axs, color_dic, psi_dic, A_dic, no_A_label):
        
        sk = int(pol_list[0])
        sd = int(pol_list[-1]) + 1
            
        
        plot_dat = xl.data['iout_data'][i_name][itername][psi_dic['st']:psi_dic['ed'], sk:sd]
        
        axs.add_artist(art_text)
        
        if no_A_label:
            
            axs.plot(ang_list, plot_dat[0, :], linestyle='-', 
                color= color_dic[itername])
        else:
            
            axs.plot(ang_list, plot_dat[0, :], linestyle='-', 
                color= color_dic[itername], label= 'A = {}'.format(A_dic[itername]))
            
        
def flux_iout_loader():
    
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
    
    
    print(res_qu_list)
    
    
    
def flux_iout_loader_2():
    
    flux_tuple = [('b2tfnb_fnbx001.dat', 'pol_flux_no_psch'), ('b2tfnb_fnby001.dat', 'rad_flux_no_psch'),
                  ('bbx.dat', 'bx'), ('bbz.dat', 'bz'), ('bb.dat', 'B')]
    
    
    flux_qu_list = []
    
    for flux_tpl in flux_tuple:
        
        flux_qu = xl.load_iout(filename = flux_tpl[0], simple_quant = flux_tpl[1])
        print(flux_qu)
        # print(flux_qu.split('_'))
        flux_qu_list.append(flux_qu)
    
    res_qu_list = []
    
    

                    
if topic == 'all_fluxes':
    
    if xl.withshift == True and xl.withseries == False:
        
        flux_iout_loader()
    
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
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    i_name = dat_name, art_text = text_list[ind], axs = axs[ind], 
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
        
        flux_iout_loader()
    
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
            
            flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
            i_name = 'flux_x_0', art_text = pol_text, axs = axs[0], 
            color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = True)
            
        axs[0].axhline(y=0, color = 'black', linestyle = '--', label= '$\Gamma_x$ = 0')
        
        
        
        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    i_name = 'flux_y_0', art_text = rad_text, axs = axs[1], 
            color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = False)
            
            
        axs[1].set_xlabel('poloidal angle')
        axs[1].legend(loc= 'upper right')
        axs[0].legend(loc= 'upper right')
        axs[0].set_title('Particle flux at separatrix')
        
        
        plt.subplots_adjust(hspace=.0)


if topic == 'fluxes_no_geo':
    
    if xl.withshift == True and xl.withseries == False:
        
        flux_iout_loader()
    
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
            
            flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
            i_name = 'poloidal_flux', art_text = pol_text, axs = axs[0], 
            color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = True)
            
        axs[0].axhline(y=0, color = 'black', linestyle = '--', label= '$\Gamma_x$ = 0')
        
        
        
        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    i_name = 'radial_flux', art_text = rad_text, axs = axs[1], 
            color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = False)
            
            
        axs[1].set_xlabel('poloidal angle')
        axs[1].legend(loc= 'upper right')
        axs[0].legend(loc= 'upper right')
        axs[0].set_title('Particle flux at separatrix')
        
        
        plt.subplots_adjust(hspace=.0)




if topic == 'neu_den':
    
    if xl.withshift == True and xl.withseries == False:

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
            
            axs[1].plot(ang_list, neuden_dat[-1, :], '-', color= color_dic[aa])


if topic == 'coe_check':
    
    if xl.withshift == True and xl.withseries == False:
        
        flux_iout_loader()
    
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
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    i_name = dat_name, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = A_label_bol)
                    
            
            axs[ind].legend(loc= 'upper right')
        
        axs[3].set_xlabel('poloidal angle')
        axs[0].set_title('coefficients at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)


if topic == 'fluxes_no_psch':
    
    if xl.withshift == True and xl.withseries == False:
        
        flux_iout_loader_2()
    
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
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    i_name = dat_name, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = A_label_bol)
                    
            
            axs[ind].legend(loc= 'upper right')
        
        axs[3].set_xlabel('poloidal angle')
        axs[0].set_title('fluxes and magnetic strength at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)


if topic == 'mag_pol':
    
    if xl.withshift == True and xl.withseries == False:
        
        flux_iout_loader_2()
    
        psi_st = 19
        psi_ed = 38
        
        psi_dic = {'st': psi_st, 'ed': psi_ed}
        
        
        pol_list_a = []
        for i in range(36):
            pol_list_a.append('{}'.format(28 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        
        flux_list = ['B', 'bx', 'bz']
        
        
        fig, axs = plt.subplots(3, 1)
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        pol_text = AnchoredText('{}'.format('Poloidal flux with no PSch $\Gamma_x$: [$m^{-1} s^{-1}$]'), 
                                     loc='upper center')
        
        B_text = AnchoredText('{}'.format('B: [T]'), 
                                     loc='upper center')
        
        bx_text = AnchoredText('{}'.format('$B_x$: [T]'), 
                                     loc='upper center')
        
        bz_text = AnchoredText('{}'.format('$B_z$: [T]'), 
                                     loc='lower center')
        
        text_list = [B_text, bx_text, bz_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                # if ind % 2 == 0:
                #     text = pol_text
                    
                # else:
                #     text = rad_text
                
                if ind == 2:
                    A_label_bol = False
                else:
                    A_label_bol = True
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    i_name = dat_name, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, psi_dic = psi_dic, A_dic = A_dic, no_A_label = A_label_bol)
                    
            
            axs[ind].legend(loc= 'upper right')
        
        axs[2].set_xlabel('poloidal angle')
        axs[0].set_title('magnetic strength at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)




if topic == 'mag_contour':
    
    if xl.withshift == True and xl.withseries == False:
        
        flux_iout_loader_2()
        
        color_dic = {'org': 'red', 'dot3': 'darkorange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        
        for aa in xl.data['dircomp']['multi_shift']:
            
            data = xl.data['iout_data']['bz'][aa]
            xl.plot_change_data(data = data, log_bar = False, bounds = None,
                    itername = aa, quant = 'magnetic field z', ma100 = False,
                    color_dic = color_dic, A_dic = A_dic)
                
        
