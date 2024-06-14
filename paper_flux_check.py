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
from scipy.optimize import curve_fit



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


['geo_coe_fit', 'varify_vpara']

topic = 'varify_vpara'


"""

('b2tfnb_bxuanax001.dat', False), ('hx.dat', True), ('hy.dat', True), 
('b2npc11_fnax001.dat', False), ('hz.dat', True)

"""
psi_st = 17
psi_ed = 38

psi_dic = {'st': psi_st, 'ed': psi_ed}


def weight_generater(pol_list):
    
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
    
    return np, psi_st, psi_ed, sk, sd, psi_w_dic






def flux_poloidal_sub(input_dat, itername, pol_list, ang_list, art_text,
                      axs, color_dic, A_dic, no_A_label, input_ls):
        """
        
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
        
        """
        
        np, psi_st, psi_ed, sk, sd, psi_w_dic = weight_generater(pol_list)
        
        
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
    
    "=============== source data loader ==========="
    
    s_tuple = [('b2npc11_sna001.dat', 'source')]
    
    s_list = general_iout_loader(tpl_list = s_tuple)
    
    
    
    
    
    qu_list_dic = {'geo_coe': fcoe_list, 'mag': mag_list, 'flux': flux_qu_list,
                   'psch': psch_list, 'derive psch': derive_psch_list, 
                   'vpara': vpara_list, 'derive vpara': derive_vpara_list,
                   'gradn': gradn_list, 'derive gradn': derive_gradn_list,
                   'ni': ni_list, 'hz': hz_list, 'vp and p': vp_and_p_list,
                   'tcoe': tcoe_list, 'derive tcoe': derive_tcoe_list,
                   'single tcoe': single_tcoe_list, 'corrdpc': corrdpc_list, 
                   'derive corrdpc': derive_corrdpc_list, 'source': s_list,
                   'derived flux': res_qu_list, 'flux no psch': flux_no_psch_list}
    
    
    xl.data['iout_load_quant'] = qu_list_dic



if topic == 'geo_coe_test':
    
    if xl.withshift == True and xl.withseries == False:
        
        
        pol_list_a = []
        for i in range(36):
            pol_list_a.append('{}'.format(26 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        
        
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





                    
if topic == 'geo_coe_fit':
    
    
    if xl.withshift == True and xl.withseries == False:
        
        flux_iout_loader()
    
        
        pol_list_a = []
        for i in range(36):
            pol_list_a.append('{}'.format(26 + i))
              
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        
        
        flux_list = ['hx', 'hy']
        
        
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
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = dat, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
                    
            axs[ind].legend(loc= 'upper right')
        
        
        axs[1].set_xlabel('poloidal angle')
        axs[0].set_title('v parallel flux at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)
    










if topic == 'varify_vpara':
    
    if xl.withshift == True and xl.withseries == False:
        
        flux_iout_loader()
    
        
        pol_list_a = []
        for i in range(36):
            pol_list_a.append('{}'.format(26 + i))
              
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
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
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
        vpara_text = AnchoredText('{}'.format('$b_x = B_{pol}/B$'), 
                                     loc='upper center')
        
        text_list = [vpara_text]
        
        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            
            A_label_bol = False
                
                
            bx = xl.data['iout_data'][mag[0]][aa]
            
            bb = xl.data['iout_data'][mag[2]][aa]
            
            m_dat = np.divide(bx, bb)
            
                       
            flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                input_dat = m_dat, art_text = text_list[0], axs = axs, 
                color_dic = color_dic, A_dic = A_dic, 
                no_A_label = A_label_bol, input_ls = '-')
                
                
        
        axs.legend(loc= 'upper right')
        
        axs.set_xlabel('poloidal angle')
        axs.set_title('$b_x$ at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)
        
        
        gradn = xl.data['iout_load_quant']['gradn']
        drgradn = xl.data['iout_load_quant']['derive gradn']
        
        flux_list = [gradn[0], drgradn[0]]
        
        
        fig, axs = plt.subplots(2, 1)
        
        polgradn_text = AnchoredText('{}'.format('$\sqrt{g}/h_x^2 D_n \partial n_i/ \partial x$: [$s^{-1}$]'), 
                                     loc='upper center')
        
        ngpolgradn_text = AnchoredText('{}'.format('$1/h_x D_n \partial n_i/ \partial x$: [$m^{-2}s^{-1}$]'), 
                                     loc='upper center')
        
        radgradn_text = AnchoredText('{}'.format('$\sqrt{g}/h_y^2 D_n \partial n_i/ \partial y$: [$s^{-1}$]'), 
                                     loc='lower center')
        
        ngradgradn_text = AnchoredText('{}'.format('$1/h_y D_n \partial n_i/ \partial y$: [$m^{-2}s^{-1}$]'), 
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
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
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
        radgradn_text = AnchoredText('{}'.format('$\sqrt{g}/h_y^2 D_n \partial n_i/ \partial y$: [$s^{-1}$]'), 
                                     loc='lower center')
        
        ngradgradn_text = AnchoredText('{}'.format('$1/h_y D_n \partial n_i/ \partial y$: [$m^{-2}s^{-1}$]'), 
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
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = dat, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
                    
            
            axs[ind].legend(loc= 'upper right')
        
        axs[1].set_xlabel('poloidal angle')
        axs[0].set_title('gradn radial flux at the separatrix')
        
        
        plt.subplots_adjust(hspace=.0)
        
        
        
        
        source = xl.data['iout_load_quant']['source']
        
        flux_list = source
        
        
        fig, axs = plt.subplots()
        
        source_text = AnchoredText('{}'.format('source: [$m^{-3} s^{-1}$]'), 
                                     loc='upper center')
        
        text_list = [source_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                

                A_label_bol = False

                    
                sn_dat = xl.data['iout_data'][dat_name][aa]
                g_dat = xl.data['iout_data']['sqrt_g'][aa]
                
                source_dat = np.divide(sn_dat, g_dat)
                
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = source_dat, art_text = text_list[0], axs = axs, 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
                    
            
            axs.legend(loc= 'upper right')
        
        axs.set_xlabel('poloidal angle')
        axs.set_title('source at the separatrix')
        
        
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
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = ni, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
                
                    
            
        for aa in xl.data['dircomp']['multi_shift']:
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            A_label_bol = True
            
            ne = xl.data['b2fstate'][aa]['ne']
            
            flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
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
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = vpp_dat, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
            
            
            axs[ind].legend(loc= 'upper right')
        
        
        
        axs[1].set_xlabel('poloidal angle')
        axs[0].set_title('parallel velocity and potential at the separatrix')
        
        plt.subplots_adjust(hspace=.0)
        
        
        fig, axs = plt.subplots(2, 1)
        
        vpp = xl.data['iout_load_quant']['vp and p']
        
        flux_list = [vpp[0]]
        
        
        vp_text = AnchoredText('{}'.format('parallel velocity: [$m/s$]'), 
                                     loc='upper center')
        
        text_list = [vp_text]
        
        for ind, dat_name in enumerate(flux_list):
            
            for aa in xl.data['dircomp']['multi_shift']:
                
                ang_list = xl.data['angle']['angle_list'][aa]
                
                vpp_dat = xl.data['iout_data'][dat_name][aa]
                

                A_label_bol = False

                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = vpp_dat, art_text = text_list[0], axs = axs[0], 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
            
            
            axs[0].legend(loc= 'upper right')
        
        
        axs[0].axhline(y=0, color = 'black', linestyle = '--', label= '$b_x v_{\parallel} n_i$ = 0')
        # axs.set_xlabel('poloidal angle')
        axs[0].set_title('parallel velocity and difference at the separatrix')
        
        # plt.subplots_adjust(hspace=.0)
        
        
        
        # fig, axs = plt.subplots(2, 1)
        
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
                
                sub_vpp_dat = vpp_dat - vpp_dat_org
                
                A_label_bol = True
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = sub_vpp_dat, art_text = text_list[0], axs = axs[1], 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
            
            
        axs[1].legend(loc= 'upper right')
        axs[1].set_xlabel('poloidal angle')
        # axs.set_title('parallel velocity difference $v_{\parallel} - v^{MAST}_{\parallel}$ at the separatrix')
        
        plt.subplots_adjust(hspace=.0)
        
        
        fig, axs = plt.subplots()
        
        R_text = AnchoredText('{}'.format('R: [m]'), 
                                     loc='upper center')
        
        text_list = [R_text]
        
            
        for aa in xl.data['dircomp']['multi_shift']:
            
                
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            R_dat = xl.data['grid']['RadLoc'][aa]
            
            
            A_label_bol = True
            
            flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
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
        nnp = 1
        
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
            

            A_label_bol = False
            
            flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
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
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
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
                
                flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
                    input_dat = dpc_dat, art_text = text_list[ind], axs = axs[ind], 
                    color_dic = color_dic, A_dic = A_dic, 
                    no_A_label = A_label_bol, input_ls = '-')
            
            
            axs[ind].legend(loc= 'upper right')
        
        
        
        axs[1].set_xlabel('poloidal angle')
        axs[0].set_title('parallel velocity and potential at the separatrix')
        
        plt.subplots_adjust(hspace=.0)


"""

'v parallel velocity try and error'



if aa == 'org':
    pass
else:
    
    vpara_derive = xl.data['iout_data'][drvpara[0]][aa]
    drvpara_org = xl.data['iout_data'][drvpara[0]]['org']
    
    drvpara_ratio = np.divide(vpara_derive, drvpara_org)
    
    dat = drvpara_ratio
    
    flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
        input_dat = dat, art_text = text_list[0], axs = axs, 
        color_dic = color_dic, A_dic = A_dic, 
        no_A_label = A_label_bol, input_ls = '-')
    

    
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
    
    
    hz = xl.data['iout_data']['hz'][aa]
    hz_org = xl.data['iout_data']['hz']['org']
    
    hz_dat = np.divide(hz, hz_org)
    

    
    sqg = xl.data['iout_data']['sqrt_g'][aa]
    sqg_org = xl.data['iout_data']['sqrt_g']['org']
    
    g_dat = np.divide(sqg, sqg_org)
    
               
    flux_poloidal_sub(itername = aa, pol_list = pol_list_a, ang_list = ang_list,
        input_dat = g_dat, art_text = text_list[0], axs = axs, 
        color_dic = color_dic, A_dic = A_dic, 
        no_A_label = A_label_bol, input_ls = '--')


"""




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



