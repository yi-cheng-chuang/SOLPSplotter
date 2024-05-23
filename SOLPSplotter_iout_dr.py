# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:51:44 2024

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

topic = 'Q3'


"""

('b2tfnb_bxuanax001.dat', False), ('hx.dat', True), ('hy.dat', True), 
('b2npc11_fnax001.dat', False), ('hz.dat', True)

"""

if topic == 'Q1-1':
    f_tuple1 = [('b2npco_sna000.dat', 'g_source_0'), ('vol.dat', 'sqrt_g')]
    f_tuple2 = [('b2npc11_sna001.dat', 'g_source_1'), ('vol.dat', 'sqrt_g')]
    f_tuple3 = [('b2npco_dnadt000.dat', 'g_dnadt_0'), ('vol.dat', 'sqrt_g')]
    f_tuple4 = [('b2npc11_dnadt001.dat', 'g_dnadt_1'), ('vol.dat', 'sqrt_g')]
    
    f_list1 = [f_tuple1, f_tuple2, f_tuple3, f_tuple4]
    single = False
    f_spc_list = [f_tuple3]
    
    
    if single:
        for tuples in f_spc_list:
            qu = xl.load_iout_ratio(file_tuple = tuples)
            print(qu)
            
            xl.iout_contour_plot(quant = qu, log_bar = False)
            xl.iout_contour_plot(quant = qu, log_bar = True)
            
    else:
        for tuples in f_list1:
            qu = xl.load_iout_ratio(file_tuple = tuples)
            print(qu)
            
            xl.iout_contour_plot(quant = qu, log_bar = False)
            xl.iout_contour_plot(quant = qu, log_bar = True)
    
    
    # f_tuple1 = [('b2npco_sna000.dat', 'g_source_0'), ('b2npco_dnadt000.dat', 'g_dnadt_0')]
    f_tuple2 = [('b2npc11_sna001.dat', 'g_source_1'), ('b2npc11_dnadt001.dat', 'g_dnadt_1')]
    
    f_list1 = [f_tuple2]

   
    for tuples in f_list1:
        qu = xl.load_iout_ratio(file_tuple = tuples)
        print(qu)
        
        xl.iout_contour_plot(quant = qu, log_bar = False)
        xl.iout_contour_plot(quant = qu, log_bar = True)

elif topic == 'Q1-2':
    
    f_tuple5 = [('hx.dat', 'hx'),('vol.dat', 'sqrt_g')]
    f_tuple6 = [('hy.dat', 'hy'), ('vol.dat', 'sqrt_g')]
    
    flux_tuple = [('b2npco_fnax000.dat', 'flux_x_0'), ('b2npc11_fnax001.dat', 'flux_x_1'),
                ('b2npco_fnay000.dat', 'flux_y_0'), ('b2npc11_fnay001.dat', 'flux_y_1')]
    
    
    f_list2 = [f_tuple5, f_tuple6]
    fcoe_list = []
    
    
    for tuples in f_list2:
        qu = xl.load_iout_ratio(file_tuple = tuples)
        print(qu)
        fcoe_list.append(qu)
    
    # print(fcoe_list)
    
    flux_qu_list = []
    
    for flux_tpl in flux_tuple:
        
        flux_qu = xl.load_iout(filename = flux_tpl[0], simple_quant = flux_tpl[1])
        print(flux_qu)
        # print(flux_qu.split('_'))
        flux_qu_list.append(flux_qu)
    
    plot_flux_ratio = False
    
    if plot_flux_ratio:
        print(flux_qu_list)
        
        for i in range(2):
            qu = xl.load_iout_name_ratio(name1 = flux_qu_list[i], name2 = flux_qu_list[i + 2])
            print(qu)
            
            xl.iout_contour_plot(quant = qu, log_bar = False)
            xl.iout_contour_plot(quant = qu, log_bar = True)
    else:
        pass
        
    
    
    
    
    res_qu_list = []
    
    for item in flux_qu_list:
        if item.split('_')[1] == 'x':
            res_qu = xl.load_iout_multi(name1 = fcoe_list[0], name2 = item)
            print(res_qu)
            res_qu_list.append(res_qu)

            xl.iout_contour_plot(quant = res_qu, log_bar = False)
            xl.iout_contour_plot(quant = res_qu, log_bar = True)
        
        elif item.split('_')[1] == 'y':
            res_qu = xl.load_iout_multi(name1 = fcoe_list[1], name2 = item)
            print(res_qu)
            res_qu_list.append(res_qu)

            xl.iout_contour_plot(quant = res_qu, log_bar = False)
            xl.iout_contour_plot(quant = res_qu, log_bar = True)
        
        else:
            print('a bug! flag is {}'.format(item.split('_')[1]))

    print(res_qu_list)
    
    for i in range(2):
        qu = xl.load_iout_name_ratio(name1 = res_qu_list[i], name2 = res_qu_list[i + 2])
        print(qu)
        
        xl.iout_contour_plot(quant = qu, log_bar = False)
        xl.iout_contour_plot(quant = qu, log_bar = True)
        
    
elif topic == 'Q2':
    
    if xl.withshift == True and xl.withseries == False:
        
        
        gcoe_tuple = [('hx.dat', 'hx'), ('hy.dat', 'hy'), ('hz.dat', 'hz')]
        
        for item in gcoe_tuple:
            
            qu = xl.load_iout(filename = item[0], simple_quant = item[1])
            print(qu)
        
        gcoe_list = ['hx', 'hy', 'hz']
        for coe in gcoe_list:
            for aa in xl.data['dircomp']['multi_shift']:
                if aa == 'org':
                    pass
                else:
                    dif_qu = xl.load_iout_differ(name = coe, itername1 = aa, itername2 = 'org')
                    print(dif_qu)
                    
                    data, per_qu = xl.load_iout_name_ratio(name1 = dif_qu, name2 = coe,
                                stdname = 'org', itername = aa)
                    print(per_qu)
                    
                    xl.plot_change_data(data = data, log_bar = False, 
                             itername = aa, quant = per_qu, ma100 = True)
                    
elif topic == 'Q3':
    
    if xl.withshift == True and xl.withseries == False:
        
        
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
        
        # for rqu in res_qu_list:
            
        #     xl.iout_contour_plot(quant = rqu, log_bar= True, ma100= False, bounds = {})
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        
        
        
        for resqu in res_qu_list:
            
            resqu_per_dic = {}
            std_cp = 'org'
            for aa in xl.data['dircomp']['multi_shift']:
                if aa == std_cp:
                    pass
                else:
                    dif_qu = xl.load_differ(name = resqu, itername1 = aa, 
                                            itername2 = std_cp, setname = 'iout_data')
                    print(dif_qu)
                    
                    data, per_qu = xl.load_iout_name_ratio(name1 = dif_qu, name2 = resqu,
                                    stdname = std_cp, itername = aa, setname = 'iout_data')
                    print(per_qu)
                    resqu_per_dic[aa] = data
                    bon = {'max': 100, 'min': -100}
                    xl.plot_change_data(data = data, log_bar = False, bounds = bon,
                            itername = aa, quant = per_qu, ma100 = True,
                            color_dic = color_dic, A_dic = A_dic)
                
            xl.data['iout_data'][per_qu] = resqu_per_dic
        
        
        
        # plt.figure(figsize=(7,7))

        # for aa in xl.data['dircomp']['multi_shift']:
        #     color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
        #                  'dot7': 'blue', 'one': 'purple'}
        #     A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
        #               'dot7': '2.8', 'one': '3.4'}
        #     st = int(poloidal_index_list[0]) -1
        #     ed = int(poloidal_index_list[-1])
        #     sol_pol_flux = xl.data['iout_data']['poloidal_flux'][aa][20:, st:ed]
        #     mean_pol_flux = np.mean(sol_pol_flux, axis=0)
        #     std_pol_flux = np.std(sol_pol_flux, axis=0)
            
        #     norm_pol_flux = mean_pol_flux/ max(mean_pol_flux)
        #     neuden_d = xl.data['opacity_poloidal'][aa]['neutral_density']
        #     y = np.transpose(neuden_d)/ max(neuden_d)
        #     plt.errorbar(mean_pol_flux, y, xerr= std_pol_flux, fmt = '-', 
        #                  color = color_dic[aa], label= 'aspect ratio = {}'.format(A_dic[aa]))
        #     plt.title('neutral density and ion poloidal flux ')
        #     plt.xlabel('ion poloidal flux')
        #     plt.legend()


        # plt.figure(figsize=(7,7))

        # for aa in xl.data['dircomp']['multi_shift']:
        #     color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
        #                  'dot7': 'blue', 'one': 'purple'}
        #     A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
        #               'dot7': '2.8', 'one': '3.4'}
        #     st = int(poloidal_index_list[0]) -1
        #     ed = int(poloidal_index_list[-1])
        #     sol_pol_flux = xl.data['iout_data']['radial_flux'][aa][20:, st:ed]
        #     mean_pol_flux = np.mean(sol_pol_flux, axis=0)
        #     std_pol_flux = np.std(sol_pol_flux, axis=0)
        #     norm_pol_flux = mean_pol_flux/ max(mean_pol_flux)
            
            
        #     neuden_d = xl.data['opacity_poloidal'][aa]['neutral_density']
        #     y = np.transpose(neuden_d)/ max(neuden_d)
        #     plt.errorbar(mean_pol_flux, y, xerr= std_pol_flux, fmt = '-', 
        #                 color = color_dic[aa], label= 'aspect ratio = {}'.format(A_dic[aa]))
        #     plt.title('neutral density and ion radial flux correlation')
        #     plt.xlabel('ion radial flux')
        #     plt.legend()
        
        
        
    
        
        
        
        anglemean_pol_flux = []
        cv_pol_flux = []
        std_pol_flux = []
        neuden_list = []
        cv_neuden = []
        std_neuden = []
        anglemean_rad_flux = []
        std_rad_flux = []
        st = int(poloidal_index_list[0]) -1
        ed = int(poloidal_index_list[-1])
        one_angle = False
        ind = ed
        for aa in xl.data['dircomp']['multi_shift']:
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            sol_pol_flux = xl.data['iout_data']['poloidal_flux'][aa][20:, st:ed]
            mean_pol_flux = np.mean(sol_pol_flux, axis=0)
            anglemean_pol_flux.append(-1*np.mean(mean_pol_flux, axis=0))
            
            cv_pol_flux.append(-1*np.std(mean_pol_flux, axis=0) / np.mean(mean_pol_flux, axis=0))
            std_pol_flux.append(np.std(mean_pol_flux, axis=0))
            
            neuden_d = xl.data['opacity_poloidal'][aa]['neutral_density']
            neuden_list.append(np.mean(neuden_d, axis=0))
            # mean_neuden = np.mean(neuden_d, axis=0)
            cv_neuden.append(np.std(neuden_d, axis=0) / np.mean(neuden_d, axis=0))
            std_neuden.append(np.std(neuden_d, axis=0))
            
            sol_rad_flux = xl.data['iout_data']['radial_flux'][aa][20:, st:ed]
            mean_rad_flux = np.mean(sol_rad_flux, axis=0)
            anglemean_rad_flux.append(np.mean(mean_rad_flux, axis=0))
            std_rad_flux.append(np.std(mean_rad_flux, axis=0))
            
            
        fig, axs = plt.subplots(2, 1, figsize = (14, 7))
        aspect_list = [1.4, 2.0, 2.4, 2.8]
        color = 'tab:red'
        
        print(cv_pol_flux)
        print(cv_neuden)
        text1 = 'neutral density, ion poloidal flux verses aspect ratio'
        text2 = 'neutral density, ion radial flux verses aspect ratio'
        anchored_text1 = AnchoredText('(a){}'.format(text1), loc= 'upper left')
        anchored_text2 = AnchoredText('(b){}'.format(text2), loc= 'lower center')
        anchored_text3 = AnchoredText('neutral density', loc=4)
        
        axs[0].errorbar(aspect_list, anglemean_pol_flux, yerr= std_pol_flux, 
                      fmt = '-', color= color, label= 'poloidal flux')
        # ax1.scatter(aspect_list, anglemean_pol_flux, color= color)
        # ax1.set_ylabel('poloidal flux', color= color)
        axs[0].tick_params(axis='y', labelcolor = color)
        
        ax2 = axs[0].twinx()
        
        color2 = 'tab:blue'
        
        ax2.errorbar(aspect_list, neuden_list, yerr= std_neuden, 
                    fmt = '-', color= color2, label= 'neutral density')
        # ax2.scatter(aspect_list, neuden_list, color= color2)
        ax2.tick_params(axis='y', labelcolor = color2)
        
        # ax2.set_ylabel('neutral density', color= color2)
        
        axs[0].add_artist(anchored_text1)
        axs[0].legend(loc= 'upper right')
        ax2.legend(loc= 'lower center')
        # axs[0].set_xlabel('aspect ratio')
        
        axs[1].errorbar(aspect_list, anglemean_rad_flux, yerr= std_rad_flux, 
                      fmt = '-', color= color, label= 'radial flux')
        # ax1.scatter(aspect_list, anglemean_rad_flux, color= 'red')
        # ax1.set_ylabel('radial flux', color= color)
        axs[1].tick_params(axis='y', labelcolor = color)
        
        ax3 = axs[1].twinx()
        
        ax3.errorbar(aspect_list, neuden_list, yerr= std_neuden, 
                    fmt = '-', color= color2, label= 'neutral density')
        # ax2.set_ylabel('neutral density', color= color2)
        ax3.tick_params(axis='y', labelcolor = color2)
        
        
        axs[1].add_artist(anchored_text2)
        axs[1].set_xlabel('aspect ratio')
        axs[1].legend(loc = 'upper center')
        
        plt.subplots_adjust(hspace=.0)
        
        
        # fig_dir  = sps.set_figdir()
        # plt.savefig('{}/{}.pdf'.format(fig_dir, 'flux&neuden'), format='pdf')
        
        
        # st = int(poloidal_index_list[0]) -1
        # ed = int(poloidal_index_list[-1])
        
        
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                          'dot7': 'blue', 'one': 'purple'}
        
        anchored_text = AnchoredText('(a){}'.format('HFS poloidal flux [$m^{-1} s^{-1}$]'), 
                                     loc= 'lower center')
        
        rad_text = AnchoredText('(b){}'.format('HFS radial flux [$m^{-1} s^{-1}$]'), 
                                     loc= 'lower right')
        
        LFS_text = AnchoredText('(c){}'.format('LFS poloidal flux [$m^{-1} s^{-1}$]'), 
                                     loc= 'upper right')
        
        pol_list_a = []
        for i in range(9):
            pol_list_a.append('{}'.format(28 + i))

        
        psi_st = 19
        psi_ed = 36
        
        
        xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        fig, axs = plt.subplots(3, 1)
            
        sk = int(pol_list_a[0])
        sd = int(pol_list_a[-1])
            
            
        for ab in xl.data['dircomp']['multi_shift']:
            
                
            pol_flux_dat = xl.data['iout_data']['poloidal_flux'][ab][psi_st:psi_ed, sk:sd]
            
            # nor_pol_flux = (pol_flux_dat - pol_flux_dat.max())/(pol_flux_dat.max() - pol_flux_dat.min())
            
            nor_pol_flux = stats.zscore(pol_flux_dat)
            
            neuden_dat = np.transpose(xl.data['ft44'][ab]['dab2'][sk:sd, psi_st:psi_ed, 0])
            
            # nor_neuden = (neuden_dat - neuden_dat.min()) /(neuden_dat.max() - neuden_dat.min())
            
            nor_neuden = stats.zscore(neuden_dat)
            
            # slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(neuden_dat.flatten(), pol_flux_dat.flatten())
            
            slope, intercept, r_value, p_value, std_err = stats.linregress(neuden_dat.flatten(), pol_flux_dat.flatten())
            
            print_regress = False
            
            if print_regress:
                
                print('case {} slope = {}'.format(ab, slope))
                print('case {} intercept = {}'.format(ab, intercept))
                print('case {} std_err = {}'.format(ab, std_err))
                print('case {} r_value = {}'.format(ab, r_value))
                print('case {} p_value = {}'.format(ab, p_value))
            
            else:
                pass
            
            
            x = np.linspace(neuden_dat.min(), neuden_dat.max(), 50)
            y = slope*x + intercept
            
            
            axs[0].plot(x, y, '-', color = color_dic[ab], 
                            label= 'A= {}'.format(A_dic[ab]))
            axs[0].scatter(neuden_dat.flatten(), pol_flux_dat.flatten(), color= color_dic[ab])
        
        
        for ab in xl.data['dircomp']['multi_shift']:
            
                
            rad_flux_dat = xl.data['iout_data']['radial_flux'][ab][psi_st:psi_ed, sk:sd]
            
            nor_rad_flux = rad_flux_dat/(rad_flux_dat.max())
            
            
            neuden_dat = np.transpose(xl.data['ft44'][ab]['dab2'][sk:sd, psi_st:psi_ed, 0])
            
            nor_neuden = neuden_dat/(neuden_dat.max())
            
            slope, intercept, r_value, p_value, std_err = stats.linregress(neuden_dat.flatten(), rad_flux_dat.flatten())
            
            print_regress = False
            
            if print_regress:
                
                print('case {} slope = {}'.format(ab, slope))
                print('case {} intercept = {}'.format(ab, intercept))
                print('case {} std_err = {}'.format(ab, std_err))
                print('case {} r_value = {}'.format(ab, r_value))
                print('case {} p_value = {}'.format(ab, p_value))
            
            else:
                pass
            
            
            x = np.linspace(neuden_dat.min(), neuden_dat.max(), 50)
            y = slope*x + intercept
            
            
            axs[1].plot(x, y, '-', color = color_dic[ab], 
                            label= 'A= {}'.format(A_dic[ab]))
            axs[1].scatter(neuden_dat.flatten(), rad_flux_dat.flatten(), color= color_dic[ab])
        
        pol_list_b = []
        for i in range(6):
            pol_list_b.append('{}'.format(58 + i))
        
        sk = int(pol_list_b[0])
        sd = int(pol_list_b[-1])
        
        for ab in xl.data['dircomp']['multi_shift']:
            
                
            pol_flux_dat = xl.data['iout_data']['poloidal_flux'][ab][psi_st:psi_ed, sk:sd]
            nor_pol_flux = pol_flux_dat/(abs(pol_flux_dat).max())
            
            
            neuden_dat = np.transpose(xl.data['ft44'][ab]['dab2'][sk:sd, psi_st:psi_ed, 0])
            
            nor_neuden = neuden_dat/(neuden_dat.max())
            
            slope, intercept, r_value, p_value, std_err = stats.linregress(neuden_dat.flatten(), pol_flux_dat.flatten())
            
            print_regress = True
            
            if print_regress:
                
                print('case {} slope = {}'.format(ab, slope))
                print('case {} intercept = {}'.format(ab, intercept))
                print('case {} std_err = {}'.format(ab, std_err))
                print('case {} r_value = {}'.format(ab, r_value))
                print('case {} p_value = {}'.format(ab, p_value))
            
            else:
                pass
            
            # print(slope)
            # print(r_value)
            
            # print(intercept)
            
            x = np.linspace(neuden_dat.min(), neuden_dat.max(), 50)
            y = slope*x + intercept
            
            
            axs[2].plot(x, y, '-', color = color_dic[ab], 
                            label= 'A= {}'.format(A_dic[ab]))
            axs[2].scatter(neuden_dat.flatten(), pol_flux_dat.flatten(), color= color_dic[ab])
            
            
        
        
        axs[0].add_artist(anchored_text)
        axs[1].add_artist(rad_text)
        axs[2].add_artist(LFS_text)
        axs[0].legend(loc= 'lower right', fontsize=10)
        
        axs[2].set_xlabel('Neutral density [$m^{-3}$]')
        axs[0].set_title('Particle flux and neutral density correlation plot')
        plt.subplots_adjust(hspace=.0)
        
        fig.savefig('correlation.pdf')
        
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                          'dot7': 'blue', 'one': 'purple'}
        
        anchored_text = AnchoredText('{}'.format('Normalized poloidal flux [$m^{-1} s^{-1}$]'), loc= 'upper left')
        
        
        pol_list_a = []
        for i in range(9):
            pol_list_a.append('{}'.format(28 + i))
        
        pol_list_b = []
        for i in range(6):
            pol_list_b.append('{}'.format(58 + i))
        
        check_target = False
        
        if check_target:
            
            pol_list_c = []
            for i in range(10):
                pol_list_c.append('{}'.format(1 + i))
            
            pol_list_d = []
            for i in range(9):
                pol_list_d.append('{}'.format(87 + i))
        else:
            pass
            

        
        psi_st = 19
        psi_ed = 38
        
        
        # xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
        
        iter_list = [pol_list_a, pol_list_b]
        
        fig_n = len(iter_list)
        
        fig, axs = plt.subplots(fig_n, 1)
        
        for i, pl in enumerate(iter_list):
            
            sk = int(pl[0])
            sd = int(pl[-1])
            
            
        
        
            for ab in xl.data['dircomp']['multi_shift']:
                
                    
                pol_flux_dat = xl.data['iout_data']['poloidal_flux'][ab][psi_st:psi_ed, sk:sd]
                nor_pol_flux = pol_flux_dat/(abs(pol_flux_dat).max())
                
                
                # print(pol_flux_dat.max())
                # print(pol_flux_dat)
                
                
                
                neuden_dat = np.transpose(xl.data['ft44'][ab]['dab2'][sk:sd, psi_st:psi_ed, 0])
                
                nor_neuden = neuden_dat/(neuden_dat.max())
                
                rad_flux_dat = xl.data['iout_data']['radial_flux'][ab][psi_st:psi_ed, sk:sd]
                
                nor_rad_flux = rad_flux_dat/(rad_flux_dat.max())
                
                cor_list = stats.pearsonr(nor_neuden.flatten(), nor_pol_flux.flatten())
                
                # print(cor_list)
                
                slope, intercept, r_value, p_value, std_err = stats.linregress(nor_neuden.flatten(), nor_pol_flux.flatten())
                
                print_regress = False
                
                if print_regress:
                    
                    print(slope)
                    print(intercept)
                    print(std_err)
                    print(r_value)
                    print(p_value)
                
                else:
                    pass
                    
                
                x = np.linspace(nor_neuden.min(), nor_neuden.max(), 50)
                y = slope*x + intercept
                
                
                axs[i].plot(x, y, '-', color = color_dic[ab], 
                                label= 'A= {}'.format(A_dic[ab]))
                axs[i].scatter(nor_neuden.flatten(), nor_pol_flux.flatten(), color= color_dic[ab])
            
            
            axs[0].add_artist(anchored_text)
            
            
            axs[fig_n -1].set_xlabel('Normalized neutral density [$m^{-3}$]')
            
            plt.subplots_adjust(hspace=.0)
        
        
        
        
        
        
        
        
        
        
        geo_list = []
        geo_tuple = ('vol.dat', 'sqrt_g')
        geo_qu = xl.load_iout(filename = geo_tuple[0], simple_quant = geo_tuple[1])
        print(geo_qu)
        # print(flux_qu.split('_'))
        geo_list.append(geo_qu)
        
        
        
        s_list = []
        source_tuple = [('b2npc11_sna001.dat', 'source'), ('b2stel_sna_ion001.dat', 'source_ion'),
                        ('b2stel_sna_rec001.dat', 'source_rec'), ('b2stcx_sna_001.dat', 'source_cx'),
                        ('b2stbr_sna_eir001.dat', 'source_eir')]
        for s_tpl in source_tuple:
            
            s_qu = xl.load_iout(filename = s_tpl[0], simple_quant = s_tpl[1])
            print(s_qu)
            # print(flux_qu.split('_'))
            s_list.append(s_qu)
        
        
        
        
        res_key = []
        
        # for s_name in s_list:
        
        for i in range(len(s_list)):
            res_qu_dic = {}
            
            for aa in xl.data['dircomp']['multi_shift']:
            
                data, res_qu = xl.load_iout_divide(name1 = s_list[i], name2 = geo_list[0], 
                                    input_name= '{}_no_jacobian'.format(s_list[i]), itername = aa)
                
                res_qu_dic[aa] = data
                
                # xl.plot_change_data(data = data, log_bar = True, itername = aa, 
                #                 quant = res_qu, ma100 = False, bounds = {})
            
            xl.data['iout_data'][res_qu] = res_qu_dic
            res_key.append(res_qu)
        
        
        
        fig_n = len(iter_list)
        
        fig, axs = plt.subplots(fig_n, 1)
        
        for i, pl in enumerate(iter_list):
            
            sk = int(pl[0])
            sd = int(pl[-1])
        
        
            for at in xl.data['dircomp']['multi_shift']:
                
                    
                pol_flux_dat = xl.data['iout_data']['poloidal_flux'][at][psi_st:psi_ed, sk:sd]
                nor_pol_flux = pol_flux_dat/(abs(pol_flux_dat).max())
                
                
                # print(pol_flux_dat.max())
                # print(pol_flux_dat)
                
                
                s_term = xl.data['iout_data']['source_no_jacobian'][at][psi_st:psi_ed, sk:sd]
                
                nor_source = s_term/(s_term.max())
                
                # print(cor_list)
                
                slope, intercept, r_value, p_value, std_err = stats.linregress(nor_source.flatten(), nor_pol_flux.flatten())
                
                # print(slope)
                # print(intercept)
                # print(std_err)
                # print(r_value)
                # print(p_value)
                
                x = np.linspace(nor_source.min(), nor_source.max(), 50)
                y = slope*x + intercept
                
                
                # axs[i].plot(x, y, '-', color = color_dic[at], 
                #                 label= 'A= {}'.format(A_dic[at]))
                axs[i].scatter(nor_source.flatten(), nor_pol_flux.flatten(), color= color_dic[at])
            
            
            # axs[0].add_artist(anchored_text)
            
            
            axs[fig_n -1].set_xlabel('Normalized source rate [$m^{-3}$]')
            
            plt.subplots_adjust(hspace=.0)
        
        
        fig_n = len(iter_list)
        
        fig, axs = plt.subplots(fig_n, 1)
        
        for i, pl in enumerate(iter_list):
            
            sk = int(pl[0])
            sd = int(pl[-1])
        
        
            for at in xl.data['dircomp']['multi_shift']:
                
                
                neuden_dat = np.transpose(xl.data['ft44'][at]['dab2'][sk:sd, psi_st:psi_ed, 0])
                
                nor_neuden = neuden_dat/(neuden_dat.max())
                
                # print(pol_flux_dat.max())
                # print(pol_flux_dat)
                
                
                s_term = xl.data['iout_data']['source_no_jacobian'][at][psi_st:psi_ed, sk:sd]
                
                nor_source = s_term/(s_term.max())
                
                # print(cor_list)
                
                slope, intercept, r_value, p_value, std_err = stats.linregress(neuden_dat.flatten(), s_term.flatten())
                
                print(slope)
                print(intercept)
                print(std_err)
                print(r_value)
                print(p_value)
                
                x = np.linspace(neuden_dat.min(), neuden_dat.max(), 50)
                y = slope*x + intercept
                
                
                axs[i].plot(x, y, '-', color = color_dic[at], 
                                label= 'A= {}'.format(A_dic[at]))
                axs[i].scatter(neuden_dat.flatten(), s_term.flatten(), color= color_dic[at])
            
            
            # axs[0].add_artist(anchored_text)
            
            
            axs[fig_n -1].set_xlabel('Normalized neutral density [$m^{-3}$]')
            
            plt.subplots_adjust(hspace=.0)
        
        
        
        
        iter_list = [pol_list_a, pol_list_b]
        
        fig_n = len(iter_list)
        
        fig, axs = plt.subplots(fig_n, 1)
        
        for i, pl in enumerate(iter_list):
            
            sk = int(pl[0])
            sd = int(pl[-1])
            
            
        
        
            for ab in xl.data['dircomp']['multi_shift']:
                
                    
                pol_flux_dat = xl.data['iout_data']['poloidal_flux'][ab][psi_st:psi_ed, sk:sd]
                nor_pol_flux = pol_flux_dat/(abs(pol_flux_dat).max())
                
                
                # print(pol_flux_dat.max())
                # print(pol_flux_dat)
                
                
                
                neuden_dat = np.transpose(xl.data['ft44'][ab]['dab2'][sk:sd, psi_st:psi_ed, 0])
                
                nor_neuden = neuden_dat/(neuden_dat.max())
                
                rad_flux_dat = xl.data['iout_data']['radial_flux'][ab][psi_st:psi_ed, sk:sd]
                
                nor_rad_flux = rad_flux_dat/(rad_flux_dat.max())
                
                cor_list = stats.pearsonr(nor_neuden.flatten(), nor_rad_flux.flatten())
                
                # print(cor_list)
                
                slope, intercept, r_value, p_value, std_err = stats.linregress(nor_neuden.flatten(), nor_rad_flux.flatten())
                
                x = np.linspace(nor_neuden.min(), nor_neuden.max(), 50)
                y = slope*x + intercept
                
                
                axs[i].plot(x, y, '-', color = color_dic[ab], 
                                label= 'A= {}'.format(A_dic[ab]))
                axs[i].scatter(nor_neuden.flatten(), nor_rad_flux.flatten(), color= color_dic[ab])
            
            axs[fig_n -1].set_xlabel('Normalized neutral density [$m^{-3}$]')
            
            plt.subplots_adjust(hspace=.0)
        
        
        
        
        
        
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
        

        for ac in xl.data['dircomp']['multi_shift']:
            
            sk = int(pol_list_a[0])
            sd = int(pol_list_a[-1]) + 1
            
            # neuden_dat = np.transpose(xl.data['ft44'][ac]['dab2'][sk:sd, psi_st, 0])
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            pol_flux_dat_a = []
            pol_flux_dat_b = []
            
            for kk in pol_list_a:
                
                pol_flux_dat = xl.data['iout_data']['poloidal_flux'][ac][psi_st:psi_ed, int(kk)]
                
                pol_flux_dat_a.append(pol_flux_dat.max())
                pol_flux_dat_b.append(pol_flux_dat.min())
                
            
            pol_flux_dat = xl.data['iout_data']['poloidal_flux'][ac][psi_st:psi_ed, sk:sd]
            
            # pol_flux_dat_a = xl.data['iout_data']['poloidal_flux'][ac][psi_st, sk:sd]
            # pol_flux_dat_b = xl.data['iout_data']['poloidal_flux'][ac][psi_ed -1, sk:sd]
            
            axs[0].add_artist(pol_text)
            
            # axs[0].fill_between(ang_list, pol_flux_dat_a, pol_flux_dat_b, 
            #         color = color_dic[ac], alpha = 0.4, label= 'A = {}'.format(A_dic[ac]))
            
            axs[0].plot(ang_list, pol_flux_dat[0, :], linestyle='-', color= color_dic[ac])
            
            # axs[0].plot(ang_list, pol_flux_dat[-1, :], '-', color= color_dic[ac])
            
        axs[0].axhline(y=0, color = 'black', linestyle = '--', label= '$\Gamma_x$ = 0')
        
        
        
        
        # for aa in xl.data['dircomp']['multi_shift']:
            
        #     neuden_data_a = []
        #     neuden_data_b = []
            
        #     for kt in pol_list_a:
                
        #         neuden_data = xl.data['ft44'][aa]['dab2'][int(kt), psi_st:psi_ed]
                
        #         neuden_data_a.append(neuden_data.max())
        #         neuden_data_b.append(neuden_data.min())
            
            
        #     ang_list = xl.data['angle']['angle_list'][aa]
        
        
        #     neuden_dat = np.transpose(xl.data['ft44'][aa]['dab2'][sk:sd, psi_st:psi_ed, 0])
            
        #     axs[1].add_artist(neu_text)
            
        #     # axs[1].fill_between(ang_list, neuden_data_a, neuden_data_b, 
        #     #                  color= color_dic[aa], alpha = 0.4)
            
        #     axs[1].plot(ang_list, neuden_dat[0, :], linestyle='-', color= color_dic[aa])
            
        #     axs[1].plot(ang_list, neuden_dat[-1, :], '-', color= color_dic[aa])
            
            
            
        
        for aa in xl.data['dircomp']['multi_shift']:
            
            rad_data_a = []
            rad_data_b = []
            
            for kt in pol_list_a:
                
                rad_data = xl.data['iout_data']['flux_y_0'][aa][psi_st:psi_ed, int(kt)]
                
                rad_data_a.append(rad_data.max())
                rad_data_b.append(rad_data.min())
            
            rad_dat = xl.data['iout_data']['radial_flux'][aa][psi_st:psi_ed, sk:sd]
            
            ang_list = xl.data['angle']['angle_list'][aa]
            
            # axs[2].fill_between(ang_list, rad_data_a, rad_data_b, 
            #                   color= color_dic[aa], alpha = 0.4)
            
            axs[1].add_artist(rad_text)
            
            axs[1].plot(ang_list, rad_dat[0, :], linestyle='-', color= color_dic[aa], 
                        label= 'A = {}'.format(A_dic[aa]))
            
            # axs[2].plot(ang_list, rad_dat[-1, :], '-', color= color_dic[aa])
        
        
        axs[1].set_xlabel('poloidal angle')
        axs[1].legend(loc= 'upper right')
        axs[0].legend(loc= 'upper right')
        axs[0].set_title('Particle flux and neutral density at separatrix')
        
        
        plt.subplots_adjust(hspace=.0)
        
            
        
        
                    
elif topic == 'Q3-1':
    
    
    coe = 'dab2'
    
    A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
              'dot7': '2.8', 'one': '3.4'}
    color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                 'dot7': 'blue', 'one': 'purple'}
    
    for aa in xl.data['dircomp']['multi_shift']:
        if aa == 'org':
            pass
        else:
            dif_qu = xl.load_differ(name = coe, itername1 = aa, 
                                itername2 = 'org', setname = 'ft44')
            print(dif_qu)
            
            data, per_qu = xl.load_iout_name_ratio(name1 = dif_qu, name2 = coe,
                        stdname = 'org', itername = aa, setname = 'ft44')
            print(per_qu)
            
            neuden_data = np.transpose(data[:, :, 0])
            print(np.shape(neuden_data))
            
            bon = {'max': 100, 'min': -100}
            
            xl.plot_change_data(data = neuden_data, log_bar = True, bounds = bon, 
        itername = aa, quant = 'neutral density increase percentage', ma100 = True,
                color_dic = color_dic, A_dic = A_dic)
          
        
        
    fig, axs = plt.subplots(1, 2, sharey= True)  
    
    data = xl.data['ft44']['org'][coe]
    neuden_data = np.transpose(data[:, :, 0])
    
    datamap = np.abs(data)
    
    log_bar = True
    plot_2dval = neuden_data
    
    if log_bar:
        if np.all(plot_2dval == 0):
            print('data_file is an zero matrix')
            
        elif np.any(plot_2dval == 0):
            plot_2dval = ma.masked_where(plot_2dval <= 0, plot_2dval)
            
            datamap = np.abs(plot_2dval)
        else:
            
            datamap = np.abs(plot_2dval)
    
    
    
    CPB = cm.viridis
    Lnorm = LogNorm(vmin = datamap.min(), vmax = datamap.max())
    k = np.linspace(np.log10(datamap.min()), np.log10(datamap.max()), 11)
    
    lev = np.power(10, k)
    print(lev)
    
    levels = lev
    
    RadLoc = xl.data['grid']['RadLoc']['org']
    VertLoc = xl.data['grid']['VertLoc']['org']
    
    R_coord = RadLoc[1:37, 1:97]
    Z_coord = VertLoc[1:37, 1:97]
    
    axs[0].contourf(R_coord, Z_coord, datamap, levels= levels, 
                 cmap = CPB, norm = Lnorm)
    
    
    
    vessel = xl.data['vessel']['org']
    axs[0].plot(vessel[:,0]/1000, vessel[:,1]/1000, color = color_dic['org'])
    
    data = xl.data['ft44']['dot7'][coe]
    neuden_data = np.transpose(data[:, :, 0])
    
    log_bar = True
    plot_2dval = neuden_data
    
    if log_bar:
        if np.all(plot_2dval == 0):
            print('data_file is an zero matrix')
            
        elif np.any(plot_2dval == 0):
            plot_2dval = ma.masked_where(plot_2dval <= 0, plot_2dval)
            
            datamap = np.abs(plot_2dval)
        else:
            
            datamap = np.abs(plot_2dval)
    
    RadLoc = xl.data['grid']['RadLoc']['dot7']
    VertLoc = xl.data['grid']['VertLoc']['dot7']
    
    R_coord = RadLoc[1:37, 1:97]
    Z_coord = VertLoc[1:37, 1:97]
    
    axs[1].contourf(R_coord, Z_coord, datamap, levels= levels, 
                 cmap = CPB, norm = Lnorm)
    
    vessel = xl.data['vessel']['dot7']
    axs[1].plot(vessel[:,0]/1000, vessel[:,1]/1000, color = color_dic['dot7'])
    
    axs[0].set_ylabel('Z [m]')
    fig.supxlabel('R [m]')
    
    fig.suptitle('Neutral density [$m^{-3}$] contour plot')
       
    smap = cm.ScalarMappable(norm = Lnorm, cmap = CPB)    
    plt.colorbar(smap, ax= axs)
    
    # orientation='horizontal'
    
    
    
    pol_list_a = []
    for i in range(9):
        pol_list_a.append('{}'.format(28 + i))
    
    
    xl.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
    
    xl.data['angle']['angle_list'][aa]
    
    psi_st = 18
    psi_ed = 35
    sk = int(pol_list_a[0])
    sd = int(pol_list_a[-1]) + 1
    
    fig, axs = plt.subplots()
    
    
    for aa in xl.data['dircomp']['multi_shift']:
        
        data = xl.data['ft44'][aa][coe]
        neuden_data_a = np.transpose(data[sk:sd, psi_st, 0])
        neuden_data_b = np.transpose(data[sk:sd, psi_ed, 0])
        
        ang_list = xl.data['angle']['angle_list'][aa]
        
        axs.fill_between(ang_list, neuden_data_a, neuden_data_b, 
                         color= color_dic[aa], alpha = 0.4)
    
    
    
        
elif topic == 'Q3-2':
    
    geo_list = []
    geo_tuple = ('vol.dat', 'sqrt_g')
    geo_qu = xl.load_iout(filename = geo_tuple[0], simple_quant = geo_tuple[1])
    print(geo_qu)
    # print(flux_qu.split('_'))
    geo_list.append(geo_qu)
    
    
    
    s_list = []
    source_tuple = [('b2npc11_sna001.dat', 'source'), ('b2stel_sna_ion001.dat', 'source_ion'),
                    ('b2stel_sna_rec001.dat', 'source_rec'), ('b2stcx_sna_001.dat', 'source_cx'),
                    ('b2stbr_sna_eir001.dat', 'source_eir')]
    for s_tpl in source_tuple:
        
        s_qu = xl.load_iout(filename = s_tpl[0], simple_quant = s_tpl[1])
        print(s_qu)
        # print(flux_qu.split('_'))
        s_list.append(s_qu)
    
    
    
    
    res_key = []
    
    # for s_name in s_list:
    
    for i in range(len(s_list)):
        res_qu_dic = {}
        
        for aa in xl.data['dircomp']['multi_shift']:
        
            data, res_qu = xl.load_iout_divide(name1 = s_list[i], name2 = geo_list[0], 
                                input_name= '{}_no_jacobian'.format(s_list[i]), itername = aa)
            
            res_qu_dic[aa] = data
            
            # xl.plot_change_data(data = data, log_bar = True, itername = aa, 
            #                 quant = res_qu, ma100 = False, bounds = {})
        
        xl.data['iout_data'][res_qu] = res_qu_dic
        res_key.append(res_qu)
    
    
    
    coe = res_key[0]
    
    for aa in xl.data['dircomp']['multi_shift']:
        if aa == 'org':
            pass
        else:
            dif_qu = xl.load_differ(name = coe, itername1 = aa, 
                                itername2 = 'org', setname = 'iout_data')
            print(dif_qu)
            
            data, per_qu = xl.load_iout_name_ratio(name1 = dif_qu, name2 = coe,
                        stdname = 'org', itername = aa, setname = 'iout_data')
            print(per_qu)
            
            bon = {'max': 300, 'min': -300}
            
            xl.plot_change_data(data = data, log_bar = False, bounds = bon, 
                     itername = aa, quant = 'source rate increase percentage', ma100 = True)
    
    
    for res_name in res_key:
        
        plt.figure(figsize=(7,7))

        for aa in xl.data['dircomp']['multi_shift']:
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            st = int(poloidal_index_list[0]) -1
            ed = int(poloidal_index_list[-1])
            s_term = xl.data['iout_data'][res_name][aa][10:, st:ed]
            mean_s_term = np.mean(s_term, axis=1)
            std_s_term = np.std(s_term, axis=1)
            # norm_pol_flux = mean_pol_flux/ max(mean_pol_flux)
            
            psi = xl.data['psi']['psival'][aa][10:-2, st-1: ed-1]
            mean_psi = np.mean(psi, axis=1)
            std_psi = np.std(psi, axis=1)
            # neuden_d = xl.data['opacity_poloidal'][aa]['neutral_density']
            # y = np.transpose(neuden_d)/ max(neuden_d)
            plt.errorbar(mean_psi, mean_s_term, xerr = std_psi, fmt ='-', 
            color = color_dic[aa], label= 'aspect ratio = {}'.format(A_dic[aa]))
            plt.title(res_name)
            plt.xlabel('psiN')
            plt.legend()
    
    
    
    
    fig, ax = plt.subplots(figsize=(7,7))
    
    anchored_text = AnchoredText('source rate $m^{-3} s^{-1}$', loc=2)
    
    
    for aa in xl.data['dircomp']['multi_shift']:
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue', 'one': 'purple'}
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8', 'one': '3.4'}
        st = int(poloidal_index_list[0]) -1
        ed = int(poloidal_index_list[-1])
        s_term = xl.data['iout_data']['source_no_jacobian'][aa][10:, st:ed]
        mean_s_term = np.mean(s_term, axis=1)
        std_s_term = np.std(s_term, axis=1)
        # norm_pol_flux = mean_pol_flux/ max(mean_pol_flux)
        
        psi = xl.data['psi']['psival'][aa][10:-2, st-1: ed-1]
        mean_psi = np.mean(psi, axis=1)
        std_psi = np.std(psi, axis=1)
        # neuden_d = xl.data['opacity_poloidal'][aa]['neutral_density']
        # y = np.transpose(neuden_d)/ max(neuden_d)
        ax.errorbar(mean_psi, mean_s_term, xerr = std_psi, fmt ='-', 
        color = color_dic[aa])
        ax.set_xlabel('$\psi_N$')
    
    ax.add_artist(anchored_text) 

    
elif topic == 'Q5': 
    
    f_tuple1 = [('b2npco_sna000.dat', 'g_source_0'), ('vol.dat', 'sqrt_g')]
    f_tuple2 = [('b2npc11_sna001.dat', 'g_source_1'), ('vol.dat', 'sqrt_g')]
    f_tuple3 = [('b2npco_dnadt000.dat', 'g_dnadt_0'), ('vol.dat', 'sqrt_g')]
    f_tuple4 = [('b2npc11_dnadt001.dat', 'g_dnadt_1'), ('vol.dat', 'sqrt_g')]
    
    f_list1 = [f_tuple1, f_tuple2, f_tuple3, f_tuple4]
    single = False
    f_spc_list = [f_tuple3]
    
    
    if single:
        for tuples in f_spc_list:
            qu = xl.load_iout_ratio(file_tuple = tuples)
            print(qu)
            
            xl.iout_contour_plot(quant = qu, log_bar = False)
            xl.iout_contour_plot(quant = qu, log_bar = True)
            
    else:
        for tuples in f_list1:
            qu = xl.load_iout_ratio(file_tuple = tuples)
            print(qu)
            
            xl.iout_contour_plot(quant = qu, log_bar = False)
            xl.iout_contour_plot(quant = qu, log_bar = True)
    
    
    

elif topic == 'test':
    
    iout_quant = [ ('b2npco_sna000.dat', False), ('b2npc11_sna001.dat', False), 
                  ('b2stel_sna_ion000.dat', False), ('b2stbr_sna_eir001.dat', False), 
                  ('b2stcx_sna_001.dat', False), ('b2npc11_dnadt001.dat', False), ('b2stel_sna_rec000.dat', False)]



    for item in iout_quant:
        
        qu = xl.load_iout(filename = item[0], simple_quant = item[1])
        print(qu)

        xl.iout_contour_plot(quant = qu)
        

else:
    print('no topic!')





# for aa in xl.data['dircomp']['multi_shift']:
#     if aa == 'org':
#         pass
#     else:
#         plt.figure(figsize=(7,7))
#         st = int(poloidal_index_list[0]) -1
#         ed = int(poloidal_index_list[-1])
#         x = xl.data['iout_data']['poloidal_flux_change_percent'][aa][19, st:ed]
#         y = xl.data['neuden_change'][aa]
#         plt.scatter(x, y)
#         plt.title('{} percent change correlation'.format(aa))












    


