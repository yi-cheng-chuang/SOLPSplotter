# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:34:00 2023

@author: user
"""
import matplotlib.pyplot as plt


def opacity_plot(pol_loc, result_dic, unit_dic, log_flag, charactor, 
                 iter_list, change_ver_dic):
    withshift = charactor['withshift']
    withseries = charactor['withseries']
    # j = 1
    adj_list = list(result_dic.keys())
    # print(adj_list)
    adj_list.append('method2_fitting_width_error')
    # adj_list.append('width_relation')
    # print(adj_list)
    for i in adj_list:
        plt.figure(figsize=(7,7))
        if log_flag:
            plt.yscale('log')
        else:
            pass
        if withshift == False and withseries == False:
            if i == 'efold_length_method2':
                plt.plot(pol_loc, result_dic['efold_length_method2'], 'o--', 
                         label= '{}'.format(i))
                plt.title('{} verses poloidal index from {} to {}'.format(i, 
                                                         pol_loc[0], pol_loc[-1]))
                plt.xlabel('poloidal angle')
                
            elif i == 'std_m2':
                plt.errorbar(pol_loc, result_dic['efold_length_method2'], yerr = result_dic[i], fmt = 'o--', 
                         label= '{} with error'.format('efold_length_method2'))
                plt.title('{} verses poloidal index from {} to {}'.format('efold_length_method2', 
                                                          pol_loc[0], pol_loc[-1]))
                plt.xlabel('poloidal angle')
                
            elif i == 'electron_pedestal_density':
                plt.scatter(result_dic['electron_pedestal_density'], result_dic['efold_length_method2'], 
                         label= '{}'.format('efold_length_method2'))
                plt.scatter(result_dic['electron_pedestal_density'], result_dic['pedestal_width'], 
                         label= '{}'.format('pedestal_width'))
                plt.xlabel('electron_pedestal_density: $n_{ped}$ (m$^{-3}$)')
            
            
            elif i == 'method2_fitting_width_error':
                plt.scatter(result_dic['method2_fitting_width'], result_dic['std_m2'], 
                         label= '{}'.format('method2_fitting_width'))
                plt.xlabel('method2_fitting_width')
                plt.title('{} verses {} from {} to {} for different modified distances'.format('method2_fitting_error', 
                                                         'method2_fitting_width', pol_loc[0], pol_loc[-1]))
                      
            else:
                plt.plot(pol_loc, result_dic[i],'o--', 
                     label= '{}'.format(i))
                plt.title('{} verses poloidal index from {} to {}'.format(i, 
                                                     pol_loc[0], pol_loc[-1]))
                plt.xlabel('poloidal angle')
            # plt.ylabel('{}'.format(unit_dic[i]))
            if max(pol_loc) > 90 and i != 'electron_pedestal_density' and i != 'method2_fitting_width_error' and i != 'width_relation':
                plt.axvline(x= 90, color='black',lw=3, ls='--')
            else:
                pass
        
        elif withshift == True and withseries == False:
            for aa in iter_list:
                if i == 'efold_length_method2':
                    plt.plot(pol_loc[aa], result_dic['efold_length_method2'][aa], 'o--', 
                             label= '{}'.format(i))
                    plt.title('{} verses poloidal index from {} to {} for different modified distances'.format(i, 
                                                             pol_loc[aa][0], pol_loc[aa][-1]))
                    plt.xlabel('poloidal angle')
                    
                    # print(pol_loc[aa])
                    # print(result_dic['efold_length_method2'][aa])
                    
                elif i == 'std_m2':
                    plt.errorbar(pol_loc[aa], result_dic['efold_length_method2'][aa], yerr = result_dic[i][aa], fmt = 'o--', 
                             label= '{} with error'.format('efold_length_method2'))
                    plt.title('{} verses poloidal index from {} to {} for different modified distances'.format('efold_length_method2', 
                                                             pol_loc[aa][0], pol_loc[aa][-1]))
                    plt.xlabel('poloidal angle')
                    
                elif i == 'electron_pedestal_density':
                    plt.scatter(result_dic['electron_pedestal_density'][aa], result_dic['efold_length_method2'][aa], 
                             label= '{} {}'.format(aa, 'efold_length_method2'))
                    plt.scatter(result_dic['electron_pedestal_density'][aa], result_dic['pedestal_width'][aa], 
                             label= '{} {}'.format(aa, 'pedestal_width'))
                    plt.xlabel('electron_pedestal_density: $n_{ped}$ (m$^{-3}$)')
                    plt.title('{} and {} verses {} from {} to {} for different modified distances'.format('efold_length_method2', 
                                                             'pedestal_width', 'electron_pedestal_density', pol_loc[aa][0], pol_loc[aa][-1]))
                
                elif i == 'method2_fitting_width_error':
                    plt.scatter(result_dic['method2_fitting_width'][aa], result_dic['std_m2'][aa], 
                             label= '{} {}'.format(aa, 'method2_fitting_width'))
                    plt.xlabel('method2_fitting_width')
                    plt.title('{} verses {} from {} to {} for different modified distances'.format('method2_fitting_error', 
                                                             'method2_fitting_width', pol_loc[aa][0], pol_loc[aa][-1]))                    
                
                
                # elif i == 'width_relation':
                #     a_w = result_dic['method2_fitting_width'][aa]
                #     x = [min(a_w), max(a_w)]
                #     plt.scatter(result_dic['method2_fitting_width'][aa], result_dic['pedestal_width'][aa], 
                #              label= '{} {}'.format(aa, 'width_relation'))
                #     # plt.plot(x, x, label= 'x=y straight line')
                #     plt.xlabel('method2_fitting_width')
                #     plt.title('{} verses {} from {} to {} for different modified distances'.format('method2_fitting_width', 
                #                                                 'pedestal_width', pol_loc[0], pol_loc[-1]))    
                    
                else:
                    plt.plot(pol_loc[aa], result_dic[i][aa], 'o--',
                         label= '{} for modified {} m case'.format(i, change_ver_dic[aa]))
                    plt.title('{} verses poloidal index from {} to {} for different modified distances'.format(i, 
                                                     pol_loc[aa][0], pol_loc[aa][-1]))
                    plt.xlabel('poloidal angle')
            if max(pol_loc[aa]) > 90 and i != 'electron_pedestal_density' and i != 'method2_fitting_width_error' and i != 'width_relation':
                plt.axvline(x= 90, color='black',lw=3, ls='--')
            else:
                pass
            # plt.ylabel('{}'.format(unit_dic[i]))
            
        elif withshift == False and withseries == True:
            for aa in iter_list:
                if i == 'efold_length_method2':
                    plt.plot(pol_loc, result_dic['efold_length_method2'][aa], 'o--', 
                             label= '{}'.format(i))
                    plt.title('{} verses poloidal index from {} to {} for different core electron density'.format(i, 
                                                             pol_loc[0], pol_loc[-1]))
                    plt.xlabel('poloidal angle')
                    
                elif i == 'std_m2':
                    plt.errorbar(pol_loc, result_dic['efold_length_method2'][aa], yerr = result_dic[i][aa], fmt = 'o--', 
                             label= '{} with error'.format('efold_length_method2'))
                    plt.title('{} verses poloidal index from {} to {} for different core electron density'.format('efold_length_method2', 
                                                             pol_loc[0], pol_loc[-1]))
                    plt.xlabel('poloidal angle')
                
                elif i == 'electron_pedestal_density':
                    plt.scatter(result_dic['electron_pedestal_density'][aa], result_dic['efold_length_method2'][aa], 
                              label= '{} {}'.format(aa, 'efold_length_method2'))
                    plt.scatter(result_dic['electron_pedestal_density'][aa], result_dic['pedestal_width'][aa], 
                              label= '{} {}'.format(aa, 'pedestal_width'))
                    plt.xlabel('electron_pedestal_density: $n_{ped}$ (m$^{-3}$)')
                    plt.title('{} and {} verses {} from {} to {} for different core electron density'.format('efold_length_method2', 
                                                             'pedestal_width', 'electron_pedestal_density', pol_loc[0], pol_loc[-1]))
                
                elif i == 'method2_fitting_width_error':
                    plt.scatter(result_dic['method2_fitting_width'][aa], result_dic['std_m2'][aa], 
                             label= '{} {}'.format(aa, 'method2_fitting_width'))
                    plt.xlabel('method2_fitting_width')
                    plt.title('{} verses {} from {} to {} for different core electron density'.format('method2_fitting_error', 
                                                             'method2_fitting_width', pol_loc[0], pol_loc[-1]))                    
                    
                else:
                    plt.plot(pol_loc, result_dic[i][aa],'o--', 
                         label= '{} and electron density is {}'.format(i, change_ver_dic[aa]))
                    plt.title('{} verses poloidal index from {} to {} for different core electron density'.format(i, 
                                                     pol_loc[0], pol_loc[-1]))
                    plt.xlabel('poloidal angle')
            # plt.ylabel('{}'.format(unit_dic[i]))
            if max(pol_loc) > 90 and i != 'electron_pedestal_density' and i != 'method2_fitting_width_error' and i != 'width_relation':
                plt.axvline(x= 90, color='black',lw=3, ls='--')
            else:
                pass
        
        elif withshift == True and withseries == True:
            print('opacity_plot is not there yet, to be continue...')
        
        else:
            print('opacity_plot has a bug')
        
        
        
        plt.legend()
        # j = j + 1


def data_reorder(iter_list, change_var, data_collect, char):
    withshift = char['withshift']
    withseries = char['withseries']
    
    if withshift == False and withseries == True:
        plt.figure(figsize=(7,7))
        color_list = ['r', 'b', 'g', 'c', 'k']
        for p in range(len(change_var)):
            for l in range(len(iter_list)):
                plt.scatter(change_var[p], data_collect[l, p], color= color_list[p], 
                            label= 'dimensionless opaqueness')
        plt.xlabel('Electron density: ${n_e}$ (m$^{-3}$)')
        # plt.ylabel('dimensionless opaqueness')
        plt.title('opaqueness verses shift core electron density')
        
    elif withshift == True and withseries == False:
        plt.figure(figsize=(7,7))
        color_list = ['r', 'b', 'g', 'c', 'k']
        for p in range(len(change_var)):
            for l in range(len(iter_list)):
                plt.scatter(change_var[p], data_collect[l, p], color= color_list[p], 
                            label= 'dimensionless opaqueness')
        plt.xlabel('modify distance: [m]')
        # plt.ylabel('dimensionless opaqueness')
        plt.title('opaqueness verses different modify distance')
    
    
    # ii = ii + 1
    # plt.figure(ii)
    # for m in iter_list:
    #     plt.plot(change_var_dic[m], efold_dic[m],'o-', color= 'r', label= 'neutral penetration')
    # plt.xlabel('shift: [m]')
    # plt.ylabel('neutral penetration length: [m]')
    # plt.title('neutral penetration length verses shift distance')
    
    # ii = ii + 1
    # plt.figure(ii)
    # for p in iter_list:
    #     plt.plot(change_var_dic[p], delta_dic[p],'o-', color= 'r', label= 'neutral penetration')
    # plt.xlabel('shift: [m]')
    # plt.ylabel('pedestal width: [m]')
    # plt.title('pedestal width verses shift distance')





def opacity_radial_method_single(result_dic, SEP, x_choice, x_coord, Nd, Ne, Te, 
                                      P, log_flag):
    
    tanh_ne_fit = result_dic['tanh_ne_fit']
    tanh_te_fit = result_dic['tanh_te_fit']
    exp_an_fit = result_dic['exp_fit']
    dn = result_dic['pedestal_width']
    dtn = result_dic['temperature_pedestal_width']
    efold = result_dic['efold_length']
    opq = result_dic['dimensionless_opaqueness']
    xcoord_cut = result_dic['x_coord_cut']
    sym_pt = result_dic['ne_symmetry_point']
    te_sym_pt = result_dic['te_symmetry_point']
    
    exp_fit_m2 = result_dic['exp_fit_m2']
    x_m2 = result_dic['x_m2']
    efold_m2 = result_dic['efold_length_method2']
    opq_m2 = result_dic['dimensionless_opaqueness_method2']
    # print(x_m2)
       
       
    x = [-efold_m2 + max(xcoord_cut), max(xcoord_cut)]
    y = [min(exp_fit_m2), min(exp_fit_m2)]
    xd = [-dn + sym_pt, dn + sym_pt]
    yd = [tanh_ne_fit[SEP] , tanh_ne_fit[SEP]]
    xt = [-dtn + te_sym_pt, dtn + te_sym_pt]
    yt = [tanh_te_fit[SEP], tanh_te_fit[SEP]]
    
    
    
    plt.figure(figsize=(7,7))
    if log_flag:
        plt.yscale('log')
    else:
        pass
    plt.plot(x_coord, Nd,'o-', color = 'green', label= 'solps neutral density')
    # plt.plot(psi_RGI, Nd,'o-', color = 'b', label= 'RGI_solps neutral density')
    plt.plot(xcoord_cut, exp_an_fit, color='r',lw= 5, label= 'exponential fit')
    # plt.plot(x_m2, exp_fit_m2, color='b',lw= 5, label= 'exponential fit m2')
    # plt.plot(x_m3, exp_fit_m3, color='cyan',lw= 5, label= 'exponential fit m3')
    # plt.plot(exp_dsa, exp_fit, color='r',lw= 3, label= 'exponential fit')
    plt.axvline(x=max(xcoord_cut), color='orange',lw=3)
    plt.plot(x,y, color='orange', lw=3, label= 'Neutral penetration length [m]: $\lambda_{n_D}$')
    plt.axvline(x=-efold_m2 + max(xcoord_cut), color='orange',lw=3)
    plt.axvline(x=dn + sym_pt, color='black',lw=3, ls='--', 
                label= 'Pedestal width [m]: $\Delta n_e$')
    plt.axvline(x=-dn + sym_pt, color='black',lw=3, ls='--')
    # plt.axvline(x= x_m2[0], color='purple',lw=3, ls='--', 
    #             label= 'exp fitting width')
    # plt.axvline(x= x_m2[-1], color='purple',lw=3, ls='--')
    if x_choice == 'psiN':
        plt.xlabel('psiN')
    elif x_choice == 'RRsep':
        plt.xlabel('Radial coordinate: $R- R_{sep}$')
    # plt.ylabel(P['NeuDen'])
    plt.title('Neutral density with fits')
    plt.legend()
        
    # plt.subplot(211, sharex= ax1)
    plt.figure(figsize=(7,7))
    if log_flag:
        plt.yscale('log')
    else:
        pass
    # plt.plot(psi_xport, Ne,'o-', color = 'r', label= 'solps_electron density')
    plt.plot(x_coord, Ne,'o-', color = 'b', label= 'solps electron density')
    # plt.plot(fitdsa, cutNe,'o-', color = 'g', label= 'experiment electron density')
    plt.plot(x_coord, tanh_ne_fit, color='r',lw= 3, label= 'tanh fit')
    plt.plot(xd, yd, color='black', lw=3, label= 'Pedestal width [m]: $\Delta n_e$')
    plt.axvline(x=dn + sym_pt, color='black',lw=3)
    plt.axvline(x=-dn + sym_pt, color='black',lw=3)
    plt.axvline(x=max(xcoord_cut), color='orange',lw=3, ls='--', 
                label = 'Neutral penetration length [m]: $\lambda_{n_D}$')
    plt.axvline(x=-efold_m2 + max(xcoord_cut), color='orange',lw=3, ls='--')
    if x_choice == 'psiN':
        plt.xlabel('psiN')
    elif x_choice == 'RRsep':
        plt.xlabel('Radial coordinate: $R- R_{sep}$')
    # plt.ylabel(P['Ne'])
    plt.title('Electron density with fits')
    # plt.title(plot_dic['ne3da.last10'][0],fontdict={"family":"Calibri","size": 20})
    plt.legend()
    
    plt.figure(figsize=(7,7))
    if log_flag:
        plt.yscale('log')
    else:
        pass
    # plt.plot(psi_xport, Ne,'o-', color = 'r', label= 'solps_electron density')
    plt.plot(x_coord, Te,'o-', color = 'b', label= 'solps electron tempurature')
    # plt.plot(fitdsa, cutNe,'o-', color = 'g', label= 'experiment electron density')
    plt.plot(x_coord, tanh_te_fit, color='r',lw= 3, label= 'tanh fit')
    plt.plot(xt, yt, color='black', lw=3, label= 'temperature pedestal width [m]: $\Delta n_e$')
    plt.axvline(x=dtn + te_sym_pt, color='black',lw=3)
    plt.axvline(x=-dtn + te_sym_pt, color='black',lw=3)
    if x_choice == 'psiN':
        plt.xlabel('psiN')
    elif x_choice == 'RRsep':
        plt.xlabel('Radial coordinate: $R- R_{sep}$')
    # plt.ylabel(P['Te'])
    plt.title('Electron temperature with fits')
    plt.legend()
    
    # plt.figure(figsize=(7,7))
    # if log_flag:
    #     plt.yscale('log')
    # else:
    #     pass
    # plt.plot(dsa, psi,'o-', color = 'green', label= 'solps neutral density')
    # plt.plot(dsa_psn, psi_psn, color='r',lw= 5, label= 'exponential fit')
    # plt.axvline(x=dn, color='black',lw=3, ls= '--', label= 'Pedestal width [m]: $\Delta n_e$')
    # plt.axvline(x=-dn, color='black',lw=3, ls= '--')
    # plt.xlabel('RR_sep')
    # plt.ylabel('Magnetic flux coordinate: ${\psi_N}$')
    # plt.title('dsa_psi_fit')
    # plt.legend()






def opacity_radial_method_multi(x_choice, result_dic, SEP, iter_list, 
                                change_var_dic, log_flag, char, P):
    
    tanh_ne_fit_dic = result_dic['tanh_ne_fit']
    tanh_te_fit_dic = result_dic['tanh_te_fit']
    exp_an_fit_dic = result_dic['exp_fit']
    delta_dic = result_dic['pedestal_width']
    tdelta_dic = result_dic['temperature_pedestal_width']
    efold_dic = result_dic['efold_length']
    opq_dic = result_dic['dimensionless_opaqueness']
    xcoord_cut_dic = result_dic['xcoord_cut']
    Nd_dic = result_dic['neutral_density']
    Ne_dic = result_dic['electron_density']
    Te_dic = result_dic['electron_temperature_density']
    psi_dic = result_dic['psiN']
    dsa_dic = result_dic['dsa']
    te_sym_pt_dic = result_dic['te_symmetry_point']
    ne_sym_pt_dic = result_dic['ne_symmetry_point']
    
    exp_fit_m2_dic = result_dic['exp_fit_m2']
    x_m2_dic = result_dic['x_m2']
    efold_m2_dic = result_dic['efold_length_method2']
    
    
    withshift = char['withshift']
    withseries = char['withseries'] 
    
    
    # ii = 1
    
    for i in iter_list:
        plt.figure(figsize=(7,7))
        xp = [-efold_m2_dic[i] + max(xcoord_cut_dic[i]), max(xcoord_cut_dic[i])]
        yp = [min(exp_an_fit_dic[i]), min(exp_an_fit_dic[i])]
        if x_choice == 'psiN':
            plt.plot(psi_dic[i], Nd_dic[i],'o-', color = 'green', label= 'solps neutral density')
            plt.plot(xcoord_cut_dic[i], exp_an_fit_dic[i], color='r',lw= 5, label= 'exponential fit')
            # plt.plot(x_m2_dic[i], exp_fit_m2_dic[i], color='blue',lw= 5, label= 'exponential fit m2')
            # plt.plot(x_m3_dic[i], exp_fit_m3_dic[i], color='cyan',lw= 5, label= 'exponential fit m3')
        elif x_choice == 'RRsep':
            plt.plot(dsa_dic[i], Nd_dic[i],'o-', color = 'green', label= 'solps neutral density')
            plt.plot(xcoord_cut_dic[i], exp_an_fit_dic[i], color='r',lw= 5, label= 'exponential fit')
        plt.axvline(x=max(xcoord_cut_dic[i]), color='orange',lw=3)
        plt.plot(xp, yp, color='orange', lw=3, label= 'Neutral penetration length: $\lambda_{n_D}$')
        plt.axvline(x = xp[0], color='orange',lw=3)
        plt.axvline(x= delta_dic[i] + ne_sym_pt_dic[i], color='black',lw=3, ls='--', label= 'Pedestal width')
        plt.axvline(x=-delta_dic[i] + ne_sym_pt_dic[i], color='black',lw=3, ls='--')
        # plt.axvline(x= x_m2_dic[i][0], color='purple',lw=3, ls='--', label= 'exp fitting width')
        # plt.axvline(x= x_m2_dic[i][-1], color='purple',lw=3, ls='--')
        if x_choice == 'psiN':
            plt.xlabel('psiN')
        elif x_choice == 'RRsep':
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
        # plt.ylabel(P['NeuDen'])
        shift_value = str(change_var_dic[i])
        plt.title('Modify {} m Neutral density with fits'.format(shift_value))
        # plt.title(plot_dic['an3da.last10'][0],fontdict={"family":"Calibri","size": 20})
        plt.legend()
        # ii = ii + 1
        if log_flag:
            plt.yscale('log')
        else:
            pass
    
    # plt.plot(psi_xport, Ne,'o-', color = 'r', label= 'solps_electron density')
    for j in iter_list:
        plt.figure(figsize=(7,7))
        if withshift == True and withseries == False:
            sep = int(SEP[j])
        elif withshift == False and withseries == True:
            sep = int(SEP)
        xd = [-delta_dic[j] + ne_sym_pt_dic[j], delta_dic[j] + ne_sym_pt_dic[j]]
        yd = [tanh_ne_fit_dic[j][sep], tanh_ne_fit_dic[j][sep]]
        if x_choice == 'psiN':
            plt.plot(psi_dic[j], Ne_dic[j],'o-', label= 'solps electron density_{}'.format(j))
            plt.plot(psi_dic[j], tanh_ne_fit_dic[j], color='r',lw= 3, label= 'exponential fit')
        elif x_choice == 'RRsep':
            plt.plot(dsa_dic[j], Ne_dic[j],'o-', label= 'solps electron density_{}'.format(j))
            plt.plot(dsa_dic[j], tanh_ne_fit_dic[j], color='r',lw= 3, label= 'exponential fit')  
        plt.plot(xd, yd, color='black', lw=3, label= 'Pedestal width [m]: $\Delta n_e$')
        plt.axvline(x= delta_dic[j] + ne_sym_pt_dic[j], color='black',lw=3)
        plt.axvline(x=-delta_dic[j] + ne_sym_pt_dic[j], color='black',lw=3)
        plt.axvline(x=max(xcoord_cut_dic[j]), color='orange',lw=3, ls='--')
        plt.axvline(x=-efold_m2_dic[j] + max(xcoord_cut_dic[j]), color='orange',lw=3, ls='--')
        if x_choice == 'psiN':
            plt.xlabel('psiN')
        elif x_choice == 'RRsep':
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
        # plt.ylabel(P['Ne'])
        plt.title('Electron density with fits')
        plt.legend()
    
    for k in iter_list:
        plt.figure(figsize=(7,7))
        if withshift == True and withseries == False:
            sep = int(SEP[k])
        elif withshift == False and withseries == True:
            sep = int(SEP)
        xt = [-tdelta_dic[k] + te_sym_pt_dic[k], tdelta_dic[k] + te_sym_pt_dic[k]]
        yt = [tanh_te_fit_dic[k][sep], tanh_te_fit_dic[k][sep]]
        # plt.plot(psi_xport, Ne,'o-', color = 'r', label= 'solps_electron density')  
        # plt.plot(fitdsa, cutNe,'o-', color = 'g', label= 'experiment electron density')
        if x_choice == 'psiN': 
            plt.plot(psi_dic[k], Te_dic[k],'o-', color = 'b', label= 'solps electron tempurature')
            plt.plot(psi_dic[k], tanh_te_fit_dic[k], color='r',lw= 3, label= 'tanh fit')
        if x_choice == 'RRsep':
            plt.plot(dsa_dic[k], Te_dic[k],'o-', color = 'b', label= 'solps electron tempurature')
            plt.plot(dsa_dic[k], tanh_te_fit_dic[k], color='r',lw= 3, label= 'tanh fit')
        plt.plot(xt, yt, color='black', lw=3, label= 'temperature pedestal width [m]: $\Delta n_e$')
        plt.axvline(x=tdelta_dic[k] + te_sym_pt_dic[k], color='black',lw=3)
        plt.axvline(x=-tdelta_dic[k] + te_sym_pt_dic[k], color='black',lw=3)
        if x_choice == 'psiN':
            plt.xlabel('psiN')
        elif x_choice == 'RRsep':
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
        # plt.ylabel(P['Te'])
        plt.title('Electron temperature with fits')
        plt.legend()
    
    
    # for m in iter_list:      
    #     plt.figure(figsize=(7,7))
    #     plt.plot(dsa_dic[m], psi_dic[m],'o-', color = 'green', label= 'solps neutral density')
    #     plt.plot(dsa_psn_dic[m], psi_psn_dic[m], color='r',lw= 5, label= 'exponential fit')
    #     plt.axvline(x=delta_dic[m], color='black',lw=3, ls= '--', label= 'Pedestal width [m]: $\Delta n_e$')
    #     plt.axvline(x=-delta_dic[m], color='black',lw=3, ls= '--')
    #     plt.xlabel('Radial coordinate: $R- R_{sep}$')
    #     plt.ylabel('Magnetic flux coordinate: ${\psi_N}$')
    #     plt.title('dsa_psi_fit')
    #     plt.legend()
    
    # plt.close()

    
    



