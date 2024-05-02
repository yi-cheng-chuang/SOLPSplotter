# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:34:00 2023

@author: user
"""
import matplotlib.pyplot as plt
import numpy as np


def opacity_study_unit():
    unit = {'efold_length_psiN': 'Neutral penetration length ($\psi_N$)',
            'pedestal_width_psiN': 'Pedestal width ($\psi_N$)',
              'dimensionless_opaqueness': 'Experimental opaqueness', 
              'neutral_density': 'Neutral density ${n_D}$ (m$^{-3}$)', 
              'electron_pedestal_density': 'Electron pedestal density: $n_{ped}$ (m$^{-3}$)',
              'temperature_pedestal_width': 'Temperature pedestal width: $\Delta T$',
              'flux_expansion': 'Flux expansion',
              'efold_length': '$\lambda_{n_D}$ [mm]',
              'pedestal_width': '$\Delta n_e$ [mm]',
              
              }
    return unit


def opacity_plot(pol_loc, result_dic, unit_dic, log_flag, charactor, 
                 iter_list, change_ver_dic, xpoint_loc):
    withshift = charactor['withshift']
    withseries = charactor['withseries']
    # j = 1
    adj_list = list(result_dic.keys())
    
    # print(adj_list)
    for i in adj_list:
        plt.figure(figsize=(7,7))
        if log_flag:
            plt.yscale('log')
        else:
            pass
        if withshift == False and withseries == False:              
            if i == 'electron_pedestal_density':
                plt.scatter(result_dic['electron_pedestal_density'], result_dic['efold_length'], 
                         label= '{}'.format('efold_length'))
                plt.scatter(result_dic['electron_pedestal_density'], result_dic['pedestal_width'], 
                         label= '{}'.format('pedestal_width'))
                plt.xlabel('electron_pedestal_density: $n_{ped}$ (m$^{-3}$)')
                
            elif i == 'pedestal_width_psiN' or i == 'temperature_pedestal_width':
                plt.errorbar(pol_loc, np.round_(result_dic[i], 2), yerr= np.std(result_dic[i]), fmt= 'o', 
                             label= '{}'.format(i))
            else:
                plt.plot(pol_loc, result_dic[i],'o', 
                     label= '{}'.format(i))
                plt.title('{} verses poloidal index from {:.2f} to {:.2f}'.format(i, 
                                                     pol_loc[0], pol_loc[-1]))
                plt.xlabel('poloidal angle')
            # plt.ylabel('{}'.format(unit_dic[i]))
            if max(pol_loc) > 90 and i != 'electron_pedestal_density' and i != 'width_relation':
                plt.axvline(x= 90, color='black',lw=3, ls='--', label= '90')
            else:
                pass
            if max(pol_loc) > 180 and i != 'electron_pedestal_density':
                plt.axvline(x= 180, color='deeppink',lw=3, ls='--', label= 'inner midplane')
            else:
                pass
            if min(pol_loc) < 0 and i != 'electron_pedestal_density':
                plt.axvline(x= 0, color='chocolate',lw=3, ls='--', label= 'outer midplane')
            else:
                pass
            if min(pol_loc) < -90 and i != 'electron_pedestal_density':
                plt.axvline(x= xpoint_loc, color='magenta',lw=3, ls='--', label= 'xpoint angle')
                plt.axvline(x= xpoint_loc + 360, color='magenta',lw=3, ls='--')
            else:
                pass
        
        elif withshift == True and withseries == False:
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            for aa in iter_list:
                if i == 'electron_pedestal_density':
                    plt.scatter(result_dic['electron_pedestal_density'][aa], result_dic['efold_length'][aa], 
                             color= color_dic[aa], label= 'aspect ratio = {} {}'.format(aa, 'efold length'))
                    plt.scatter(result_dic['electron_pedestal_density'][aa], result_dic['pedestal_width'][aa], 
                             color= color_dic[aa], label= 'aspect ratio = {} {}'.format(aa, 'pedestal width'))
                    plt.xlabel('electron_pedestal_density: $n_{ped}$ (m$^{-3}$)')
                    plt.title('{} and {} verses {} from {:.2f} to {:.2f}'.format('efold length', 
                                                              'pedestal width', 'electron pedestal density', max(pol_loc[aa]), min(pol_loc[aa])))
                    # plt.title('{} and {} verses {}'.format('efold length', 
                    #                 'pedestal width', 'electron pedestal density'))
                    plt.legend()
                    
                elif i == 'pedestal_width_psiN':
                    # plt.scatter(pol_loc[aa], np.round_(result_dic[i][aa], 2), 
                    #              label= 'modified {} m case'.format(change_ver_dic[aa]))
                    # plt.errorbar(pol_loc[aa], np.round_(result_dic[i][aa], 2), yerr= np.std(result_dic[i][aa]), fmt= 'o', 
                    #               label= '{}'.format(i))
                    plt.plot(pol_loc[aa], np.round_(result_dic[i][aa], 2), '-', color= color_dic[aa],
                          label= 'aspect ratio = {}'.format(change_ver_dic[aa]))
                    plt.title('{} verses poloidal angle from {:.2f} to {:.2f}'.format(unit_dic[i], 
                                                      max(pol_loc[aa]), min(pol_loc[aa])))
                    plt.xlabel('poloidal angle')
                    plt.legend()

                elif i == 'pedestal_width' or i == 'efold_length':
                    new_r = result_dic[i][aa]*1000
                    plt.plot(pol_loc[aa], new_r, '-', color= color_dic[aa], 
                             label= 'aspect ratio = {}'.format(change_ver_dic[aa]))
                    plt.title('{} verses poloidal angle from {:.2f} to {:.2f}'.format(unit_dic[i], 
                                                      max(pol_loc[aa]), min(pol_loc[aa])))                   
                    plt.xlabel('poloidal angle')
                    plt.legend()
                    
                else:
                    plt.plot(pol_loc[aa], result_dic[i][aa], '-', color= color_dic[aa] ,
                          label= 'aspect ratio = {}'.format(change_ver_dic[aa]))
                    # plt.plot(pol_loc[aa], result_dic[i][aa], '-', color= color_dic[aa])
                    plt.title('{} verses poloidal angle from {:.2f} to {:.2f}'.format(unit_dic[i], 
                                                      max(pol_loc[aa]), min(pol_loc[aa])))
                    # plt.title('{} verses poloidal angle'.format(unit_dic[i]))
                    # plt.title('{}'.format(unit_dic[i]))                   
                    plt.xlabel('poloidal angle')
                    plt.legend()

            if max(pol_loc['org']) > 90 and i != 'electron_pedestal_density' and i != 'width_relation':
                plt.axvline(x= 90, color='black',lw=3, ls='--', label= '90')
            else:
                pass
            if max(pol_loc['org']) > 180 and i != 'electron_pedestal_density':
                plt.axvline(x= 180, color='seagreen',lw=3, ls='--', label= 'inner midplane')
            else:
                pass
            if min(pol_loc['org']) < 0 and i != 'electron_pedestal_density':
                plt.axvline(x= 0, color='darkorange',lw=3, ls='--', label= 'outer midplane')
            else:
                pass
            
            if min(pol_loc['org']) < -70 and i != 'electron_pedestal_density':
                plt.axvline(x= xpoint_loc[aa], color='darkblue',lw=3, ls='--', label= 'xpoint angle')
                plt.axvline(x= xpoint_loc[aa] + 360, color='darkblue',lw=3, ls='--')
            else:
                pass
            # plt.ylabel('{}'.format(unit_dic[i]))
            plt.legend()
            
        elif withshift == False and withseries == True:
            for aa in iter_list:
                if i == 'electron_pedestal_density':
                    plt.scatter(result_dic['electron_pedestal_density'][aa], result_dic['efold_length'][aa], 
                              label= '{} {}'.format(aa, 'efold_length'))
                    plt.scatter(result_dic['electron_pedestal_density'][aa], result_dic['pedestal_width'][aa], 
                              label= '{} {}'.format(aa, 'pedestal_width'))
                    plt.xlabel('electron_pedestal_density: $n_{ped}$ (m$^{-3}$)')
                    plt.title('{} and {} verses {} from {:.2f} to {:.2f}'.format('efold_length', 
                                                             'pedestal_width', 'electron_pedestal_density', pol_loc[0], pol_loc[-1]))
                    
                else:
                    plt.plot(pol_loc, result_dic[i][aa],'o', 
                         label= '{} and electron density is {}'.format(i, change_ver_dic[aa]))
                    plt.title('{} verses poloidal index from {:.2f} to {:.2f}'.format(unit_dic[i], 
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
        
        

def data_reorder(iter_list, change_var, data_collect, char):
    withshift = char['withshift']
    withseries = char['withseries']
    series_flag = char['series_flag']
    
    if withshift == False and withseries == True:
        plt.figure(figsize=(7,7))
        color_list = ['red', 'salmon', 'orange', 'lime', 'green', 'darkgreen', 
                      'cyan', 'deepskyblue', 'navy', 'purple']
        if series_flag == 'change_den':
            for p in range(len(change_var)):
                x_cor = np.ones(len(iter_list))*change_var[p]
                plt.scatter(x_cor, data_collect[:, p], color= color_list[p], 
                                label= 'efold length for different core electron density {}'.format(change_var[p]))
            plt.xlabel('Electron density: ${n_e}$ (m$^{-3}$)')
            plt.title('dimensionless opaqueness verses shift core electron density')
            plt.legend()
            
        elif series_flag == 'eireneN':
            for p in range(len(change_var)):
                x_cor = np.ones(len(iter_list))*change_var[p]
                plt.scatter(x_cor, data_collect[:, p], color= color_list[p], 
                                label= '{}'.format(change_var[p]))
            plt.xlabel('Eirene SOL particle number')
            # plt.yscale('log')
            plt.xscale('log')
            # plt.title('dimensionless opaqueness verses different SOL eirene particle number')
            plt.title('dimensionless opaqueness')
            # plt.legend()
        
    elif withshift == True and withseries == False:
        plt.figure(figsize=(7,7))
        A_list = [1.4, 2.0, 2.4, 2.8, 3.4]
        color_list = ['red', 'orange', 'green', 'blue', 'purple']
        for p in range(len(change_var)):
            x_cor = np.ones(len(iter_list))*A_list[p]
            plt.scatter(x_cor, data_collect[:, p], color= color_list[p], 
                        label= 'aspect ratio= {}'.format(str(A_list[p])))
        plt.axvline(x= 0.7/0.5, color='salmon',lw=3, ls='--', 
                    label= 'MAST aspect ratio')
        plt.axvline(x= 1.67/0.67, color='darkgreen',lw=3, ls='--', label= 'D3D aspect ratio')
        plt.axvline(x= 0.68/0.22, color='cyan',lw=3, ls='--', label= 'C-Mod/ ITER aspect ratio')
        # plt.axvline(x= 3.4, color='black',lw=3, ls='--', label= 'JT-60 aspect ratio')
        # plt.xlabel('aspect ratio')
        # plt.ylabel('dimensionless opaqueness')
        # plt.title('dimensionless opaqueness verses different modify distance')
        plt.title('dimensionless opaqueness')
        plt.legend()
    






def opacity_radial_method_single(result_dic, SEP, x_coord, Nd, Ne, Te, 
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
    
       
       
    x = [-efold + max(xcoord_cut), max(xcoord_cut)]
    y = [min(exp_an_fit), min(exp_an_fit)]
    xd = [-dn + sym_pt, dn + sym_pt]
    yd = [tanh_ne_fit[SEP] , tanh_ne_fit[SEP]]
    xt = [-dtn + te_sym_pt, dtn + te_sym_pt]
    yt = [tanh_te_fit[SEP], tanh_te_fit[SEP]]
    
    
    
    plt.figure(figsize=(7,7))
    if log_flag:
        plt.yscale('log')
    else:
        pass
    plt.plot(x_coord, Nd,'-', color = 'green', label= 'solps neutral density')
    # plt.plot(psi_RGI, Nd,'-', color = 'b', label= 'RGI_solps neutral density')
    plt.plot(xcoord_cut, exp_an_fit, color='r',lw= 5, ls='-', label= 'exponential fit')
    plt.axvline(x=max(xcoord_cut), color='orange',lw=3)
    plt.plot(x,y, color='orange', lw=3, ls='-', label= 'Neutral penetration length : $\lambda_{n_D}$')
    plt.axvline(x=-efold + max(xcoord_cut), color='orange',lw=3)
    plt.axvline(x= max(xcoord_cut), color='black',lw=3, ls='--', 
                label= 'fit range : $\Delta n_e$')
    plt.axvline(x= min(xcoord_cut), color='black',lw=3, ls='--')
    # plt.axvline(x= x_m2[0], color='purple',lw=3, ls='--', 
    #             label= 'exp fitting width')
    # plt.axvline(x= x_m2[-1], color='purple',lw=3, ls='--')
    plt.xlabel('psiN')
    # plt.ylabel(P['NeuDen'])
    plt.title('Neutral density with fits')
    plt.legend()
        
    # plt.subplot(211, sharex= ax1)
    plt.figure(figsize=(7,7))
    if log_flag:
        plt.yscale('log')
    else:
        pass
    # plt.plot(psi_xport, Ne,'-', color = 'r', label= 'solps_electron density')
    plt.plot(x_coord, Ne,'-', color = 'b', label= 'solps electron density')
    # plt.plot(fitdsa, cutNe,'-', color = 'g', label= 'experiment electron density')
    plt.plot(x_coord, tanh_ne_fit, ls='-', color='r',lw= 3, label= 'tanh fit')
    plt.plot(xd, yd, color='black', ls='-', lw=3, label= 'Pedestal width : $\Delta n_e$')
    plt.axvline(x=dn + sym_pt, color='black',lw=3)
    plt.axvline(x=-dn + sym_pt, color='black',lw=3)
    plt.axvline(x=max(xcoord_cut), color='orange',lw=3, ls='--', 
                label = 'Neutral penetration length : $\lambda_{n_D}$')
    plt.axvline(x=-efold + max(xcoord_cut), color='orange',lw=3, ls='--')
    plt.xlabel('psiN')
    # plt.ylabel(P['Ne'])
    plt.title('Electron density with fits')
    # plt.title(plot_dic['ne3da.last10'][0],fontdict={"family":"Calibri","size": 20})
    plt.legend()
    
    plt.figure(figsize=(7,7))
    if log_flag:
        plt.yscale('log')
    else:
        pass
    plt.plot(x_coord, Te, ls='-', color = 'b', label= 'solps electron tempurature')
    plt.plot(x_coord, tanh_te_fit, ls='-', color='r',lw= 3, label= 'tanh fit')
    plt.plot(xt, yt, color='black', ls='-', lw=3, label= 'temperature pedestal width [m]: $\Delta n_e$')
    plt.axvline(x=dtn + te_sym_pt, color='black', lw=3, ls='-')
    plt.axvline(x=-dtn + te_sym_pt, color='black', lw=3, ls='-')
    plt.xlabel('psiN')
    # plt.ylabel(P['Te'])
    plt.title('Electron temperature with fits')
    plt.legend()
    
   

def opacity_radial_method_multi(result_dic, SEP, iter_list, 
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
    

    
    
    withshift = char['withshift']
    withseries = char['withseries']
    series_flag = char['series_flag']
    
    
    # ii = 1
    
    for i in iter_list:
        plt.figure(figsize=(7,7))
        xp = [-efold_dic[i] + max(xcoord_cut_dic[i]), max(xcoord_cut_dic[i])]
        yp = [min(exp_an_fit_dic[i]), min(exp_an_fit_dic[i])]
        plt.plot(psi_dic[i], Nd_dic[i], ls='-', color = 'green', label= 'solps neutral density')
        plt.plot(xcoord_cut_dic[i], exp_an_fit_dic[i], ls='-', color='r',lw= 5, label= 'exponential fit')
        plt.axvline(x=max(xcoord_cut_dic[i]), color='orange',lw=3)
        plt.plot(xp, yp, color='orange', lw=3, ls='-', label= 'Neutral penetration length: $\lambda_{n_D}$')
        plt.axvline(x = xp[0], color='orange', lw=3, ls='-')
        plt.axvline(x= delta_dic[i] + ne_sym_pt_dic[i], color='black',lw=3, ls='--', label= 'Pedestal width')
        plt.axvline(x=-delta_dic[i] + ne_sym_pt_dic[i], color='black',lw=3, ls='--')
        plt.xlabel('psiN')
        # plt.ylabel(P['NeuDen'])
        
        if withshift == True and withseries == False:
            shift_value = str(change_var_dic[i])
            plt.title('Modify {} m Neutral density with fits'.format(shift_value))
        elif withshift == False and withseries == True:
            if series_flag == 'change_den':
                shift_value = str(change_var_dic[i])
                plt.title('Neutral density with fits for electron core density {} m^3'.format(shift_value))
            elif series_flag == 'eireneN':
                shift_value = str(change_var_dic[i])
                plt.title('Neutral density with fits for SOL eirene particle as {}'.format(shift_value))
            else:
                print('There is a bug')
        plt.legend()
        # ii = ii + 1
        if log_flag:
            plt.yscale('log')
        else:
            pass
    
    for j in iter_list:
        plt.figure(figsize=(7,7))
        if withshift == True and withseries == False:
            sep = int(SEP[j])
        elif withshift == False and withseries == True:
            sep = int(SEP)
        xd = [-delta_dic[j] + ne_sym_pt_dic[j], delta_dic[j] + ne_sym_pt_dic[j]]
        yd = [tanh_ne_fit_dic[j][sep], tanh_ne_fit_dic[j][sep]]
        plt.plot(psi_dic[j], Ne_dic[j], ls='-', label= 'solps electron density_{}'.format(j))
        plt.plot(psi_dic[j], tanh_ne_fit_dic[j], ls='-', color='r',lw= 3, label= 'exponential fit')     
        plt.plot(xd, yd, ls='-', color='black', lw=3, label= 'Pedestal width [m]: $\Delta n_e$')
        plt.axvline(x= delta_dic[j] + ne_sym_pt_dic[j], color='black', lw=3, ls='-')
        plt.axvline(x=-delta_dic[j] + ne_sym_pt_dic[j], color='black', lw=3, ls='-')
        plt.axvline(x=max(xcoord_cut_dic[j]), color='orange',lw=3, ls='--')
        plt.axvline(x=-efold_dic[j] + max(xcoord_cut_dic[j]), color='orange',lw=3, ls='--')
        plt.xlabel('psiN')
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
        plt.plot(psi_dic[k], Te_dic[k], ls='-', color = 'b', label= 'solps electron tempurature')
        plt.plot(psi_dic[k], tanh_te_fit_dic[k], ls='-', color='r', lw= 3, label= 'tanh fit')
        plt.plot(xt, yt, color='black', lw=3, ls='-', label= 'temperature pedestal width [m]: $\Delta n_e$')
        plt.axvline(x=tdelta_dic[k] + te_sym_pt_dic[k], color='black', lw=3)
        plt.axvline(x=-tdelta_dic[k] + te_sym_pt_dic[k], color='black', lw=3)
        plt.xlabel('psiN')
        # plt.ylabel(P['Te'])
        plt.title('Electron temperature with fits')
        plt.legend()
    
    if withshift == False and withseries == True:
        density_dic = {}
        if series_flag == 'change_den':
            for k in iter_list:
                kk = float(k)*pow(10, 19)
                density_dic[k] = kk
            
            plt.figure(figsize=(7,7))
            for l in iter_list:
                plt.plot(psi_dic[l], Ne_dic[l], ls='-', label= 'solps electron density {}'.format(density_dic[l]))
                # plt.plot(psi_dic[j], tanh_ne_fit_dic[j], color='r',lw= 3, label= 'exponential fit')
            
            plt.xlabel('psiN')
            # plt.ylabel(P['Ne'])
            plt.title('Electron density with fits')
            plt.legend()
    


    
    



