# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 15:48:40 2025

@author: ychuang
"""


import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.offsetbox import AnchoredText
from twscan_module.twinscan_prepare import twscan_assist
from correlation_plot.twscan_bin_plot import bin_plot


class core_edge:
    
    def __init__(self, DF, data, twa: twscan_assist, bp: bin_plot):
        
        self.DF = DF
        self.data = data
        self.twa = twa
        self.bp = bp


    
    
    def twscan_CE(self, scan_var, format_option, pol_list, plot_case):
        
        if plot_case == 'all25':
            
            if format_option == 'detachment':
                fig, axs = plt.subplots(1, 2)
                
            elif format_option == 'neuden_peak':
                fig, axs = plt.subplots(1, 2)
            
            elif format_option == 'nd_sep':
                fig, axs = plt.subplots(1, 2)
            
            elif format_option == 'efold':
                fig, axs = plt.subplots(1, 2)
            
            elif format_option == 'source_peak':
                fig, axs = plt.subplots(1, 2)
            
            elif format_option == 'opaqueness':
                fig, axs = plt.subplots(1, 2)
            
            elif format_option == 'Te':
                fig, axs = plt.subplots(1, 2)
        
        else:
            pass
            
            
        
        
        
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        dat_struc = {'nx': nx, 'ny': ny}
            
        
        plt.subplots_adjust(hspace=.0)
        anchored_text_1 = AnchoredText('{}'.format('Electron temperature at inner target [eV]'), loc='upper left')
        anchored_text_2 = AnchoredText('{}'.format('Electron temperature at outer target [eV]'), loc='upper left')
        # print('this is 201:')
        
        
        if self.series_flag == 'twin_scan':
            
            dircomp = self.data['dircomp']
            
            scan_style = 'denscan'
                
            key_a = 'tempscan_list'
            key_b = 'denscan_list'
            
            keylist_a = []
            
            
            for x in dircomp[key_a]:
                keylist_a.append('{:.3f}'.format(x))
            
            for ta in keylist_a:
                
                keylist_b = []
                
                for x in dircomp[key_b]:
                    keylist_b.append('{:.3f}'.format(x))
                
                color_list = ['red', 'orange', 'green', 'blue', 'purple']
                
                color_dic = self.pair_dic(keys = keylist_a, values = color_list)
                
                scan_list = []
                # print('scan_list after initial:')
                # print(scan_list)
                iter_key = []
                
                
                for tb in keylist_b:
                    
                    it_in = (tb, ta)
                    
                    iter_key.append(it_in)
                
                print('input key is:')
                print(iter_key)
                
                
                
                
                
                if plot_case == 'fivescan':
                    
                    if format_option == 'detachment':
                        fig, axs = plt.subplots(1, 2)
                        
                    elif format_option == 'neuden_peak':
                        fig, axs = plt.subplots(1, 2)
                    
                    elif format_option == 'nd_sep':
                        fig, axs = plt.subplots(1, 2)
                    
                    elif format_option == 'efold':
                        fig, axs = plt.subplots(1, 2)
                    
                    elif format_option == 'source_peak':
                        fig, axs = plt.subplots(1, 2)
                    
                    elif format_option == 'opaqueness':
                        fig, axs = plt.subplots(1, 2)
                    
                    elif format_option == 'Te':
                        fig, axs = plt.subplots(1, 2)
                
                else:
                    pass
                    
                
                
                
                
                te_inner_list = []
                te_outer_list = []
                ne_ped_list = []
                ne_sep_list = []
                ndpeak_inner_list = []
                ndpeak_outer_list = []
                nd_sep_list = []
                efold_list = []
                multi_efold_list = []
                multi_nesep_list = []
                multi_neped_list = []
                multi_ndsep_list = []
                multi_spk_list = []
                multi_opq_list = []
                
                
                
                if '59' in pol_list:
                    
                    index = pol_list.index('59')
                    print('the index of OMP is {}'.format(index))
                
                else:
                    
                    print('OMP is not in pol list')
                
                
                for aa in iter_key:
                    
                    
                    psi_dic = self.data['target_profile'][aa[0]][aa[1]]['psiN']
                    ne_dic = self.data['target_profile'][aa[0]][aa[1]]['ne']
                    te_dic = self.data['target_profile'][aa[0]][aa[1]]['te']
                    sx_dic = self.data['target_profile'][aa[0]][aa[1]]['source']
                    neuden_dic = self.data['target_profile'][aa[0]][aa[1]]['neuden']
                    
                    ne_ped = self.data['opacity_poloidal'][aa]['electron_pedestal_density']
                    ne_sep = self.data['opacity_poloidal'][aa]['electron_density_separatrix']
                    te_ped = self.data['opacity_poloidal'][aa]['electron_pedestal_temperature']
                    te_sep = self.data['opacity_poloidal'][aa]['electron_temperature_separatrix']
                    efold_length = self.data['opacity_poloidal'][aa]['efold_length']*pow(10, 3)
                    
                    opq = self.data['opacity_poloidal'][aa]['dimensionless_opaqueness']
                    
                    neuden_dat = self.data['ft44'][aa[0]][aa[1]]['dab2'][:, :, 0]
                    
                    st = int(pol_list[0])
                    ed = int(pol_list[-1]) + 1
                    
                    nd_sep = neuden_dat[st:ed, 18]
                    
                    
                    source = self.data['b2wdat'][aa[0]][aa[1]]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                    vol = self.data['b2wdat'][aa[0]][aa[1]]['vol'][1:nx+1, 1:ny+1]
                    sx = np.divide(source, vol)
                    
                    psi_coord = self.data['psi']['psival'][1:ny+1, 1:nx+1]
                    

                    int_pol_list = [int(x) for x in pol_list]
                    
                    s_peak_list = []
                    
                    for pol in int_pol_list:
                        
                        sx_list = list(sx[pol, :])
                        psi_list = psi_coord[:, pol]
                        
                        
                        sx_ind = sx_list.index(max(sx_list))
                        s_peak_list.append(psi_list[sx_ind])
                        
                        
                    
                    if '59' in pol_list:
                        
                        ne_ped_list.append(ne_ped[index])
                        ne_sep_list.append(ne_sep[index])
                    
                    else:
                        
                        print('OMP is not in pol list')
                        
                    
                    
                    if format_option == 'detachment':
                        
                        te_inner = te_dic['inner target']
                        te_outer = te_dic['outer target']
                        
                        te_inner_list.append(max(te_inner))
                        te_outer_list.append(max(te_outer))
                    
                    elif format_option == 'neuden_peak':
                        
                        nd_inner = neuden_dic['inner target']
                        nd_outer = neuden_dic['outer target']
                        
                        ndpeak_inner_list.append(max(nd_inner))
                        ndpeak_outer_list.append(max(nd_outer))
                    
                    elif format_option == 'efold':
                        
                        efold_list.append(efold_length)
                                          
                    elif format_option == 'nd_sep':
                                              
                        nd_sep_list.append(nd_sep)
                        
                    
                    
                    if self.series_flag == 'twin_scan':
                        
                        if scan_style == 'tempscan':
                            
                            ad = aa[1]
                            ap = aa[0]
                        
                        elif scan_style == 'denscan':
                            
                            ad = aa[0]
                            ap = aa[1]
                        
                        else:
                            print('neteTSplot_method, please check scan_style')
                    
                    else:
                        ad = aa
                    
                    label_ap = float(ap)*pow(10, 5)
                    
                    if plot_case == 'fivescan':
                        
                        for ee in efold_length:
                            
                            multi_efold_list.append(ee)
                        
                        
                        for nep in ne_ped:
                            
                            multi_neped_list.append(nep)
                        
                        
                        for nes in ne_sep:
                            
                            multi_nesep_list.append(nes)
                        
                        for nds in nd_sep:
                            
                            multi_ndsep_list.append(nds)
                        
                        
                        for spk in s_peak_list:
                            
                            multi_spk_list.append(spk)
                        
                        for op in opq:
                            
                            multi_opq_list.append(op)
                    
                    else:
                        pass
                        
                            
                    
                        
                    
                    
                    
                    if plot_case == 'single':
                        
                        if format_option == 'detachment':
                            fig, axs = plt.subplots(1, 2)
                            
                        elif format_option == 'neuden_peak':
                            fig, axs = plt.subplots(1, 2)
                        
                        elif format_option == 'nd_sep':
                            fig, axs = plt.subplots(1, 2)
                        
                        elif format_option == 'efold':
                            fig, axs = plt.subplots(1, 2)
                        
                        elif format_option == 'source_peak':
                            fig, axs = plt.subplots(1, 2)
                        
                        elif format_option == 'opaqueness':
                            fig, axs = plt.subplots(1, 2)
                        
                        elif format_option == 'Te':
                            fig, axs = plt.subplots(1, 2)
                    
                    else:
                        pass
                        
                    
                    
                    
                    
                    k_list = self.bin_with_angle()
                    
                    kp = len(k_list)
                    
                    # print(f'the length of k_list is: {kp}')
                    
                    
                    twneped_list = []
                    twnesep_list = []
                    twefold_list = []
                    twopq_list = []
                    twteped_list = []
                    twtesep_list = []
                    
                    
                    for em in range(7):
                        twneped_list.append([])
                        twnesep_list.append([])
                        twefold_list.append([])
                        twopq_list.append([])
                        twteped_list.append([])
                        twtesep_list.append([])
                        
                        
                    for i, kp in enumerate(k_list):
                        
                        for kc in kp:
                            
                            twneped_list[i].append(ne_ped[kc])
                            twnesep_list[i].append(ne_sep[kc])
                            twefold_list[i].append(efold_length[kc])
                            twopq_list[i].append(opq[kc])
                            twteped_list[i].append(te_ped[kc])
                            twtesep_list[i].append(te_sep[kc])
                            
                            
                            
                    
                    # print(ne_ped)
                    
                    marker_list = ["o", "v", "s", "p", "x", "D", "1"]
                    
                    
                    
                    
                        
                    if format_option == 'efold':
                        
                        
                        for i in range(len(k_list)):
                            
                            neped_ar = np.array(twneped_list[i])
                            nesep_ar = np.array(twnesep_list[i])
                            efold_ar = np.array(twefold_list[i])
                            
                            if i == 3:
                                axs[0].scatter(neped_ar, efold_ar, marker = marker_list[i], color = color_dic[ap], label= '{:.3E} W'.format(label_ap))
                                axs[1].scatter(nesep_ar, efold_ar, marker = marker_list[i], color = color_dic[ap], label= '{:.3E} W'.format(label_ap))
                            else:
                                pass
                            
                            # axs[0].scatter(neped_ar, efold_ar, marker = marker_list[i])
                            # axs[1].scatter(nesep_ar, efold_ar, marker = marker_list[i])
                            
                        
                        # print(neped_ar) 
                        # print(efold_ar)
                        # axs[0].scatter(neped_ar, efold_ar)
                        # axs[1].scatter(nesep_ar, efold_ar)
                            
                        axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                        axs[1].set_xlabel('$electron pedestal density at OMP$')
                        fig.suptitle('efold length')

                        
                        
                    
                    
                    elif format_option == 'nd_sep':
                        
                        axs[0].scatter(multi_nesep_list, multi_ndsep_list, color = color_dic[ap], label= '{:.3E} W'.format(label_ap))
                        axs[1].scatter(multi_neped_list, multi_ndsep_list, color = color_dic[ap], label= '{:.3E} W'.format(label_ap))           
                        axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                        axs[1].set_xlabel('$electron pedestal density at OMP$')
                        # axs[0].add_artist(anchored_text_1)
                        # axs[1].add_artist(anchored_text_2)
                        # axs[2].add_artist(anchored_text_3)
                        fig.suptitle('neutral density at the separatrix')
                        
                        axs[0].legend(loc= 'lower right')
                    
                    
                    elif format_option == 'source_peak':
                        
                        axs[0].scatter(multi_nesep_list, multi_spk_list, color = color_dic[ap], label= '{:.3E} W'.format(label_ap))
                        axs[1].scatter(multi_neped_list, multi_spk_list, color = color_dic[ap], label= '{:.3E} W'.format(label_ap))           
                        axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                        axs[1].set_xlabel('$electron pedestal density at OMP$')
                        fig.suptitle('source peak location in psiN')
                        
                        axs[0].legend(loc= 'center right')
                    
                    
                    elif format_option == 'opaqueness':
                        
                        for i in range(len(k_list)):
                            
                            neped_ar = np.array(twneped_list[i])
                            nesep_ar = np.array(twnesep_list[i])
                            opq_ar = np.array(twopq_list[i])
                            
                            # if i == 0:
                            #     axs[0].scatter(nesep_ar, opq_ar, marker = marker_list[i], color = color_dic[ap], label= '{:.3E} W'.format(label_ap))
                            #     axs[1].scatter(neped_ar, opq_ar, marker = marker_list[i], color = color_dic[ap], label= '{:.3E} W'.format(label_ap))
                            # else:
                            #     pass
                            
                            axs[0].scatter(nesep_ar, opq_ar, marker = marker_list[i])
                            axs[1].scatter(neped_ar, opq_ar, marker = marker_list[i])
                            

                            
                        axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                        axs[1].set_xlabel('$electron pedestal density at OMP$')
                        fig.suptitle('opaqueness')
                        
                        # axs[0].legend(loc= 'center right')
                    
                    
                    elif format_option == 'Te':
                        
                        for i in range(len(k_list)):
                            
                            neped_ar = np.array(twneped_list[i])
                            nesep_ar = np.array(twnesep_list[i])
                            teped_ar = np.array(twteped_list[i])
                            tesep_ar = np.array(twtesep_list[i])
                            
                            # if i == 0:
                            #     axs[0].scatter(nesep_ar, opq_ar, marker = marker_list[i], color = color_dic[ap], label= '{:.3E} W'.format(label_ap))
                            #     axs[1].scatter(neped_ar, opq_ar, marker = marker_list[i], color = color_dic[ap], label= '{:.3E} W'.format(label_ap))
                            # else:
                            #     pass
                            if plot_case == 'single':
                                
                                axs[0].scatter(nesep_ar, tesep_ar, marker = marker_list[i])
                                axs[1].scatter(neped_ar, teped_ar, marker = marker_list[i])
                            
                            elif plot_case == 'fivescan' or plot_case == 'all25':
                                
                                axs[0].scatter(nesep_ar, tesep_ar, marker = marker_list[i], color = color_dic[ap], label= '{:.3E} W'.format(label_ap))
                                axs[1].scatter(neped_ar, teped_ar, marker = marker_list[i], color = color_dic[ap], label= '{:.3E} W'.format(label_ap))

                            
                        axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                        axs[1].set_xlabel('$electron pedestal density at OMP$')
                        axs[0].set_title('$electron temperature at the separatrix$')
                        axs[1].set_title('$electron pedestal temperature at OMP$')
                        
                        # fig.suptitle('Te')
              
                # print(color_dic)
                
                if format_option == 'detachment':
                    
                    if scan_var == 'ne_sep':
                        
                        axs[0].plot(ne_sep_list, te_inner_list,'-', color = color_dic[ap], label= 'HFS {:.3E} W'.format(label_ap))
                        axs[1].plot(ne_sep_list, te_outer_list,'-', color = color_dic[ap], label= 'LFS {:.3E} W'.format(label_ap))           
                        axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                        axs[1].set_xlabel('$electron density at the separatrix on OMP$')
                        # axs[0].add_artist(anchored_text_1)
                        # axs[1].add_artist(anchored_text_2)
                        # axs[2].add_artist(anchored_text_3)
                        fig.suptitle('detachment cliff')
                        
                        axs[0].legend(loc= 'lower right')
                    
                    elif scan_var == 'ne_ped':
                        
                        axs[0].plot(ne_ped_list, te_inner_list,'-', color = color_dic[ap], label= 'HFS {:.3E} W'.format(label_ap))
                        axs[1].plot(ne_ped_list, te_outer_list,'-', color = color_dic[ap], label= 'LFS {:.3E} W'.format(label_ap))           
                        axs[0].set_xlabel('$electron pedestal density at OMP$')
                        axs[1].set_xlabel('$electron pedestal density at OMP$')
                        axs[0].add_artist(anchored_text_1)
                        axs[1].add_artist(anchored_text_2)
                        # axs[2].add_artist(anchored_text_3)
                        fig.suptitle('detachment cliff')
                        
                        axs[0].legend(loc= 'lower right')
                    
                    
                elif format_option == 'neuden_peak':
                    
                    
                    if scan_var == 'ne_sep':
                        
                        axs[0].plot(ne_sep_list, ndpeak_inner_list,'-', color = color_dic[ap], label= 'HFS {:.3E} W'.format(label_ap))
                        axs[1].plot(ne_sep_list, ndpeak_outer_list,'-', color = color_dic[ap], label= 'LFS {:.3E} W'.format(label_ap))           
                        axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                        axs[1].set_xlabel('$electron_density at the separatrix on OMP$')
                        # axs[0].add_artist(anchored_text_1)
                        # axs[1].add_artist(anchored_text_2)
                        # axs[2].add_artist(anchored_text_3)
                        fig.suptitle('neuden peak')
                        
                        axs[0].legend(loc= 'lower right')
                    
                    elif scan_var == 'ne_ped':
                        
                        axs[0].plot(ne_ped_list, ndpeak_inner_list,'-', color = color_dic[ap], label= 'HFS {:.3E} W'.format(label_ap))
                        axs[1].plot(ne_ped_list, ndpeak_outer_list,'-', color = color_dic[ap], label= 'LFS {:.3E} W'.format(label_ap))           
                        axs[0].set_xlabel('$electron pedestal density at OMP$')
                        axs[1].set_xlabel('$electron pedestal density at OMP$')
                        axs[0].add_artist(anchored_text_1)
                        axs[1].add_artist(anchored_text_2)
                        # axs[2].add_artist(anchored_text_3)
                        fig.suptitle('neuden peak')
                        
                        axs[0].legend(loc= 'lower right')
                
                
                

                    
                    

                        
                        
    
    
    
    
    
    
    
    
    