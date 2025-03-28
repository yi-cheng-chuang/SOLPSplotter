# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 15:48:40 2025

@author: ychuang
"""

from SOLPSplotter_ndscan import neuden_scan
from twscan_target import twinscan_showflow
import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.offsetbox import AnchoredText
from scipy import interpolate


class core_edge(neuden_scan):
    
    def __init__(self, DefaultSettings, loadDS):
        neuden_scan.__init__(self, DefaultSettings, loadDS)
    
    
    
    
    def twscan_tarNTdata(self, iter_index, data_struc):
        
            
        if self.series_flag == 'twin_scan':
            
            
            nf = iter_index[0]
            tf = iter_index[1]
            
            # b2fstate = self.data['b2fstate'][nf][tf]
            
            # nx = data_struc['nx']
            # ny = data_struc['ny']
            
            if data_struc['size'] == 'full':
                ne_dat = self.data['outputdata']['Ne'][nf][tf]
                te_dat = self.data['outputdata']['Te'][nf][tf]
                psi_coord = self.data['psi']['psival']
                
                              
            elif data_struc['size'] == 'small':
                nx = data_struc['nx']
                ny = data_struc['ny']
                b2fstate = self.data['b2fstate'][nf][tf]
                ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
                Te_J = b2fstate['te'][1:nx+1, 1:ny+1]
                
                ev = 1.6021766339999999 * pow(10, -19)
                te_dat = Te_J / ev
                                
                psi_coord = self.data['psi']['psival'][1:ny+1, 1:nx+1]
                
                source = self.data['b2wdat'][nf][tf]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                vol = self.data['b2wdat'][nf][tf]['vol'][1:nx+1, 1:ny+1]
                sx = np.divide(source, vol)
                
                neuden_dat = self.data['ft44'][nf][tf]['dab2'][:, :, 0]
                # neuden_dat = np.transpose(data[:, :, 0])
            
        else:
            
            if data_struc['size'] == 'full':
                ne_dat = self.data['outputdata']['Ne'][iter_index]
                te_dat = self.data['outputdata']['Te'][iter_index]
                psi_coord = self.data['psi']['psival']
                
            elif data_struc['size'] == 'small':
                nx = data_struc['nx']
                ny = data_struc['ny']
                b2fstate = self.data['b2fstate'][iter_index]
                ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
                Te_J = b2fstate['te'][1:nx+1, 1:ny+1]
                
                ev = 1.6021766339999999 * pow(10, -19)
                te_dat = Te_J / ev
                
                source = self.data['b2wdat'][iter_index]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                vol = self.data['b2wdat'][iter_index]['vol'][1:nx+1, 1:ny+1]
                sx = np.divide(source, vol)
                
                neuden_dat = self.data['ft44'][iter_index]['dab2'][:, :, 0]
                # neuden_dat = np.transpose(data[:, :, 0])
                              
                psi_coord = self.data['psi']['psival'][1:ny+1, 1:nx+1]
            
            
        psi_dic = {'inner target': psi_coord[:, 0], 'outer target': psi_coord[:, nx-1]}
        ne_dic = {'inner target': ne_dat[0, :], 'outer target': ne_dat[nx-1, :]}
        te_dic = {'inner target': te_dat[0, :], 'outer target': te_dat[nx-1, :]}
        sx_dic = {'inner target': sx[0, :], 'outer target': sx[nx-1, :]}
        neuden_dic = {'inner target': neuden_dat[0, :], 'outer target': neuden_dat[nx-1, :]}
        
        
        
        return psi_dic, ne_dic, te_dic, sx_dic, neuden_dic
    
    
    
    def twscan_CE(self, dat_size, scan_var, format_option, pol_list):
        
        if format_option == 'detachment':
            fig, axs = plt.subplots(1, 2)
            
        elif format_option == 'neuden_peak':
            fig, axs = plt.subplots(1, 2)
        
        elif format_option == 'nd_sep':
            fig, axs = plt.subplots(1, 2)
        
        elif format_option == 'efold':
            fig, axs = plt.subplots(1, 2)
        
        
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        
        
        if dat_size == 'full':
    
            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        elif dat_size == 'small':
            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        plt.subplots_adjust(hspace=.0)
        anchored_text_1 = AnchoredText('{}'.format('Electron temperature at inner target [eV]'), loc='upper left')
        anchored_text_2 = AnchoredText('{}'.format('Electron temperature at outer target [eV]'), loc='upper left')
        # print('this is 201:')
        # print(dat_size)
        
        
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
                
                
                
                
                if '59' in pol_list:
                    
                    index = pol_list.index('59')
                    print('the index of OMP is {}'.format(index))
                
                else:
                    
                    print('OMP is not in pol list')
                
                
                
            
                for aa in iter_key:
                    
                    
                    psi_dic, ne_dic, te_dic, sx_dic, neuden_dic = self.twscan_tarNTdata(iter_index = aa, 
                                                            data_struc = dat_struc)
                    
                    ne_ped = self.data['opacity_poloidal'][aa]['electron_pedestal_density']
                    ne_sep = self.data['opacity_poloidal'][aa]['electron_density_separatrix']
                    efold_length = self.data['opacity_poloidal'][aa]['efold_length']*pow(10, 3)
                    
                    neuden_dat = self.data['ft44'][aa[0]][aa[1]]['dab2'][:, :, 0]
                    
                    st = int(pol_list[0])
                    ed = int(pol_list[-1]) + 1
                    
                    nd_sep = neuden_dat[st:ed, 18]
                    
                    
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
                    
                    
                    for ee in efold_length:
                        
                        multi_efold_list.append(ee)
                    
                    
                    for nep in ne_ped:
                        
                        multi_neped_list.append(nep)
                    
                    
                    for nes in ne_sep:
                        
                        multi_nesep_list.append(nes)
                    
                    for nds in nd_sep:
                        
                        multi_ndsep_list.append(nds)
                        
                    
                    
                        
                    
                if format_option == 'efold':
                    
                    axs[0].scatter(multi_nesep_list, multi_efold_list, color = color_dic[ap], label= '{:.3E} W'.format(label_ap))
                    axs[1].scatter(multi_neped_list, multi_efold_list, color = color_dic[ap], label= '{:.3E} W'.format(label_ap))           
                    axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                    axs[1].set_xlabel('$electron pedestal density at OMP$')
                    axs[0].add_artist(anchored_text_1)
                    axs[1].add_artist(anchored_text_2)
                    # axs[2].add_artist(anchored_text_3)
                    axs[0].set_title('efold length')
                    
                    axs[0].legend(loc= 'lower right')
                
                
                elif format_option == 'nd_sep':
                    
                    axs[0].scatter(multi_nesep_list, multi_ndsep_list, color = color_dic[ap], label= '{:.3E} W'.format(label_ap))
                    axs[1].scatter(multi_neped_list, multi_ndsep_list, color = color_dic[ap], label= '{:.3E} W'.format(label_ap))           
                    axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                    axs[1].set_xlabel('$electron pedestal density at OMP$')
                    axs[0].add_artist(anchored_text_1)
                    axs[1].add_artist(anchored_text_2)
                    # axs[2].add_artist(anchored_text_3)
                    axs[0].set_title('neutral density at the separatrix')
                    
                    axs[0].legend(loc= 'lower right')
                    
                    
                    
                
                print(color_dic)
                
                if format_option == 'detachment':
                    
                    if scan_var == 'ne_sep':
                        
                        axs[0].plot(ne_sep_list, te_inner_list,'-', color = color_dic[ap], label= 'HFS {:.3E} W'.format(label_ap))
                        axs[1].plot(ne_sep_list, te_outer_list,'-', color = color_dic[ap], label= 'LFS {:.3E} W'.format(label_ap))           
                        axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                        axs[1].set_xlabel('$electron density at the separatrix on OMP$')
                        axs[0].add_artist(anchored_text_1)
                        axs[1].add_artist(anchored_text_2)
                        # axs[2].add_artist(anchored_text_3)
                        axs[0].set_title('detachment cliff')
                        
                        axs[0].legend(loc= 'lower right')
                    
                    elif scan_var == 'ne_ped':
                        
                        axs[0].plot(ne_ped_list, te_inner_list,'-', color = color_dic[ap], label= 'HFS {:.3E} W'.format(label_ap))
                        axs[1].plot(ne_ped_list, te_outer_list,'-', color = color_dic[ap], label= 'LFS {:.3E} W'.format(label_ap))           
                        axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                        axs[1].set_xlabel('$electron density at the separatrix on OMP$')
                        axs[0].add_artist(anchored_text_1)
                        axs[1].add_artist(anchored_text_2)
                        # axs[2].add_artist(anchored_text_3)
                        axs[0].set_title('detachment cliff')
                        
                        axs[0].legend(loc= 'lower right')
                    
                    
                elif format_option == 'neuden_peak':
                    
                    
                    if scan_var == 'ne_sep':
                        
                        axs[0].plot(ne_sep_list, ndpeak_inner_list,'-', color = color_dic[ap], label= 'HFS {:.3E} W'.format(label_ap))
                        axs[1].plot(ne_sep_list, ndpeak_outer_list,'-', color = color_dic[ap], label= 'LFS {:.3E} W'.format(label_ap))           
                        axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                        axs[1].set_xlabel('$electron_density at the separatrix on OMP$')
                        axs[0].add_artist(anchored_text_1)
                        axs[1].add_artist(anchored_text_2)
                        # axs[2].add_artist(anchored_text_3)
                        axs[0].set_title('neuden peak')
                        
                        axs[0].legend(loc= 'lower right')
                    
                    elif scan_var == 'ne_ped':
                        
                        axs[0].plot(ne_ped_list, ndpeak_inner_list,'-', color = color_dic[ap], label= 'HFS {:.3E} W'.format(label_ap))
                        axs[1].plot(ne_ped_list, ndpeak_outer_list,'-', color = color_dic[ap], label= 'LFS {:.3E} W'.format(label_ap))           
                        axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                        axs[1].set_xlabel('$electron density at the separatrix on OMP$')
                        axs[0].add_artist(anchored_text_1)
                        axs[1].add_artist(anchored_text_2)
                        # axs[2].add_artist(anchored_text_3)
                        axs[0].set_title('neuden peak')
                        
                        axs[0].legend(loc= 'lower right')
                    
                    
                    
                        
                
                
                    
                    

                        
                        
                        
                        
    
                        
    
    
    
    
    
    
    def twscan_CEtest_spare(self, dat_size, scan_style, format_option):
        
        if format_option == 'detachment':
            fig, axs = plt.subplots(1, 2)
            
        elif format_option == 'neuden_peak':
            fig, axs = plt.subplots(1, 2)
        
        elif format_option == 'nd/ne_ped' or format_option == 'nd/ne_sep':
            fig, axs = plt.subplots()
        
        
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        
        
        if dat_size == 'full':
    
            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        elif dat_size == 'small':
            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        plt.subplots_adjust(hspace=.0)
        anchored_text_1 = AnchoredText('{}'.format('Electron temperature at inner target [eV]'), loc='upper left')
        anchored_text_2 = AnchoredText('{}'.format('Electron temperature at outer target [eV]'), loc='upper left')
        # print('this is 201:')
        # print(dat_size)
        
        
        if self.series_flag == 'twin_scan':
            
            dircomp = self.data['dircomp']
            
            if scan_style == 'tempscan':
                
                key_a = 'denscan_list'
                key_b = 'tempscan_list'
            
            elif scan_style == 'denscan':
                
                key_a = 'tempscan_list'
                key_b = 'denscan_list'
            
            else:
                print('twinscan_plot_method, please check the scan_style!')
            
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
                    
                    if scan_style == 'tempscan':
                        
                        it_in = (ta, tb)
                    
                    elif scan_style == 'denscan':
                        
                        it_in = (tb, ta)
                    
                    else:
                        print('twinscan_plot_method, please check the scan_style!')
                    
                    
                    iter_key.append(it_in)
                
                print('input key is:')
                print(iter_key)
                
                te_inner_list = []
                te_outer_list = []
                ne_ped_list = []
                ne_sep_list = []
                ndpeak_inner_list = []
                ndpeak_outer_list = []
            
                for aa in iter_key:
                    
                    
                    psi_dic, ne_dic, te_dic, sx_dic, neuden_dic = self.twscan_tarNTdata(iter_index = aa, 
                                                            data_struc = dat_struc)
                    
                    ne_ped = self.data['opacity_poloidal'][aa]['electron_pedestal_density']
                    ne_sep = self.data['opacity_poloidal'][aa]['electron_density_separatrix']
                    
                    ne_ped_list.append(ne_ped)
                    ne_sep_list.append(ne_sep)
                    
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
                    
                
                print(color_dic)
                
                if format_option == 'detachment':
                    
                    axs[0].plot(ne_sep_list, te_inner_list,'-', color = color_dic[ap], label= 'HFS {:.3E} W'.format(label_ap))
                    axs[1].plot(ne_sep_list, te_outer_list,'-', color = color_dic[ap], label= 'LFS {:.3E} W'.format(label_ap))           
                    axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                    axs[1].set_xlabel('$electron density at the separatrix on OMP$')
                    axs[0].add_artist(anchored_text_1)
                    axs[1].add_artist(anchored_text_2)
                    # axs[2].add_artist(anchored_text_3)
                    axs[0].set_title('detachment cliff')
                    
                    axs[0].legend(loc= 'lower right')
                
                
                elif format_option == 'neuden_peak':
                                        
                    axs[0].plot(ne_sep_list, ndpeak_inner_list,'-', color = color_dic[ap], label= 'HFS {:.3E} W'.format(label_ap))
                    axs[1].plot(ne_sep_list, ndpeak_outer_list,'-', color = color_dic[ap], label= 'LFS {:.3E} W'.format(label_ap))           
                    axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                    axs[1].set_xlabel('$electron density at the separatrix on OMP$')
                    axs[0].add_artist(anchored_text_1)
                    axs[1].add_artist(anchored_text_2)
                    # axs[2].add_artist(anchored_text_3)
                    axs[0].set_title('neuden peak')
                    
                    axs[0].legend(loc= 'lower right')

                                        
                    
        
        
        
    
    
                    
                
                
                
                
            
        
        
        
        
        
        
    
    
    
    
    
    
    



