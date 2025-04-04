# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 10:22:01 2025

@author: ychuang
"""


from twscan_core_edge import core_edge
import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.offsetbox import AnchoredText
from scipy import interpolate
import matplotlib.lines as mlines


class twscan_polcompare(core_edge):
    
    def __init__(self, DefaultSettings, loadDS):
        core_edge.__init__(self, DefaultSettings, loadDS)
    
    
    def twscan_pol(self, dat_size, scan_var, format_option, pol_list, plot_case, scan_style, logscale):
        
        if plot_case == 'all25':
            
            opq_quant_list = ['efold', 'opaqueness', 'nd_sep']
            
            for aq in opq_quant_list:
                
                if format_option == aq:
                    
                    fig, axs = plt.subplots(1, 2)
                
                else:
                    pass
            
                
            if format_option == 'source_peak':
                fig, axs = plt.subplots(1, 2)
            
            elif format_option == 'Te':
                fig, axs = plt.subplots(1, 3)
        
        else:
            pass
        

    
        
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        
        
        if dat_size == 'full':

            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        elif dat_size == 'small':
            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        plt.subplots_adjust(hspace=.0)
        # print('this is 201:')
        # print(dat_size)
    
        
        
        if self.series_flag == 'twin_scan':
            
            dircomp = self.data['dircomp']
            
            if scan_style == 'temperature':
                
                key_a = 'denscan_list'
                key_b = 'tempscan_list'
            
            elif scan_style == 'density':
                
                key_a = 'tempscan_list'
                key_b = 'denscan_list'
            
            else:
                print('twinscan_plot_method, please check the scan_style!')
            
            keylist_a = []
            
            
            for x in dircomp[key_a]:
                keylist_a.append('{:.3f}'.format(x))
            
            
            if plot_case == 'all25':
                
                color_list = ['red', 'orange', 'green', 'blue', 'purple']
                
                color_dic = self.pair_dic(keys = keylist_a, values = color_list)
            
            else:
                pass
            
            for ta in keylist_a:
                
                keylist_b = []
                
                for x in dircomp[key_b]:
                    keylist_b.append('{:.3f}'.format(x))
                
                
                if plot_case == 'fivescan' or plot_case == 'single':
                    
                    color_list = ['red', 'orange', 'green', 'blue', 'purple']
                    
                    color_dic = self.pair_dic(keys = keylist_b, values = color_list)
                
                else:
                    pass
                    
                
                scan_list = []
                # print('scan_list after initial:')
                # print(scan_list)
                iter_key = []
                
                
                for tb in keylist_b:
                    
                    if scan_style == 'temperature':
                        
                        it_in = (ta, tb)
                    
                    elif scan_style == 'density':
                        
                        it_in = (tb, ta)
                    
                    else:
                        print('twinscan_plot_method, please check the scan_style!')
                    
                    
                    iter_key.append(it_in)
                    
                    
                    
                
                print('input key is:')
                print(iter_key)
                
                
                if plot_case == 'fivescan':
                    
                    opq_quant_list = ['efold', 'opaqueness', 'nd_sep']
                    
                    for aq in opq_quant_list:
                        
                        if format_option == aq:
                            
                            fig, axs = plt.subplots(1, 2)
                        
                        else:
                            pass
                    
                        
                    if format_option == 'source_peak':
                        fig, axs = plt.subplots(1, 2)
                    
                    elif format_option == 'Te':
                        fig, axs = plt.subplots(1, 3)
                    
                    
                
                else:
                    pass
                    
                
                

                ne_ped_list = []
                ne_sep_list = []
                nd_sep_list = []
                efold_list = []
                multi_efold_list = []
                multi_nesep_list = []
                multi_neped_list = []
                multi_ndsep_list = []
                multi_spk_list = []
                multi_opq_list = []
                datasets = {}
                nesep_dat = {}
                neped_dat = {}
                tesep_dat = {}
                teped_dat = {}
                teedge_dat = {}
                needge_dat = {}
                
                
                
                if '59' in pol_list:
                    
                    index = pol_list.index('59')
                    print('the index of OMP is {}'.format(index))
                
                else:
                    
                    print('OMP is not in pol list')
                
                
                for aa in iter_key:
                    
                    nf = aa[0]
                    tf = aa[1]
                    
                    ne_ped = self.data['opacity_poloidal'][aa]['electron_pedestal_density']
                    ne_sep = self.data['opacity_poloidal'][aa]['electron_density_separatrix']
                    te_ped = self.data['opacity_poloidal'][aa]['electron_pedestal_temperature']
                    te_sep = self.data['opacity_poloidal'][aa]['electron_temperature_separatrix']
                    efold_length = self.data['opacity_poloidal'][aa]['efold_length']*pow(10, 3)
                    
                    opq = self.data['opacity_poloidal'][aa]['dimensionless_opaqueness']
                    
                    neuden_dat = self.data['ft44'][aa[0]][aa[1]]['dab2'][:, :, 0]
                    b2fstate = self.data['b2fstate'][nf][tf]
                    ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
                    Te_J = b2fstate['te'][1:nx+1, 1:ny+1]
                    
                    ev = 1.6021766339999999 * pow(10, -19)
                    te_dat = Te_J / ev
                    
                    st = int(pol_list[0])
                    ed = int(pol_list[-1]) + 1
                    
                    nd_plot = neuden_dat[st:ed, :]
                    
                    
                    source = self.data['b2wdat'][aa[0]][aa[1]]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                    vol = self.data['b2wdat'][aa[0]][aa[1]]['vol'][1:nx+1, 1:ny+1]
                    sx = np.divide(source, vol)
                    
                    sx_plot = sx[st:ed, :]
                    
                    psi_coord = self.data['psi']['psival'][1:ny+1, 1:nx+1]
                    
                    
                    psi_plot = psi_coord[:, st:ed]
                    

                    int_pol_list = [int(x) for x in pol_list]
                    
                    s_peak_list = []
                    te_edge = []
                    ne_edge = []
                    
                    for pol in int_pol_list:
                        
                        sx_list = list(sx[pol, :])
                        psi_list = psi_coord[:, pol]
                        
                        
                        sx_ind = sx_list.index(max(sx_list))
                        s_peak_list.append(psi_list[sx_ind])
                        
                        te_edge.append(te_dat[pol, -1])
                        ne_edge.append(ne_dat[pol, -1])
                        
                    
                    if '59' in pol_list:
                        
                        ne_ped_list.append(ne_ped[index])
                        ne_sep_list.append(ne_sep[index])
                    
                    else:
                        
                        print('OMP is not in pol list')
                        
                    

                    
                    if scan_style == 'temperature':
                        
                        ad = aa[1]
                        ap = aa[0]
                    
                    elif scan_style == 'density':
                        
                        ad = aa[0]
                        ap = aa[1]
                    
                    else:
                        print('neteTSplot_method, please check scan_style')
                    
                    
                    
                    if scan_style == 'temperature':
                        
                        if plot_case == 'all25':
                            
                            st_val = float(iter_key[0][0])
                            ed_val = float(iter_key[-1][0])
                            
                        
                        if plot_case == 'single' or plot_case == 'fivescan':
                            
                            
                            st_val = float(iter_key[0][1])
                            ed_val = float(iter_key[-1][1])
                        
                        
                    
                    elif scan_style == 'density':
                        
                        if plot_case == 'all25':
                            
                            st_val = float(iter_key[0][1])
                            ed_val = float(iter_key[-1][1])
                            
                        
                        if plot_case == 'single' or plot_case == 'fivescan':
                            
                            
                            st_val = float(iter_key[0][0])
                            ed_val = float(iter_key[-1][0])

                    
                    
                    
                    if plot_case == 'fivescan':
                        
                        for ee in efold_length:
                            
                            multi_efold_list.append(ee)
                                                
                        for nep in ne_ped:
                            
                            multi_neped_list.append(nep)
                                                
                        for nes in ne_sep:
                            
                            multi_nesep_list.append(nes)
                        
                        
                        for spk in s_peak_list:
                            
                            multi_spk_list.append(spk)
                        
                        for op in opq:
                            
                            multi_opq_list.append(op)
                    
                    else:
                        pass
                    
                    

                    if plot_case == 'single':
                        
                        opq_quant_list = ['efold', 'opaqueness', 'nd_sep']
                        
                        for aq in opq_quant_list:
                            
                            if format_option == aq:
                                
                                fig, axs = plt.subplots(1, 2)
                            
                            else:
                                pass
                        
                            
                        if format_option == 'source_peak':
                            fig, axs = plt.subplots(1, 2)
                        
                        elif format_option == 'Te':
                            fig, axs = plt.subplots(1, 3)
                        
                        elif format_option == '':
                            fig, axs = plt.subplots()
                        
                        elif format_option == 'neuden':
                            fig, axs = plt.subplots()
                    
                    else:
                        pass
                        
                    
                    
                    
                    
                    k_list = self.bin_with_angle()
                    
                    kp = len(k_list)
                    
                    # print(f'the length of k_list is: {kp}')
                    
                    
                   
                    
                    

                        
                    psiplot = []
                    ndplot = []
                    
                    
                    for em in range(7):
                        psiplot.append([])
                        ndplot.append([])
                        
                        
                    for i, kp in enumerate(k_list):
                        
                        for kc in kp:
                            
                            psiplot[i].append(psi_plot[:, kc])
                            ndplot[i].append(nd_plot[kc, :])
                    
                    
                    # Sum using a loop
                    
                    
                    psi_total = []
                    nd_total = []
                    
                    for pii in range(7):
                        
                        total_topsi = np.zeros_like(psiplot[0][0])
                        total_tond = np.zeros_like(ndplot[0][0])
                                             
                        for arr in psiplot[pii]:
                            total_topsi += arr
                        
                        average_psi = total_topsi / len(psiplot[pii])
                        psi_total.append(average_psi)
                    
                        for arnd in ndplot[pii]:
                            total_tond += arnd
                        
                        average_nd = total_tond / len(ndplot[pii])
                        nd_total.append(average_nd)
                        
                        
                    
                    # average = total / len(arrays)
                    
                    # print("Sum:", total)
                    # print("Average:", average)
                            
                            
                    
                    
                    
                    # print(ne_ped)
                    
                    marker_list = ["o", "v", "s", "p", "x", "D", "1"]
                    
                    
                    angsep_list = np.linspace(250, -100, 8)
                    marker_label = []
                    for ai in range(len(angsep_list) -1):
                        
                        aup = angsep_list[ai]
                        adown = angsep_list[ai + 1]
                        marker_label.append(f'{aup}~{adown}')
                    
                    # print(marker_label)
                        
                    
                    # Create a custom legend handle
                    
                    marker_handle = []
                    
                    for imk, kk in enumerate(marker_list):
                        
                        if plot_case == 'single':
                            
                            marker_set = mlines.Line2D([], [], color= 'black', marker= kk,
                                                linestyle='None', markersize=10, label= marker_label[imk])
                        
                        elif plot_case == 'fivescan':
                            
                            marker_set = mlines.Line2D([0], [0], color= 'black', marker= kk,
                                                linestyle='None', markersize=10, label= marker_label[imk])
                        
                        elif plot_case == 'all25':
                            
                            marker_set = mlines.Line2D([], [], color= color_dic[ap], marker= kk,
                                                linestyle='None', markersize=10, label= marker_label[imk])
                        
                        
                        
                        marker_handle.append(marker_set)
                    
                    
                    # Create dummy handles for color legend (dataset)
                    color_handles = []
                    
                    for ic, cl in enumerate(color_list):
                        
                        if scan_style == 'density':
                            
                            val = keylist_b[ic]
                            # print(val)
                            
                            color_set = mlines.Line2D([0], [0], color= cl, marker='o', label=f'{val} $10^{{{20}}}$ 1/s')
                        
                        elif scan_style == 'temperature':
                            
                            val = keylist_b[ic]
                            
                            
                            color_set = mlines.Line2D([0], [0], color= cl, marker='o', label=f'{val} $10^5$ W')
                        
                        color_handles.append(color_set)

            
                    if format_option == 'neuden':
                        
                        
                        if plot_case == 'single':
                            
                            color_set = ['red', 'orange', 'green', 'blue', 'purple', 'brown', 'black']
                            
                            if logscale:
                                axs.set_yscale('log')
                            
                            else:
                                pass
                            
                            for i in range(len(k_list)):  # Two sublists
                                
                                
                                psi_ar = np.array(psi_total[i])
                                nd_ar = np.array(nd_total[i])
                                ml = marker_label[i]
                                
                                axs.plot(psi_ar, nd_ar, color= color_set[i], label = marker_label[i],
                                         linestyle='-')
                                
                            
    
                            # Labels and title
                            axs.set_xlabel('$\psi_N$')
                            axs.legend(loc = 'upper left')
                            if scan_style == 'temperature':
                                
                                fig.suptitle(f'Neutral density versus $\psi_N$')
                            
                            elif scan_style == 'density':
                                
                                fig.suptitle(f'Neutral density versus $\psi_N$')
                
                
                
                
                
                

                    
                    








