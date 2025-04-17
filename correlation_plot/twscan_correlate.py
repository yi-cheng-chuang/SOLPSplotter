# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 19:49:55 2025

@author: ychuang
"""



import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.offsetbox import AnchoredText
import matplotlib.lines as mlines
from twscan_module.twinscan_prepare import twscan_assist
from correlation_plot.twscan_bin_plot import bin_plot


class twcorrelate:
    
    def __init__(self, DF, data, twa: twscan_assist, bp: bin_plot):
        
        self.DF = DF
        self.data = data
        self.twa = twa
        self.bp = bp
    
    
    def twscan_corr(self, format_option, pol_list, plot_case, scan_style):
        
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
        dat_struc = {'nx': nx, 'ny': ny}
        

        
        plt.subplots_adjust(hspace=.0)
        anchored_text_1 = AnchoredText('{}'.format('Electron temperature at inner target [eV]'), loc='upper left')
        anchored_text_2 = AnchoredText('{}'.format('Electron temperature at outer target [eV]'), loc='upper left')
        # print('this is 201:')
    
        
        
        if self.DF.series_flag == 'twin_scan':
            
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
                
                color_dic = self.twa.pair_dic(keys = keylist_a, values = color_list)
            
            else:
                pass
            
            for ta in keylist_a:
                
                keylist_b = []
                
                for x in dircomp[key_b]:
                    keylist_b.append('{:.3f}'.format(x))
                
                
                if plot_case == 'fivescan' or plot_case == 'single':
                    
                    color_list = ['red', 'orange', 'green', 'blue', 'purple']
                    
                    color_dic = self.twa.pair_dic(keys = keylist_b, values = color_list)
                
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
                opq_dat = {}
                
                
                
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
                    
                    nd_sep = neuden_dat[st:ed, 18]
                    
                    
                    source = self.data['b2wdat'][aa[0]][aa[1]]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                    vol = self.data['b2wdat'][aa[0]][aa[1]]['vol'][1:nx+1, 1:ny+1]
                    sx = np.divide(source, vol)
                    
                    psi_coord = self.data['psi']['psival'][1:ny+1, 1:nx+1]
                    

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
                        
                        for nds in nd_sep:
                            
                            multi_ndsep_list.append(nds)
                        
                        
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
                    
                    else:
                        pass
                        
                    
                    
                    
                    
                    k_list = self.bp.bin_with_angle()
                    
                    kp = len(k_list)
                    
                    # print(f'the length of k_list is: {kp}')
                    
                    
                    twneped_list = []
                    twnesep_list = []
                    twefold_list = []
                    twopq_list = []
                    twteped_list = []
                    twtesep_list = []
                    twteedge_list = []
                    twneedge_list = []
                    twspk_list = []
                    
                    
                    for em in range(7):
                        twneped_list.append([])
                        twnesep_list.append([])
                        twefold_list.append([])
                        twopq_list.append([])
                        twteped_list.append([])
                        twtesep_list.append([])
                        twteedge_list.append([])
                        twneedge_list.append([])
                        twspk_list.append([])
                        
                        
                    for i, kp in enumerate(k_list):
                        
                        for kc in kp:
                            
                            twneped_list[i].append(ne_ped[kc])
                            twnesep_list[i].append(ne_sep[kc])
                            twefold_list[i].append(efold_length[kc])
                            twopq_list[i].append(opq[kc])
                            twteped_list[i].append(te_ped[kc])
                            twtesep_list[i].append(te_sep[kc])
                            twteedge_list[i].append(te_edge[kc])
                            twneedge_list[i].append(ne_edge[kc])
                            twspk_list[i].append(s_peak_list[kc])
                            
                    
                    
                    
                    
                    if plot_case == 'fivescan':
                        
                        datasets[ad] = (twspk_list, color_dic[ad])
                        nesep_dat[ad] = (twnesep_list, color_dic[ad])
                        neped_dat[ad] = (twneped_list, color_dic[ad])
                        needge_dat[ad] = (twneedge_list, color_dic[ad])
                        tesep_dat[ad] = (twtesep_list, color_dic[ad])
                        teped_dat[ad] = (twteped_list, color_dic[ad])
                        teedge_dat[ad] = (twteedge_list, color_dic[ad])
                        opq_dat[ad] = (twopq_list, color_dic[ad])
                        
                    else:
                        pass
                            
                    
                    
                    
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
  
                        
                    if format_option == 'efold':
                        
                        
                        for i in range(len(k_list)):
                            
                            neped_ar = np.array(twneped_list[i])
                            nesep_ar = np.array(twnesep_list[i])
                            efold_ar = np.array(twefold_list[i])
                            
                            if i == 3:
                                axs[0].scatter(neped_ar, efold_ar, marker = marker_list[i], color = color_dic[ap])
                                axs[1].scatter(nesep_ar, efold_ar, marker = marker_list[i], color = color_dic[ap])
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
                        
                        axs[0].scatter(multi_nesep_list, multi_ndsep_list, color = color_dic[ap])
                        axs[1].scatter(multi_neped_list, multi_ndsep_list, color = color_dic[ap])           
                        axs[0].set_xlabel('$electron_density at the separatrix on OMP$')
                        axs[1].set_xlabel('$electron pedestal density at OMP$')
                        # axs[0].add_artist(anchored_text_1)
                        # axs[1].add_artist(anchored_text_2)
                        # axs[2].add_artist(anchored_text_3)
                        fig.suptitle('neutral density at the separatrix')
                        
                        axs[0].legend(loc= 'lower right')
                    
                    
                    elif format_option == 'source_peak':
                        
                        
                        if plot_case == 'single':
                            
                            for i in range(len(k_list)):
                                
                                neped_ar = np.array(twneped_list[i])
                                nesep_ar = np.array(twnesep_list[i])
                                spk_ar = np.array(twspk_list[i])
                            
                                
                                axs[0].scatter(nesep_ar, spk_ar, marker = marker_list[i], color = color_dic[ad], 
                                               label= marker_label[i])
                                axs[1].scatter(neped_ar, spk_ar, marker = marker_list[i], color = color_dic[ad], 
                                               label= marker_label[i])
                                
                                axs[0].legend(handles= marker_handle)
                                
                
                        elif plot_case == 'all25':
                            
                                neped_ar = np.array(twneped_list[i])
                                nesep_ar = np.array(twnesep_list[i])
                                spk_ar = np.array(twspk_list[i])
                                
                                
                                axs[0].scatter(nesep_ar, spk_ar, marker = marker_list[i], color = color_dic[ap], label= marker_label[i])
                                axs[1].scatter(neped_ar, spk_ar, marker = marker_list[i], color = color_dic[ap], label= marker_label[i])
                                
                                
                                axs[0].legend(handles= marker_handle)
                        
                        else:
                            pass
                                
                            
                        if plot_case == 'single' or plot_case == 'all25':
                            
                            axs[0].set_xlabel('$electron density at the separatrix$')
                            axs[1].set_xlabel('$electron pedestal density$')
                            if scan_style == 'temperature':
                                
                                fig.suptitle(f'source peak in psiN space {scan_style} scan from heat flux {st_val} to {ed_val} $10^5$ W')
                            
                            elif scan_style == 'density':
                                
                                fig.suptitle(f'source peak in psiN space {scan_style} scan from particle flux {st_val} to {ed_val} $10^20$ 1/s')
                        
                        else:
                            pass
                            
                
                    elif format_option == 'opaqueness':
                        
                        if plot_case == 'single':
                            
                            for i in range(len(k_list)):
                                
                                neped_ar = np.array(twneped_list[i])
                                nesep_ar = np.array(twnesep_list[i])
                                opq_ar = np.array(twopq_list[i])
                            
                                
                                axs[0].scatter(nesep_ar, opq_ar, marker = marker_list[i], color = color_dic[ad], 
                                               label= marker_label[i])
                                axs[1].scatter(neped_ar, opq_ar, marker = marker_list[i], color = color_dic[ad], 
                                               label= marker_label[i])
                                
                                axs[0].legend(handles= marker_handle)
                                
                
                        elif plot_case == 'all25':
                            
                                neped_ar = np.array(twneped_list[i])
                                nesep_ar = np.array(twnesep_list[i])
                                opq_ar = np.array(twopq_list[i])
                                
                                
                                axs[0].scatter(nesep_ar, opq_ar, marker = marker_list[i], color = color_dic[ap], label= marker_label[i])
                                axs[1].scatter(neped_ar, opq_ar, marker = marker_list[i], color = color_dic[ap], label= marker_label[i])
                                
                                
                                axs[0].legend(handles= marker_handle)
                        
                        else:
                            pass
                                
                            
                        if plot_case == 'single' or plot_case == 'all25':
                            
                            axs[0].set_xlabel('$n_{e, sep}$')
                            axs[1].set_xlabel('$n_{e, ped}$')
                            if scan_style == 'temperature':
                                
                                fig.suptitle(f'Neutral opaqueness {scan_style} scan from heat flux {st_val} to {ed_val} $10^5$ W')
                            
                            elif scan_style == 'density':
                                
                                fig.suptitle(f'Neutral opaqueness {scan_style} scan from particle flux {st_val} to {ed_val} $10^{{{20}}}$ 1/s')
                        
                        else:
                            pass
                    
                    
                    elif format_option == 'Te':
                        
                        for i in range(len(k_list)):
                            
                            neped_ar = np.array(twneped_list[i])
                            nesep_ar = np.array(twnesep_list[i])
                            teped_ar = np.array(twteped_list[i])
                            tesep_ar = np.array(twtesep_list[i])
                            teedge_ar = np.array(twteedge_list[i])
                            needge_ar = np.array(twneedge_list[i])
                            
                            # if i == 0:
                            #     axs[0].scatter(nesep_ar, opq_ar, marker = marker_list[i], color = color_dic[ap], label= '{:.3E} W'.format(label_ap))
                            #     axs[1].scatter(neped_ar, opq_ar, marker = marker_list[i], color = color_dic[ap], label= '{:.3E} W'.format(label_ap))
                            # else:
                            #     pass
                            if plot_case == 'single':
                                
                                axs[0].scatter(nesep_ar, tesep_ar, marker = marker_list[i], color = color_dic[ad], 
                                               label= marker_label[i])
                                axs[1].scatter(neped_ar, teped_ar, marker = marker_list[i], color = color_dic[ad],
                                               label= marker_label[i])
                                axs[2].scatter(needge_ar, teedge_ar, marker = marker_list[i], color = color_dic[ad],
                                               label= marker_label[i])
                                
                                axs[0].legend(handles= marker_handle)
                                
                                
                            
                            elif plot_case == 'all25':
                                
                                axs[0].scatter(nesep_ar, tesep_ar, marker = marker_list[i], color = color_dic[ap], label= marker_label[i])
                                axs[1].scatter(neped_ar, teped_ar, marker = marker_list[i], color = color_dic[ap], label= marker_label[i])
                                axs[2].scatter(needge_ar, teedge_ar, marker = marker_list[i], color = color_dic[ap], label= marker_label[i])
                                
                                axs[0].legend(handles= marker_handle)

                            
                        # Labels and title
                        axs[0].set_xlabel('$n_{e, sep}$')
                        axs[1].set_xlabel('$n_{e, ped}$')
                        axs[2].set_xlabel('$n_{e, edge}$')
                        axs[0].set_title('$T_{e, sep}$')
                        axs[1].set_title('$T_{e, ped}$')
                        axs[2].set_title('$T_{e, edge}$')
                        if scan_style == 'temperature':
                            
                            fig.suptitle(f'Electron temperature with heat flux {ad} $10^5$ W and particle flux {ap} $10^{{{20}}}$ 1/s')
                        
                        elif scan_style == 'density':
                            
                            fig.suptitle(f'Electron temperature with heat flux {ad} $10^5$ W and particle flux {ap} $10^{{{20}}}$ 1/s')
                        
                        
            
            
                if format_option == 'source_peak':
                    
                    
                    if plot_case == 'fivescan':
                        
                        for label, (data, color) in datasets.items():
                            for i in range(len(k_list)):  # Two sublists
                                
                                twnesep_dat = nesep_dat[label]
                                twneped_dat = neped_dat[label]
                                
                                # print(twnesep_dat)
                            
                                
                                neped_ar = np.array(twneped_dat[0][i])
                                nesep_ar = np.array(twnesep_dat[0][i])
                                ml = marker_label[i]
                                
                                axs[0].plot(nesep_ar, data[i], marker = marker_list[i], color=color, linestyle='')
                                axs[1].plot(neped_ar, data[i], marker = marker_list[i], color=color, linestyle='')

                        

                        # Add the two legends separately
                        legend1 = axs[0].legend(handles=color_handles, title='Color', loc='center left')
                        legend2 = axs[0].legend(handles=marker_handle, title='Marker', loc='center right')

                        # Add both legends to the plot
                        axs[0].add_artist(legend1)
                        # axs[0].add_artist(legend2)

                        # Labels and title
                        axs[0].set_xlabel('$n_{e, sep}$')
                        axs[1].set_xlabel('$n_{e, ped}$')
                        if scan_style == 'temperature':
                            
                            fig.suptitle(f'source peak in psiN space {scan_style} scan from heat flux {st_val} to {ed_val} $10^5$ W')
                        
                        elif scan_style == 'density':
                            
                            fig.suptitle(f'source peak in psiN space {scan_style} scan from particle flux {st_val} to {ed_val} $10^{{{20}}}$ 1/s')
                
                
                
                
                elif format_option == 'Te':
                    
                    
                    if plot_case == 'fivescan':
                        
                        for label, (data, color) in needge_dat.items():
                            for i in range(len(k_list)):  # Two sublists
                                
                                twnesep_dat = nesep_dat[label]
                                twneped_dat = neped_dat[label]
                                twtesep_dat = tesep_dat[label]
                                twteped_dat = teped_dat[label]
                                twteedge_dat = teedge_dat[label]
                                # print(twnesep_dat)
                            
                                
                                neped_ar = np.array(twneped_dat[0][i])
                                nesep_ar = np.array(twnesep_dat[0][i])
                                teped_ar = np.array(twteped_dat[0][i])
                                tesep_ar = np.array(twtesep_dat[0][i])
                                teedge_ar = np.array(twteedge_dat[0][i])
                                ml = marker_label[i]
                                
                                axs[0].plot(nesep_ar, tesep_ar, marker = marker_list[i], color=color, linestyle='')
                                axs[1].plot(neped_ar, teped_ar, marker = marker_list[i], color=color, linestyle='')
                                axs[2].plot(data[i], teedge_ar, marker = marker_list[i], color=color, linestyle='')

                        

                        # Add the two legends separately
                        legend1 = axs[0].legend(handles=color_handles, title='Color', loc='upper left')
                        legend2 = axs[1].legend(handles=marker_handle, title='Marker', loc='upper right')

                        # Add both legends to the plot
                        axs[0].add_artist(legend1)
                        # axs[0].add_artist(legend2)

                        # Labels and title
                        axs[0].set_xlabel('$n_{e, sep}$')
                        axs[1].set_xlabel('$n_{e, ped}$')
                        axs[2].set_xlabel('$n_{e, edge}$')
                        axs[0].set_title('$T_{e, sep}$')
                        axs[1].set_title('$T_{e, ped}$')
                        axs[2].set_title('$T_{e, edge}$')
                        if scan_style == 'temperature':
                            
                            fig.suptitle(f'Electron temperature {scan_style} scan from heat flux {st_val} to {ed_val} $10^5$ W')
                        
                        elif scan_style == 'density':
                            
                            fig.suptitle(f'Electron temperature {scan_style} scan from particle flux {st_val} to {ed_val} $10^{{{20}}}$ 1/s')
                    
                    
                    
                    
                
                elif format_option == 'opaqueness':
                    
                    
                    if plot_case == 'fivescan':
                        
                        for label, (data, color) in opq_dat.items():
                            for i in range(len(k_list)):  # Two sublists
                                
                                twnesep_dat = nesep_dat[label]
                                twneped_dat = neped_dat[label]

                                # print(twnesep_dat)
                            
                                
                                neped_ar = np.array(twneped_dat[0][i])
                                nesep_ar = np.array(twnesep_dat[0][i])
                                # opq_ar = np.array(twopq_dat[0][i])

                                ml = marker_label[i]
                                
                                axs[0].plot(nesep_ar, data[i], marker = marker_list[i], color=color, linestyle='')
                                axs[1].plot(neped_ar, data[i], marker = marker_list[i], color=color, linestyle='')

                        

                        # Add the two legends separately
                        legend1 = axs[0].legend(handles=color_handles, title='Color', loc='upper left')
                        legend2 = axs[1].legend(handles=marker_handle, title='Marker', loc='upper right')

                        # Add both legends to the plot
                        axs[0].add_artist(legend1)
                        # axs[0].add_artist(legend2)

                        # Labels and title
                        axs[0].set_xlabel('$n_{e, sep}$')
                        axs[1].set_xlabel('$n_{e, ped}$')
                        if scan_style == 'temperature':
                            
                            fig.suptitle(f'Neutral opaqueness {scan_style} scan from heat flux {st_val} to {ed_val} $10^5$ W')
                        
                        elif scan_style == 'density':
                            
                            fig.suptitle(f'Neutral opaqueness {scan_style} scan from particle flux {st_val} to {ed_val} $10^{{{20}}}$ 1/s')
                    
                    
                    
                        
                    
                
                    
                    
                    
                
               












            
            
            

                
                

                    
                    








