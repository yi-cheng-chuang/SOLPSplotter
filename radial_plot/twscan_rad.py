# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 06:24:13 2025

@author: ychuang
"""




from twscan_module.twinscan_prepare import twscan_assist
import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.offsetbox import AnchoredText
from scipy import interpolate
import matplotlib.lines as mlines


class twscan_radial:
    
    
    def __init__(self, DF, data, twa: twscan_assist):
        
        self.DF = DF
        self.data = data
        self.twa = twa
    
    
    
    def scan_rad(self, format_option, pol_loc, plot_case, scan_style, log_scale):
        
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
                        fig, axs = plt.subplots()
                    
                    elif format_option == 'neuden':
                        fig, axs = plt.subplots()
                    
                    elif format_option == 'Te':
                        fig, axs = plt.subplots(1, 3)
                
                else:
                    pass


                multi_spk_list = []
                multi_opq_list = []
                datasets = {}
                nesep_dat = {}
                neped_dat = {}
               
                
                for aa in iter_key:
                    
                    nf = aa[0]
                    tf = aa[1]
                    
                    
                    neuden_dat = self.data['ft44'][aa[0]][aa[1]]['dab2'][:, :, 0]
                    b2fstate = self.data['b2fstate'][nf][tf]
                    ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
                    Te_J = b2fstate['te'][1:nx+1, 1:ny+1]
                    
                    ev = 1.6021766339999999 * pow(10, -19)
                    te_dat = Te_J / ev
                    
                    
                    source = self.data['b2wdat'][aa[0]][aa[1]]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                    vol = self.data['b2wdat'][aa[0]][aa[1]]['vol'][1:nx+1, 1:ny+1]
                    sx = np.divide(source, vol)
                    
                    psi_coord = self.data['psi']['psival'][1:ny+1, 1:nx+1]
                    
                    pol = int(pol_loc)
                    
                    
                    sx_list = sx[pol, :]
                    psi_list = psi_coord[:, pol]
                    neuden_list = neuden_dat[pol, :]
                    
                    
                    psi_cut = []
                    S_cut = []
                    nu_cut = []
                    
                    
                    for ind, coord in enumerate(psi_list):
                        
                        if coord >= 0.95 and coord <= 1.1:
                            psi_cut.append(coord)
                            S_cut.append(sx_list[ind])
                            nu_cut.append(neuden_list[ind])
                            
                    
                    S_plot = S_cut/max(S_cut)
                    
                    
                    # psi_nu = []
                    # nu_cut = []
                    
                    # for ind, coord in enumerate(psi_list):
                        
                    #     if coord >= 0.95 and coord <= 1.05:
                    #         psi_nu.append(coord)
                    #         nu_cut.append(neuden_list[ind])
                    
                            
                    target = 1
                    
                    # Compute the closest value and its index
                    closest_value = min(psi_cut, key=lambda x: abs(x - target))
                    closest_index = psi_cut.index(closest_value)
                    
                    
                    if closest_value > 1:
                        
                        w = (target - psi_cut[closest_index -1])/(closest_value - psi_cut[closest_index -1])
                        nu_sep = w*nu_cut[closest_index] + (1- w)*nu_cut[closest_index -1]
                    
                    elif closest_value < 1:
                        
                        w = (target - closest_value)/(psi_cut[closest_index + 1] -closest_value)
                        nu_sep = w*nu_cut[closest_index + 1] + (1- w)*nu_cut[closest_index]
                    
                    else:
                        print('check scan_rad function in twscan_rad.py')
                        
                    
                    nu_plot = nu_cut/nu_sep

                    
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
                    
                    
                    if format_option == 'source_peak':
                        
                        
                        if plot_case == 'fivescan':
                            
                            
                                    
                            axs.plot(psi_cut, S_plot, color= color_dic[ad], label = f'{ad} $10^{{{20}}}$ 1/s')


                            # Labels and title
                            axs.set_xlabel('$\psi_N$')
                            axs.legend(loc = 'upper center')
                            if scan_style == 'temperature':
                                
                                fig.suptitle(f'Source {scan_style} scan from heat flux {st_val} to {ed_val} $10^5$ W')
                            
                            elif scan_style == 'density':
                                
                                fig.suptitle(f'Source {scan_style} scan from particle flux {st_val} to {ed_val} $10^{{{20}}}$ 1/s')
                    
                    
                    if format_option == 'neuden':
                        
                        
                        if plot_case == 'fivescan':
                            
                            if log_scale:
                                axs.set_yscale('log')
                            else:
                                pass
                                    
                            axs.plot(psi_cut, nu_plot, color= color_dic[ad], label = f'{ad}')


                            # Labels and title
                            axs.set_xlabel('$\psi_N$')
                            if scan_style == 'temperature':
                                
                                fig.suptitle(f'Neutral density {scan_style} scan from heat flux {st_val} to {ed_val} $10^5$ W')
                            
                            elif scan_style == 'density':
                                
                                fig.suptitle(f'Neutral density {scan_style} scan from particle flux {st_val} to {ed_val} $10^{{{20}}}$ 1/s')
                        
                    
                    
                    
                    
                    
                    
                
                    
                   
            
            
                
                        
                        
                        
                        
            
                    
                        
                    
                
                    
                    
                    
                
               












            
            
            

                
                

                    
                    


















