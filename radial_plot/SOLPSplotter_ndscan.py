# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 00:45:40 2024

@author: ychuang
"""


from SOLPSplotter_NTplot import NT_plot
from SOLPSplotter_ioutflux import iout_flux
import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.offsetbox import AnchoredText
from scipy import interpolate


class neuden_scan(NT_plot, iout_flux):
    
    def __init__(self, DefaultSettings, loadDS):
        NT_plot.__init__(self, DefaultSettings, loadDS)
        iout_flux.__init__(self, DefaultSettings, loadDS)

    

    def twinscan_ndrad_method(self, log_flag, cl_dic, A_dic, xcoord_type,
                                   scandetail, iterlist, dat_size, scan_style, format_option):
        
        if format_option == '2x1':
            fig, axs = plt.subplots(2, 1)
        
        elif format_option == 'neuden':
            fig, axs = plt.subplots()
        
        elif format_option == 'ionize':
            fig, axs = plt.subplots()
        
            
        midplane_psi = self.data['midplane_calc']['psi_solps_mid']
        r_rsep = self.data['midplane_calc']['R_Rsep']
             
        psi_to_dsa_func = interpolate.interp1d(midplane_psi, r_rsep, fill_value = 'extrapolate')
        
        
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        
        
        if dat_size == 'full':
    
            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        elif dat_size == 'small':
            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        plt.subplots_adjust(hspace=.0)
        anchored_text = AnchoredText('{}'.format('Neutral density [$m^{-3}$]'), loc='upper left')
        anchored_text_2 = AnchoredText('{}'.format('Source [$s^{-1}$]'), loc='upper left')
        # print('this is 201:')
        # print(dat_size)
        
        for aa in iterlist:
            
            psi_list, nd_list = self.twinscan_ndmid(iter_index = aa, 
                                                            data_struc = dat_struc)
            
            psi_list, S_list = self.twinscan_Smid(iter_index = aa, 
                                                            data_struc = dat_struc)
            
            rrsep_solps = psi_to_dsa_func(psi_list)
            
            
            nd_array = np.array(nd_list)
            S_array = np.array(S_list)
            S_norm = max(S_list)
            # print(f'the length of nd_array is: {len(nd_array)}')
            
            target = 1
            
            # Find the minimum distance to the target
            min_diff = min(abs(x - target) for x in psi_list)
            
            # Find all values and indices that have this minimum difference
            closest_elements = [(i, x) for i, x in enumerate(psi_list) if abs(x - target) == min_diff]
            
            # Print results
            print(f"Closest values to {target:.4f} (within ±{min_diff:.4f}):")
            for index, value in closest_elements:
                print(f"Value: {value:.4f}, Index: {index}")
            
            if value < 1:
                
                nd_norm = 0.5*(nd_list[index] + nd_list[index + 1])
                
                
            elif value > 1:
                
                nd_norm = 0.5*(nd_list[index] + nd_list[index - 1])
            
            else:
                print('please check twinscan_ndrad_method function in SOLPSplotter_ndscan.py')
                
            

            norm_nd_array = nd_array / nd_norm
            norm_S_array = S_array / S_norm
            
            
            """
            label= 'core density {} $10^{19}$'.format(aa)
            
            """
            
            
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
                
                
            if scan_style == 'denscan':
                
                title_ap = float(ap)*pow(10, 5)
                label_ad = float(ad)*pow(10, 20)
                
                if format_option == '2x1':
                    
                    if log_flag:
                        axs[0].set_yscale('log')
                    else:
                        pass
                    
                    
                    axs[0].plot(psi_list, norm_nd_array,'-', color = cl_dic[ad], label= '{:.3E} (1/s)'.format(label_ad))
                    axs[1].plot(psi_list, norm_S_array,'-', color = cl_dic[ad], label= '{:.3E} (1/s)'.format(label_ad))
                    axs[1].set_xlabel('$\psi_N$')
                    axs[0].axvline(x= 1, color='black', lw=3, ls='--')
                    axs[1].axvline(x= 1, color='black', lw=3, ls='--')
                    axs[0].add_artist(anchored_text)
                    axs[1].add_artist(anchored_text_2)
                    axs[0].set_title('Particle flux scan with heat flux = {:.3E} W'.format(title_ap))
                    
                    axs[0].legend(loc= 'lower right')
                
                elif format_option == 'neuden':
                    
                    if log_flag:
                        axs.set_yscale('log')
                    else:
                        pass
                    
                    if xcoord_type == 'psi':
                        
                        
                        axs.plot(psi_list, norm_nd_array,'-', color = cl_dic[ad], label= '{:.3E} (1/s)'.format(label_ad))
                        axs.set_xlabel('$\psi_N$')
                    
                    elif xcoord_type == 'rrsep':
                        
                        axs.plot(rrsep_solps, norm_nd_array,'-', color = cl_dic[ad], label= '{:.3E} (1/s)'.format(label_ad))
                        axs.set_xlabel('$R - R_{sep}$')
                    
                    
                    
                    axs.add_artist(anchored_text)
                    axs.set_title('Particle flux scan with heat flux = {:.3E} W'.format(title_ap))
                
                elif format_option == 'ionize':
                    
                    if log_flag:
                        axs.set_yscale('log')
                    else:
                        pass
                    
                    if xcoord_type == 'psi':
                        
                        axs.plot(psi_list, norm_S_array,'-', color = cl_dic[ad], label= '{:.3E} (1/s)'.format(label_ad))
                        axs.set_xlabel('$\psi_N$')
                    
                    elif xcoord_type == 'rrsep':
                        
                        axs.plot(rrsep_solps, norm_S_array,'-', color = cl_dic[ad], label= '{:.3E} (1/s)'.format(label_ad))
                        axs.set_xlabel('$R - R_{sep}$')
                    
                    
                    
                    # axs.add_artist(anchored_text_2)
                    axs.set_title('Particle flux scan with heat flux = {:.3E} W'.format(title_ap))
                    
                            
            elif scan_style == 'tempscan':
                
                title_ap = float(ap)*pow(10, 20)
                label_ad = float(ad)*pow(10, 5)
                # exp_an_fit = fit_dat['exp_fit']
                # xcoord_cut = fit_dat['x_coord_cut']
                if format_option == '2x1':
                    
                    if log_flag:
                        axs[0].set_yscale('log')
                    else:
                        pass
                    
                    axs[0].plot(psi_list, norm_nd_array,'-', color = cl_dic[ad], label= '{:.3E} W'.format(label_ad))
                    axs[1].plot(psi_list, norm_S_array,'-', color = cl_dic[ad], label= '{:.3E} W'.format(label_ad))
                    axs[0].axvline(x= 1, color='black', lw=3, ls='--')
                    axs[1].axvline(x= 1, color='black', lw=3, ls='--')
                    axs[0].add_artist(anchored_text)
                    axs[1].add_artist(anchored_text_2)
                    axs[0].set_title('Heat flux scan with particle flux = {:.3E} (1/s)'.format(title_ap))
                    axs[1].set_xlabel('$\psi_N$')
                    axs[0].legend(loc= 'lower right')
                
                
                elif format_option == 'neuden':
                    
                    if log_flag:
                        axs.set_yscale('log')
                    else:
                        pass
                    
                    if xcoord_type == 'psi':
                        
                        axs.plot(psi_list, norm_nd_array,'-', color = cl_dic[ad], label= '{:.3E} W'.format(label_ad))
                        axs.set_xlabel('$\psi_N$')
                    
                    elif xcoord_type == 'rrsep':
                        
                        axs.plot(rrsep_solps, norm_nd_array,'-', color = cl_dic[ad], label= '{:.3E} W'.format(label_ad))
                        axs.set_xlabel('$R - R_{sep}$')
                    
                    
                    axs.add_artist(anchored_text)
                    axs.set_title('Heat flux scan with particle flux = {:.3E} (1/s)'.format(title_ap))
                
                elif format_option == 'ionize':
                    
                    if log_flag:
                        axs.set_yscale('log')
                    else:
                        pass
                    
                    if xcoord_type == 'psi':
                        
                        axs.plot(psi_list, norm_S_array,'-', color = cl_dic[ad], label= '{:.3E} W'.format(label_ad))
                        axs.set_xlabel('$\psi_N$')
                    
                    elif xcoord_type == 'rrsep':
                        
                        axs.plot(rrsep_solps, norm_S_array,'-', color = cl_dic[ad], label= '{:.3E} W'.format(label_ad))
                        axs.set_xlabel('$R - R_{sep}$')
                    
                    
                    # axs.add_artist(anchored_text_2)
                    axs.set_title('Heat flux scan with particle flux = {:.3E} (1/s)'.format(title_ap))
                    
                    
            else:
                print('neteTSplot_structure, please check the scan parameter')
        
    
        if format_option == 'neuden':
            
            
            if xcoord_type == 'psi':
                
                n_sym = self.data['radial_fit_data'][aa]['ne_symmetry_point']
                dn = self.data['opacity_poloidal'][aa]['pedestal_width_psiN']
                axs.axvline(x= n_sym + dn, color='gray',lw=3, label= '$\Delta n_e$')
                axs.axvline(x= n_sym - dn, color='gray',lw=3)
            
            elif xcoord_type == 'rrsep':
                
                n_sym_psi = self.data['radial_fit_data'][aa]['ne_symmetry_point']
                n_sym = psi_to_dsa_func(n_sym_psi)
                dn = self.data['opacity_poloidal'][aa]['pedestal_width']
                axs.axvline(x= n_sym + dn, color='gray',lw=3, label= '$\Delta n_e$')
                axs.axvline(x= n_sym - dn, color='gray',lw=3)
            
            
            axs.legend(loc= 'lower right')
        
        else:
            pass
            
    
    
    
    
    def twinscan_ndrad_plot(self, scan_style, dat_size, log_flag, format_option, xcoord_type):
        
        
        
        if self.withshift == True and self.withseries == False:
            
            
            for aa in self.data['dircomp']['multi_shift']:
                
                psi_list, nd_list, fit_dat = self.neuden_midprof(iter_index = aa)
                
    
                self.paper_neuden_radial_method(fit_dat= fit_dat, itername = aa,
                    x_coord = psi_list, Nd = nd_list, log_flag = True)
                
                
        
        elif self.withshift == False and self.withseries == True:
            print('Opacity_study_radial_plot is not there yet, to be continue...')  
            
            
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
                    
                    color_dic = self.pair_dic(keys = keylist_b, values = color_list)
                    
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
                        
                        
                        nx = self.data['b2fgeo']['nx']
                        ny = self.data['b2fgeo']['ny']
                        
                        
                        if dat_size == 'full':
            
                            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
                        
                        elif dat_size == 'small':
                            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
                            
                            
                        psi_coord, mid_ne_pro, mid_te_pro, mid_neu_pro= self.nete_midprof(itername = it_in, 
                                                                data_struc = dat_struc)
                        
                        
                        
                        if scan_style == 'tempscan':
                            
                            scan_add = '{:.1f} eV'.format(mid_te_pro[0])
                        
                        elif scan_style == 'denscan':
                            
                            scan_add = '{:.2E} '.format(mid_ne_pro[0])
                        
                        else:
                            print('twinscan_plot_method, please check the scan_style!')
                        
                        scan_list.append(scan_add)
                        iter_key.append(it_in)
                    
                    
                    print('NT scan list: {}'.format(ta))
                    print(scan_list)
                    
                    
                    if scan_style == 'tempscan':
                        psi_coord, mid_ne_pro, mid_te_pro, mid_neu_pro= self.nete_midprof(itername = (ta, '4.115'),
                                                                   data_struc = dat_struc)
                        scan_title = '{:.2E}'.format(mid_ne_pro[0])
                    
                    elif scan_style == 'denscan':
                        psi_coord, mid_ne_pro, mid_te_pro, mid_neu_pro= self.nete_midprof(itername = ('5.512', ta), 
                                                            data_struc = dat_struc)
                        scan_title = '{:.1f}'.format(mid_te_pro[0])
                    
                    else:
                        print('twinscan_plot_method, please check the scan_style!')
                    
                    label_dic = self.pair_dic(keys = keylist_b, values = scan_list)
                    
                    
                    self.twinscan_ndrad_method(iterlist = iter_key, cl_dic = color_dic, xcoord_type = xcoord_type,
                                A_dic = label_dic, scan_style = scan_style, log_flag = log_flag,
                                scandetail = scan_title, dat_size = dat_size,format_option = format_option)
            
            
        else:
            print('Opacity_study_radial_plot has a bug')



