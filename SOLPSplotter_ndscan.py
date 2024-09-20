# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 00:45:40 2024

@author: ychuang
"""




from SOLPSplotter_NTplot import NT_plot
import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.offsetbox import AnchoredText


class neuden_scan(NT_plot):
    
    def __init__(self, DefaultSettings, loadDS):
        NT_plot.__init__(self, DefaultSettings, loadDS)

    
    def twinscan_ndmid(self, iter_index, data_struc):
        
        result_dic = self.data['radial_fit_data'][iter_index]
        
        
        if self.series_flag == 'twin_scan':
            
            nf = iter_index[0]
            tf = iter_index[1]
            
            # b2fstate = self.data['b2fstate'][nf][tf]
            
            # nx = data_struc['nx']
            # ny = data_struc['ny']
            
            if data_struc['size'] == 'full':
                neu_pro = self.data['outputdata']['NeuDen'][nf][tf]
                weight = self.data['midplane_calc']['weight']
                psi_coord = self.data['midplane_calc']['psi_solps_mid']
            elif data_struc['size'] == 'small':
                nx = data_struc['nx']
                ny = data_struc['ny']
                data = self.data['ft44'][nf][tf]['dab2']
                neu_pro = np.transpose(data[:, :, 0])
                weight = self.data['midplane_calc']['weight'][1:ny+1]
                psi_coord = self.data['midplane_calc']['psi_solps_mid'][1:ny+1]
            
        else:
            
            if data_struc['size'] == 'full':
                neu_pro = self.data['outputdata']['NeuDen'][iter_index]
                weight = self.data['midplane_calc'][iter_index]['weight']
                psi_coord = self.data['midplane_calc'][iter_index]['psi_solps_mid']
                
            elif data_struc['size'] == 'small':
                nx = data_struc['nx']
                ny = data_struc['ny']
                data = self.data['ft44'][iter_index]['dab2']
                neu_pro = np.transpose(data[:, :, 0])
                weight = self.data['midplane_calc'][iter_index]['weight'][1:ny+1]
                psi_coord = self.data['midplane_calc'][iter_index]['psi_solps_mid'][1:ny+1]
            
            
    
        
        # neu_pro = self.data['outputdata']['NeuDen'][iter_index]
        
        # weight = self.data['midplane_calc'][iter_index]['weight']
        weight_B = np.ones(len(weight))- weight
        
        
        mid_neu_pro = np.multiply(neu_pro[:, 58], weight) + np.multiply(neu_pro[:, 60], weight_B)
        
        psi_list = []
        nd_list = []
        
        for ind, coord in enumerate(psi_coord):
            
            if coord >= 0.95 and coord <= 1.05:
                psi_list.append(coord)
                nd_list.append(mid_neu_pro[ind])
        
        
        
        return psi_list, nd_list

    

    def twinscan_ndrad_method(self, log_flag, cl_dic, A_dic
                                   ,scandetail, iterlist, dat_size, scan_style):
        
        
        fig, axs = plt.subplots()
        
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        
        
        if dat_size == 'full':
    
            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        elif dat_size == 'small':
            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        
        print('this is 201:')
        print(dat_size)
        
        for aa in iterlist:
            
            psi_list, nd_list = self.twinscan_ndmid(iter_index = aa, 
                                                            data_struc = dat_struc)
            
            
            """
            label= 'core density {} $10^{19}$'.format(aa)
            
            """
            
            axs.legend(loc= 'lower left', fontsize=10)
            
            if self.series_flag == 'twin_scan':
                
                if scan_style == 'tempscan':
                    
                    ad = aa[1]
                
                elif scan_style == 'denscan':
                    
                    ad = aa[0]
                
                else:
                    print('neteTSplot_method, please check scan_style')
            
            else:
                ad = aa
                
            
    
            if scan_style == 'denscan':
                
                # exp_an_fit = fit_dat['exp_fit']
                # xcoord_cut = fit_dat['x_coord_cut']
                
                if log_flag:
                    axs.set_yscale('log')
                else:
                    pass
                
                axs.plot(psi_list, nd_list,'-', color = cl_dic[ad], label= '{}'.format(A_dic[ad]))
                # axs.plot(xcoord_cut, exp_an_fit, color='r',lw= 5, ls='-', label= 'exponential fit')
                # axs.axvline(x= max(xcoord_cut), color='black',lw=3, ls='--', 
                #             label= 'fit range : $\Delta n_e$')
                # axs.axvline(x= min(xcoord_cut), color='black',lw=3, ls='--')
                axs.set_title('Density scan with Te = {} eV'.format(scandetail))
                axs.set_xlabel('$\psi_N$')
                
                axs.legend(loc= 'lower right')
    
                            
            elif scan_style == 'tempscan':
                
    
                # exp_an_fit = fit_dat['exp_fit']
                # xcoord_cut = fit_dat['x_coord_cut']
                
                if log_flag:
                    axs.set_yscale('log')
                else:
                    pass
                
                axs.plot(psi_list, nd_list,'-', color = cl_dic[ad], label= '{}'.format(A_dic[ad]))
                # axs.plot(xcoord_cut, exp_an_fit, color='r',lw= 5, ls='-', label= 'exponential fit')
                # axs.axvline(x= max(xcoord_cut), color='black',lw=3, ls='--', 
                #             label= 'fit range : $\Delta n_e$')
                # axs.axvline(x= min(xcoord_cut), color='black',lw=3, ls='--')
                axs.set_title('Temperature scan with Ne = {}'.format(scandetail))
                axs.set_xlabel('$\psi_N$')
                
                axs.legend(loc= 'lower right')
            
            else:
                print('neteTSplot_structure, please check the scan parameter')
        
    
    
    
    
    
    
    def twinscan_ndrad_plot(self, scan_style, dat_size, log_flag):
        
        
        
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
                    
                    
                    self.twinscan_ndrad_method(iterlist = iter_key, cl_dic = color_dic, 
                                A_dic = label_dic, scan_style = scan_style, log_flag = log_flag,
                                scandetail = scan_title, dat_size = dat_size)
            
            
        else:
            print('Opacity_study_radial_plot has a bug')