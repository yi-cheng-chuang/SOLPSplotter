# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 14:02:21 2024

@author: ychuang
"""

from SOLPSplotter_radial import radial_plot
import matplotlib.pyplot as plt
import fitting_method as fm 
from scipy import interpolate
from scipy.optimize import curve_fit
import numpy as np
from matplotlib.offsetbox import AnchoredText
import fitting_method as fm


class NT_plot(radial_plot):
    
    def __init__(self, DefaultSettings, loadDS):
        radial_plot.__init__(self, DefaultSettings, loadDS)
        
    
    
    def nete_midprof(self, itername, data_struc):
        


        if self.withshift == True and self.withseries == False:
            
            b2fstate = self.data['b2fstate'][itername]
        
        elif self.withshift == False and self.withseries == True:
            
            if self.series_flag == 'twin_scan':
                
                nf = itername[0]
                tf = itername[1]
                
                b2fstate = self.data['b2fstate'][nf][tf]
                
                nx = data_struc['nx']
                ny = data_struc['ny']
                
                if data_struc['size'] == 'full':
                    neu_pro = self.data['outputdata']['NeuDen'][nf][tf]
                elif data_struc['size'] == 'small':
                    data = self.data['ft44'][nf][tf]['dab2']
                    neu_pro = np.transpose(data[:, :, 0])
                
                
            
            else:
                
                b2fstate = self.data['b2fstate'][itername]
                neu_pro = self.data['outputdata']['NeuDen'][itername]
        
        nx = data_struc['nx']
        ny = data_struc['ny']
        
        if data_struc['size'] == 'full':
            ne_pro = b2fstate['ne'].transpose()
            Te_J = b2fstate['te'].transpose()
            
        elif data_struc['size'] == 'small':
            ne_pro = b2fstate['ne'][1:nx+1, 1:ny+1].transpose()
            Te_J = b2fstate['te'][1:nx+1, 1:ny+1].transpose()
            
        
            
        
        ev = 1.6021766339999999 * pow(10, -19)
        te_pro = Te_J / ev
        
        # neu_pro = np.transpose(data[:, :, 0])
        
        if self.withshift == True and self.withseries == False:
        
            leftcut = self.data['b2fgeo'][itername]['leftcut'][0]
            rightcut = self.data['b2fgeo'][itername]['rightcut'][0]
            weight = self.data['midplane_calc'][itername]['weight']
            psi_coord = self.data['midplane_calc'][itername]['psi_solps_mid']
            
            
        elif self.withshift == False and self.withseries == True:
        
            leftcut = self.data['b2fgeo']['leftcut'][0]
            rightcut = self.data['b2fgeo']['rightcut'][0]
            
            
            if data_struc['size'] == 'full':
                psi_coord = self.data['midplane_calc']['psi_solps_mid']
                weight = self.data['midplane_calc']['weight']
                
            elif data_struc['size'] == 'small':
                psi_coord = self.data['midplane_calc']['psi_solps_mid'][1:ny+1]
                weight = self.data['midplane_calc']['weight'][1:ny+1]
        
        else:
            print('nete_TSplotmethod geo cut has a bug!')
        
        weight_B = np.ones(len(weight))- weight
             
        mid_ne_pro = np.multiply(ne_pro[:, 58], weight) + np.multiply(ne_pro[:, 60], weight_B)
        mid_te_pro = np.multiply(te_pro[:, 58], weight) + np.multiply(te_pro[:, 60], weight_B)
        mid_neu_pro = np.multiply(neu_pro[:, 58], weight) + np.multiply(neu_pro[:, 60], weight_B)
        
        return psi_coord, mid_ne_pro, mid_te_pro
    
    
    
    def plot_neteTSdat(self):
        
        
        """     
        # if self.data['outputdata'].any() == None or self.data['outputdata']['Te'].any() == None:
        if 'Ne' and 'Te' and 'NeuDen' in self.data['outputdata']:
            pass
        else:
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
            self.load_output_data(param= 'NeuDen')
        
        ne_pro = self.data['outputdata']['Ne']
        te_pro = self.data['outputdata']['Te']
        neu_pro = self.data['outputdata']['NeuDen']
        
        """
        
        psiN = self.data['experimental_fit']['psiN']
        ne = self.data['experimental_fit']['ne']*pow(10, 20)
        te = self.data['experimental_fit']['te']*pow(10, 3)
        
        exp = self.data['ExpDict']
        # psi = exp['psi_normal']
        
        
        psi = []
        exp_ne = []
        ne_er = []
        exp_te = []
        te_er = []
        for ep in range(len(exp['psi_normal'])):
            
            if  exp['psi_normal'][ep] >= min(psiN):
                psi.append(exp['psi_normal'][ep])
                exp_ne.append(exp['electron_density(10^20/m^3)'][ep]*pow(10, 20))
                ne_er.append(exp['density error(10^20/m^3)'][ep]*pow(10, 20))
                exp_te.append(exp['electron_temperature(KeV)'][ep]*pow(10, 3))
                te_er.append(exp['temperature error(10^20/m^3)'][ep]*pow(10, 3))
        
        TS_dic = {'psi': psi, 'neTS': exp_ne, 'errne': ne_er,
                  'teTS': exp_te, 'errte': te_er}
        
        return TS_dic
    
    
    
    def neteTSplot_method(self, iterlist, cl_dic, A_dic, scan_style, scandetail, dat_size):
        
        TS_dic = self.plot_neteTSdat()
        
        psi = TS_dic['psi']
        exp_ne = TS_dic['neTS']
        ne_er = TS_dic['errne']
        exp_te = TS_dic['teTS']
        te_er = TS_dic['errte']
        
        
        fig, axs = plt.subplots(2, 1)
        
        anchored_text = AnchoredText('(a){}'.format('$n_e$ [$m^{-3}$]'), loc='upper right')
        axs[0].errorbar(psi, exp_ne, yerr= ne_er, fmt = 'o', color = 'black', label= '$n_e$ TS data')
        axs[0].add_artist(anchored_text)
        axs[0].legend(loc='lower left', fontsize=10)
        
        
        
        anchored_text2 = AnchoredText('(b){}'.format('$t_e$ [eV]'), loc= 'upper right')
        axs[1].errorbar(psi, exp_te, yerr= te_er, fmt = 'o', color = 'black', label= '$t_e$ TS data')
        axs[1].set_xlabel('$\psi_N$')
        axs[1].add_artist(anchored_text2)
        axs[1].legend(loc='lower left', fontsize=10)
        
        plt.subplots_adjust(hspace=.0)
        
        
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        
        
        if dat_size == 'full':

            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        elif dat_size == 'small':
            dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
        
        
        
        for aa in iterlist:
            
            psi_coord, mid_ne_pro, mid_te_pro = self.nete_midprof(itername = aa, 
                                                    data_struc = dat_struc)
            
            
            """
            label= 'core density {} $10^{19}$'.format(aa)
            
            """
            
            axs[0].legend(loc= 'lower left', fontsize=10)
            
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
                
                axs[0].plot(psi_coord, mid_ne_pro, color = cl_dic[ad])
                axs[0].set_title('Density scan with Te = {} eV'.format(scandetail))
                axs[1].plot(psi_coord, mid_te_pro, color = cl_dic[ad],
                            label= '{}'.format(A_dic[ad]))
                axs[1].legend()

                            
            elif scan_style == 'tempscan':
                
                axs[0].plot(psi_coord, mid_ne_pro, color = cl_dic[ad], 
                            label= '{}'.format(A_dic[ad]))
                axs[0].set_title('Temperature scan with Ne = {}'.format(scandetail))
                axs[1].plot(psi_coord, mid_te_pro, color = cl_dic[ad])
                axs[0].legend()
            
            else:
                print('neteTSplot_structure, please check the scan parameter')
            
            

        
        # fig.savefig('profiles.pdf')
    
    
    def pair_dic(self, keys, values):
        
        # Use zip() to pair the keys with the values
        zipped_pairs = zip(keys, values)
        
        # Convert the zipped pairs into a dictionary
        result_dic = dict(zipped_pairs)
        
        return result_dic
    
    
        
    def neteTS_plot(self, scan_style, dat_size):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            label_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            asp_ch = self.data['dircomp']['multi_shift']
            
            self.neteTSplot_structure(iterlist = asp_ch, 
                                      cl_dic = color_dic, A_dic = label_dic, scan = 'not')
        
        elif self.withshift == False and self.withseries == True:
            
            # series_flag = self.DefaultSettings['series_flag']
            
            
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
                
                keylist_a = [str(x) for x in dircomp[key_a]]
                
                for ta in keylist_a:
                    
                    keylist_b = [str(x) for x in dircomp[key_b]]
                    color_list = ['red', 'orange', 'green', 'blue', 'purple']
                    
                    color_dic = self.pair_dic(keys = keylist_b, values = color_list)
                    
                    scan_list = []
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
                            
                            
                            
                        psi_coord, mid_ne_pro, mid_te_pro = self.nete_midprof(itername = it_in, 
                                                                data_struc = dat_struc)
                        
                        
                        if scan_style == 'tempscan':
                            
                            scan_add = '{:.1f} eV'.format(mid_te_pro[0])
                        
                        elif scan_style == 'denscan':
                            
                            scan_add = '{:.2E}'.format(mid_ne_pro[0])
                        
                        else:
                            print('twinscan_plot_method, please check the scan_style!')
                        
                        scan_list.append(scan_add)
                        iter_key.append(it_in)
                        
                    if scan_style == 'tempscan':
                        psi_coord, mid_ne_pro, mid_te_pro = self.nete_midprof(itername = (ta, '3.73'),
                                                                   data_struc = dat_struc)
                        scan_title = '{:.2E}'.format(mid_ne_pro[0])
                    
                    elif scan_style == 'denscan':
                        psi_coord, mid_ne_pro, mid_te_pro = self.nete_midprof(itername = ('5.02', ta), 
                                                            data_struc = dat_struc)
                        scan_title = '{:.1f}'.format(mid_te_pro[0])
                    
                    else:
                        print('twinscan_plot_method, please check the scan_style!')
                    
                    label_dic = self.pair_dic(keys = keylist_b, values = scan_list)
                    
                    self.neteTSplot_method(iterlist = iter_key, scandetail = scan_title,
                            cl_dic = color_dic, A_dic = label_dic, 
                            scan_style = scan_style, dat_size = dat_size)
                        
                    # print(scan_list)
             
            else:
                print('neteTS_plot, please check the series flag')
                
    
        
    def AM_NT_midprof(self, itername, AM_flag):
        
        if AM_flag == 'atom':
            
            den = 'dab2'
            temp = 'tab2'
            
        elif AM_flag == 'mol':
            
            den = 'dmb2'
            temp = 'tmb2'
        
        ev = 1.6021766339999999 * pow(10, -19)
        
        if self.withshift == True and self.withseries == False:
            
            neu_data = self.data['ft44'][itername][den]
            neu_pro = np.transpose(neu_data[:, :, 0])
            atom_temp_data = self.data['ft44'][itername][temp]
            
            atom_temp = np.transpose(atom_temp_data[:, :, 0])
            atom_temp_pro = atom_temp / ev
        
        elif self.withshift == False and self.withseries == True:
            
            if self.series_flag == 'twin_scan':
                
                nf = itername[0]
                tf = itername[1]
                
                neu_data = self.data['ft44'][nf][tf][den]
                neu_pro = np.transpose(neu_data[:, :, 0])
                atom_temp_data = self.data['ft44'][nf][tf][temp]
                atom_temp = np.transpose(atom_temp_data[:, :, 0])
                atom_temp_pro = atom_temp / ev
                              
            else:
                
                neu_data = self.data['ft44'][itername][den]
                neu_pro = np.transpose(neu_data[:, :, 0])
                atom_temp_data = self.data['ft44'][itername][temp]
                atom_temp = np.transpose(atom_temp_data[:, :, 0])
                atom_temp_pro = atom_temp / ev
        
        if self.withshift == True and self.withseries == False:
        
            leftcut = self.data['b2fgeo'][itername]['leftcut'][0]
            rightcut = self.data['b2fgeo'][itername]['rightcut'][0]
            weight = self.data['midplane_calc'][itername]['weight']
            psi_coord = self.data['midplane_calc'][itername]['psi_solps_mid']
        
        elif self.withshift == False and self.withseries == True:
        
            leftcut = self.data['b2fgeo']['leftcut'][0]
            rightcut = self.data['b2fgeo']['rightcut'][0]
            weight = self.data['midplane_calc']['weight'][1:37]
            psi_coord = self.data['midplane_calc']['psi_solps_mid'][1:37]
        
        else:
            print('NeuDen_plotmethod, please check withshift and withseries flag')
        
        weight_B = np.ones(len(weight))- weight
        
        mid_neu_pro = np.multiply(neu_pro[:, 58], weight) + np.multiply(neu_pro[:, 60], weight_B)
        mid_atom_temp_pro = np.multiply(atom_temp_pro[:, 58], weight) + np.multiply(atom_temp_pro[:, 60], weight_B)
        
        return psi_coord, mid_neu_pro, mid_atom_temp_pro
    
    
    
    def AM_NTplot_method(self, iterlist, cl_dic, A_dic, AM_flag, scandetail, scan_style):
        
        
        fig, axs = plt.subplots(2, 1)
        
        for aa in iterlist:
            
            psi_coord, mid_neu_pro, mid_atom_temp_pro = self.AM_NT_midprof(itername = aa, 
                                                                AM_flag= AM_flag)
            
            if AM_flag == 'atom':
                
                anchored_text = AnchoredText('(a){}'.format('atomic density [$m^{-3}$]'), loc='upper center')
                anchored_text2 = AnchoredText('(b){}'.format('atomic temperature [eV]'), loc= 'upper center')
            
            elif AM_flag == 'mol':
                
                anchored_text = AnchoredText('(a){}'.format('molecular density [$m^{-3}$]'), loc= 'upper center')
                anchored_text2 = AnchoredText('(b){}'.format('molecular temperature [eV]'), loc= 'upper center')
                

            
            if self.series_flag == 'twin_scan':
                
                if scan_style == 'tempscan':
                    
                    ad = aa[1]
                
                elif scan_style == 'denscan':
                    
                    ad = aa[0]
                
                else:
                    print('AM_NTplot_method, please check scan_style')
                
                
            
            else:
                ad = aa

            if scan_style == 'denscan':
                
                axs[0].set_title('Density scan with Te = {} eV'.format(scandetail))
                axs[0].plot(psi_coord, mid_neu_pro, color = cl_dic[ad], 
                            label= '{}'.format(A_dic[ad]))
                axs[1].plot(psi_coord, mid_atom_temp_pro, color = cl_dic[ad])
                axs[0].add_artist(anchored_text)
                axs[1].add_artist(anchored_text2)
                axs[1].set_xlabel('$\psi_N$')
                axs[0].legend()

                            
            elif scan_style == 'tempscan':
                
                
                axs[0].set_title('Temperature scan with Ne = {}'.format(scandetail))
                axs[0].plot(psi_coord, mid_neu_pro, color = cl_dic[ad])
                axs[1].plot(psi_coord, mid_atom_temp_pro, color = cl_dic[ad], 
                            label= '{}'.format(A_dic[ad]))
                axs[0].add_artist(anchored_text)
                axs[1].add_artist(anchored_text2)
                axs[1].set_xlabel('$\psi_N$')
                axs[1].legend()
            
            else:
                print('AM_NTplot_method, please check the scan parameter')
    
            plt.subplots_adjust(hspace=.0)
    
    
    def twinscan_plot_method(self, dircomp, scan_style, AM_flag):
        
        
        if scan_style == 'tempscan':
            
            key_a = 'denscan_list'
            key_b = 'tempscan_list'
        
        elif scan_style == 'denscan':
            
            key_a = 'tempscan_list'
            key_b = 'denscan_list'
        
        else:
            print('twinscan_plot_method, please check the scan_style!')
        
        keylist_a = [str(x) for x in dircomp[key_a]]
        
        for ta in keylist_a:
            
            keylist_b = [str(x) for x in dircomp[key_b]]
            color_list = ['red', 'orange', 'green', 'blue', 'purple']
            
            color_dic = self.pair_dic(keys = keylist_b, values = color_list)
            
            scan_list = []
            iter_key = []
            
            for tb in keylist_b:
                
                
                if scan_style == 'tempscan':
                    
                    it_in = (ta, tb)
                
                elif scan_style == 'denscan':
                    
                    it_in = (tb, ta)
                
                else:
                    print('twinscan_plot_method, please check the scan_style!')
                    
                    
                psi_coord, mid_ne_pro, mid_te_pro = self.nete_midprof(itername = it_in)
                
                if scan_style == 'tempscan':
                    
                    scan_add = '{:.1f} eV'.format(mid_te_pro[0])
                
                elif scan_style == 'denscan':
                    
                    scan_add = '{:.2E}'.format(mid_ne_pro[0])
                
                else:
                    print('twinscan_plot_method, please check the scan_style!')
                    
                
                
                scan_list.append(scan_add)
                iter_key.append(it_in)
            
            
            if scan_style == 'tempscan':
                psi_coord, mid_ne_pro, mid_te_pro = self.nete_midprof(itername = (ta, '3.73'))
                scan_title = '{:.2E}'.format(mid_ne_pro[0])
            
            elif scan_style == 'denscan':
                psi_coord, mid_ne_pro, mid_te_pro = self.nete_midprof(itername = ('5.02', ta))
                scan_title = '{:.1f}'.format(mid_te_pro[0])
            
            else:
                print('twinscan_plot_method, please check the scan_style!')
            
            label_dic = self.pair_dic(keys = keylist_b, values = scan_list)
            
            self.AM_NTplot_method(iterlist = iter_key, AM_flag = AM_flag,
                    cl_dic = color_dic, A_dic = label_dic,
                    scandetail = scan_title, scan_style = scan_style)
                
            print(scan_list)
    
    
    
    
    def AtomNT_plot(self, AM_flag, scan_style):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            label_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            asp_ch = self.data['dircomp']['multi_shift']
            
            self.AtomNTplot_method(iterlist = asp_ch, 
                                      cl_dic = color_dic, A_dic = label_dic, scan = 'not')
        
        elif self.withshift == False and self.withseries == True:
            
            
            if self.series_flag == 'twin_scan':
                
                dircomp = self.data['dircomp']
                
                self.twinscan_plot_method(dircomp = dircomp, scan_style = scan_style, 
                                    AM_flag = AM_flag)
            
             
            else:
                print('AtomNT_plot, please check the series flag')
    



"""

ds_key = [str(x) for x in dircomp['denscan_list']]

for td in ds_key:
    
    ts_key = [str(x) for x in dircomp['tempscan_list']]
    color_list = ['red', 'orange', 'green', 'blue', 'purple']
    
    color_dic = self.pair_dic(keys = ts_key, values = color_list)
    
    temp_list = []
    iter_key = []
    
    for ts in ts_key:
        
        it_in = (td, ts)
        psi_coord, mid_ne_pro, mid_te_pro = self.nete_midprof(itername = it_in)
        
        
        temp_add = '{:.1f} eV'.format(mid_te_pro[0])
        
        temp_list.append(temp_add)
        iter_key.append(it_in)
        
    scan_title = '{:.2E}'.format(mid_ne_pro[0])
    
    label_dic = self.pair_dic(keys = ts_key, values = temp_list)
    
    self.AM_NTplot_method(iterlist = iter_key, AM_flag = AM_flag,
            cl_dic = color_dic, A_dic = label_dic, scan = 'tempscan',
            scandetail = scan_title)
        
    print(temp_list)


if self.series_flag == 'change_den':
    
    color_dic = {'3.4': 'red', '4.0': 'orange', '5.0': 'green',
                 '6.0': 'blue', '7.0': 'purple'}
    
    label_dic = {'3.4': '4.15*$10^{19} m^{-3}$', 
'4.0': '4.37*$10^{19} m^{-3}$', '5.0': '4.62*$10^{19} m^{-3}$',
'6.0': '4.9*$10^{19} m^{-3}$', '7.0': '5.16*$10^{19} m^{-3}$'}
    
    denscan = list(self.data['dircomp']['Attempt'].keys())
    
    self.neteTSplot_structure(iterlist = denscan, 
                cl_dic = color_dic, A_dic = label_dic, scan = 'den')

elif self.series_flag == 'change_temp':
    
    color_dic = {'2.5': 'red', '3.0': 'orange', '4.0': 'green',
                 '5.0': 'blue', '6.0': 'purple'}
    
    label_dic = {'2.5': '477 eV', '3.0': '571 eV', '4.0': '732 eV',
        '5.0': '884 eV', '6.0': '1037 eV'}
    
    denscan = list(self.data['dircomp']['Attempt'].keys())
    
    self.neteTSplot_structure(iterlist = denscan, 
                cl_dic = color_dic, A_dic = label_dic, scan = 'temp')

elif self.series_flag == 'terminal_test':
    
    color_dic = {'4.15': 'red', '5.0': 'orange', '6.0': 'green',
                 '7.0': 'blue', '8.0': 'purple', '9.0': 'brown'}
    
    label_dic = {'4.15': '4.15*$10^{19} m^{-3}$', 
'5.0': '5.0*$10^{19} m^{-3}$', '6.0': '6.0*$10^{19} m^{-3}$',
'7.0': '7.0*$10^{19} m^{-3}$', '8.0': '8.0*$10^{19} m^{-3}$',
'9.0': '9.0*$10^{19} m^{-3}$'}
    
    dencan = list(self.data['dircomp']['Attempt'].keys())
    
    
    
    self.neteTSplot_structure(iterlist = denscan, 
                cl_dic = color_dic, A_dic = label_dic, scan = 'den')

"""
    