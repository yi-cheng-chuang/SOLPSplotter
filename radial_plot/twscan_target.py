# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 11:00:25 2024

@author: ychuang
"""


from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
from twscan_module.twinscan_prepare import twscan_assist



class twinscan_target_radial:
    
    
    def __init__(self, DF, data, twa: twscan_assist):
        
        self.DF = DF
        self.data = data
        self.twa = twa
 
    
    def twscan_targetNT_method(self, log_flag, cl_dic, A_dic, scandetail, iterlist, scan_style, side):
        
        fig, axs = plt.subplots(1, 3)
            
        
        plt.subplots_adjust(hspace=.0)
        if side == 'inner target':
            
            anchored_text_1 = AnchoredText('{}'.format('inner target electron density [$m^{-3}$]'), loc='upper left')
            anchored_text_2 = AnchoredText('{}'.format('inner target electron temperature [eV]'), loc='upper left')
            anchored_text_3 = AnchoredText('{}'.format('inner target source [$m^{-3}*s^{-1}$]'), loc='upper left')
            
        elif side == 'outer target':
            
            anchored_text_1 = AnchoredText('{}'.format('outer target electron density [$m^{-3}$]'), loc='upper left')
            anchored_text_2 = AnchoredText('{}'.format('outer target electron temperature [eV]'), loc='upper left')
            anchored_text_3 = AnchoredText('{}'.format('outer target source [$m^{-3}*s^{-1}$]'), loc='upper left')
        

        # print('this is 201:')
        
        for aa in iterlist:
            
            psi_dic = self.data['target_profile'][aa[0]][aa[1]]['psiN']
            ne_dic = self.data['target_profile'][aa[0]][aa[1]]['ne']
            te_dic = self.data['target_profile'][aa[0]][aa[1]]['te']
            sx_dic = self.data['target_profile'][aa[0]][aa[1]]['source']
            neuden_dic = self.data['target_profile'][aa[0]][aa[1]]['neuden']
            
            psi_list = psi_dic[side]
            ne_list = ne_dic[side]
            te_list = te_dic[side]
            sx_list = sx_dic[side]
           
            
            if self.DF.series_flag == 'twin_scan':
                
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
                
                    
                if log_flag:
                    axs[0].set_yscale('log')
                    axs[2].set_yscale('log')
                else:
                    pass
                
                
                axs[0].plot(psi_list, ne_list,'-', color = cl_dic[ad], label= '{:.3E} (1/s)'.format(label_ad))
                axs[1].plot(psi_list, te_list,'-', color = cl_dic[ad], label= '{:.3E} (1/s)'.format(label_ad))
                axs[2].plot(psi_list, sx_list,'-', color = cl_dic[ad], label= '{:.3E} (1/s)'.format(label_ad))
                axs[0].set_xlabel('$\psi_N$')
                axs[1].set_xlabel('$\psi_N$')
                axs[2].set_xlabel('$\psi_N$')
                axs[0].axvline(x= 1, color='black', lw=3, ls='--')
                axs[1].axvline(x= 1, color='black', lw=3, ls='--')
                axs[2].axvline(x= 1, color='black', lw=3, ls='--')
                axs[0].add_artist(anchored_text_1)
                axs[1].add_artist(anchored_text_2)
                axs[2].add_artist(anchored_text_3)
                axs[0].set_title('Particle flux scan with heat flux = {:.3E} W'.format(title_ap))
                
                axs[0].legend(loc= 'lower right')
                


            elif scan_style == 'tempscan':
                
                title_ap = float(ap)*pow(10, 20)
                label_ad = float(ad)*pow(10, 5)
                # exp_an_fit = fit_dat['exp_fit']
                # xcoord_cut = fit_dat['x_coord_cut']

                    
                if log_flag:
                    axs[0].set_yscale('log')
                    axs[2].set_yscale('log')
                else:
                    pass
                
                axs[0].plot(psi_list, ne_list,'-', color = cl_dic[ad], label= '{:.3E} W'.format(label_ad))
                axs[1].plot(psi_list, te_list,'-', color = cl_dic[ad], label= '{:.3E} W'.format(label_ad))
                axs[2].plot(psi_list, sx_list,'-', color = cl_dic[ad], label= '{:.3E} (1/s)'.format(label_ad))
                axs[0].axvline(x= 1, color='black', lw=3, ls='--')
                axs[1].axvline(x= 1, color='black', lw=3, ls='--')
                axs[2].axvline(x= 1, color='black', lw=3, ls='--')
                axs[0].add_artist(anchored_text_1)
                axs[1].add_artist(anchored_text_2)
                axs[2].add_artist(anchored_text_3)
                axs[0].set_title('Heat flux scan with particle flux = {:.3E} (1/s)'.format(title_ap))
                axs[0].set_xlabel('$\psi_N$')
                axs[1].set_xlabel('$\psi_N$')
                axs[2].set_xlabel('$\psi_N$')
                axs[0].legend(loc= 'lower right')
                    
                    
            else:
                print('neteTSplot_structure, please check the scan parameter')
        
    
    def twscan_targetNT_combine_method(self, log_flag, cl_dic, A_dic,
                                   scandetail, iterlist, scan_style):
        
        fig, axs = plt.subplots(1, 2)


        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        dat_struc = {'nx': nx, 'ny': ny}
        
        
        
        plt.subplots_adjust(hspace=.0)
            
        anchored_text_1 = AnchoredText('{}'.format('Electron density at target [$m^{-3}$]'), loc='upper left')
        anchored_text_2 = AnchoredText('{}'.format('Electron temperature at target [eV]'), loc='upper left')
        # anchored_text_3 = AnchoredText('{}'.format('Source at target [$m^{-3}*s^{-1}$]'), loc='upper left')
        
        
        for aa in iterlist:
            
            
            psi_dic = self.data['target_profile'][aa[0]][aa[1]]['psiN']
            ne_dic = self.data['target_profile'][aa[0]][aa[1]]['ne']
            te_dic = self.data['target_profile'][aa[0]][aa[1]]['te']
            
            psi_list_in = psi_dic['inner target']
            ne_list_in = ne_dic['inner target']
            te_list_in = te_dic['inner target']
            # sx_list_in = sx_dic['inner target']
            
            psi_list_out = psi_dic['outer target']
            ne_list_out = ne_dic['outer target']
            te_list_out = te_dic['outer target']
            # sx_list_out = sx_dic['outer target']
            
            
            if self.DF.series_flag == 'twin_scan':
                
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
                
                    
                if log_flag:
                    axs[0].set_yscale('log')
                else:
                    pass
                
                
                axs[0].plot(psi_list_in, ne_list_in,'-', color = cl_dic[ad], label= 'HFS {:.3E} (1/s)'.format(label_ad))
                axs[1].plot(psi_list_in, te_list_in,'-', color = cl_dic[ad], label= 'HFS {:.3E} (1/s)'.format(label_ad))
                # axs[2].plot(psi_list_in, sx_list_in,'-', color = cl_dic[ad], label= 'inner {:.3E} (1/s)'.format(label_ad))
                axs[0].plot(psi_list_out, ne_list_out,'--', color = cl_dic[ad], label= 'LFS {:.3E} (1/s)'.format(label_ad))
                axs[1].plot(psi_list_out, te_list_out,'--', color = cl_dic[ad], label= 'LFS {:.3E} (1/s)'.format(label_ad))
                # axs[2].plot(psi_list_out, sx_list_out,'--', color = cl_dic[ad], label= 'outer {:.3E} (1/s)'.format(label_ad))
                axs[0].set_xlabel('$\psi_N$')
                axs[1].set_xlabel('$\psi_N$')
                # axs[2].set_xlabel('$\psi_N$')
                axs[0].axvline(x= 1, color='black', lw=3, ls='--')
                axs[1].axvline(x= 1, color='black', lw=3, ls='--')
                # axs[2].axvline(x= 1, color='black', lw=3, ls='--')
                axs[0].add_artist(anchored_text_1)
                axs[1].add_artist(anchored_text_2)
                # axs[2].add_artist(anchored_text_3)
                axs[0].set_title('Particle flux scan with heat flux = {:.3E} W'.format(title_ap))
                
                axs[0].legend(loc= 'lower right')
                


            elif scan_style == 'tempscan':
                
                title_ap = float(ap)*pow(10, 20)
                label_ad = float(ad)*pow(10, 5)
                # exp_an_fit = fit_dat['exp_fit']
                # xcoord_cut = fit_dat['x_coord_cut']

                    
                if log_flag:
                    axs[0].set_yscale('log')
                else:
                    pass
                
                axs[0].plot(psi_list_in, ne_list_in,'-', color = cl_dic[ad], label= 'HFS {:.3E} W'.format(label_ad))
                axs[1].plot(psi_list_in, te_list_in,'-', color = cl_dic[ad], label= 'HFS {:.3E} W'.format(label_ad))
                # axs[2].plot(psi_list_in, sx_list_in,'-', color = cl_dic[ad], label= 'inner {:.3E} W'.format(label_ad))
                axs[0].plot(psi_list_out, ne_list_out,'--', color = cl_dic[ad], label= 'LFS {:.3E} W'.format(label_ad))
                axs[1].plot(psi_list_out, te_list_out,'--', color = cl_dic[ad], label= 'LFS {:.3E} W'.format(label_ad))
                # axs[2].plot(psi_list_out, sx_list_out,'--', color = cl_dic[ad], label= 'outer {:.3E} W'.format(label_ad))
                axs[0].axvline(x= 1, color='black', lw=3, ls='--')
                axs[1].axvline(x= 1, color='black', lw=3, ls='--')
                # axs[2].axvline(x= 1, color='black', lw=3, ls='--')
                axs[0].add_artist(anchored_text_1)
                axs[1].add_artist(anchored_text_2)
                # axs[2].add_artist(anchored_text_3)
                axs[0].set_title('Heat flux scan with particle flux = {:.3E} (1/s)'.format(title_ap))
                axs[0].set_xlabel('$\psi_N$')
                axs[1].set_xlabel('$\psi_N$')
                # axs[2].set_xlabel('$\psi_N$')
                axs[0].legend(loc= 'lower right')
                    
                    
            else:
                print('neteTSplot_structure, please check the scan parameter')

    
    
    def twscan_targetnd_combine_method(self, log_flag, cl_dic, A_dic,
                                   scandetail, iterlist, scan_style):
        
        fig, axs = plt.subplots(1, 2)
            
        
        plt.subplots_adjust(hspace=.0)
            
        
        anchored_text_1 = AnchoredText('{}'.format('Source at target [$m^{-3}*s^{-1}$]'), loc='upper left')
        anchored_text_2 = AnchoredText('{}'.format('Atomic neutral density at target [$m^{-3}$]'), loc='upper left')
        
        
        for aa in iterlist:
            
            
            psi_dic = self.data['target_profile'][aa[0]][aa[1]]['psiN']
            sx_dic = self.data['target_profile'][aa[0]][aa[1]]['source']
            neuden_dic = self.data['target_profile'][aa[0]][aa[1]]['neuden']
            
            psi_list_in = psi_dic['inner target']
            
            sx_list_in = sx_dic['inner target']
            nd_list_in = neuden_dic['inner target']
            
            psi_list_out = psi_dic['outer target']
            
            sx_list_out = sx_dic['outer target']
            nd_list_out = neuden_dic['outer target']
            
            
            if self.DF.series_flag == 'twin_scan':
                
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
                
                    
                if log_flag:
                    axs[0].set_yscale('log')
                else:
                    pass
                
                
                
                axs[0].plot(psi_list_in, sx_list_in,'-', color = cl_dic[ad], label= 'HFS {:.3E} (1/s)'.format(label_ad))
                axs[1].plot(psi_list_in, nd_list_in,'-', color = cl_dic[ad], label= 'HFS {:.3E} (1/s)'.format(label_ad))
                axs[0].plot(psi_list_out, sx_list_out,'--', color = cl_dic[ad], label= 'LFS {:.3E} (1/s)'.format(label_ad))
                axs[1].plot(psi_list_out, nd_list_out,'--', color = cl_dic[ad], label= 'LFS {:.3E} (1/s)'.format(label_ad))
                axs[0].set_xlabel('$\psi_N$')
                axs[1].set_xlabel('$\psi_N$')
                axs[0].axvline(x= 1, color='black', lw=3, ls='--')
                axs[1].axvline(x= 1, color='black', lw=3, ls='--')
                axs[0].add_artist(anchored_text_1)
                axs[1].add_artist(anchored_text_2)
                axs[0].set_title('Particle flux scan with heat flux = {:.3E} W'.format(title_ap))
                
                axs[0].legend(loc= 'lower right')
                


            elif scan_style == 'tempscan':
                
                title_ap = float(ap)*pow(10, 20)
                label_ad = float(ad)*pow(10, 5)
                # exp_an_fit = fit_dat['exp_fit']
                # xcoord_cut = fit_dat['x_coord_cut']

                    
                if log_flag:
                    axs[0].set_yscale('log')
                else:
                    pass
            
                
                axs[0].plot(psi_list_in, sx_list_in,'-', color = cl_dic[ad], label= 'HFS {:.3E} W'.format(label_ad))
                axs[1].plot(psi_list_in, nd_list_in,'-', color = cl_dic[ad], label= 'HFS {:.3E} W'.format(label_ad))
                axs[0].plot(psi_list_out, sx_list_out,'--', color = cl_dic[ad], label= 'LFS {:.3E} W'.format(label_ad))
                axs[1].plot(psi_list_out, nd_list_out,'--', color = cl_dic[ad], label= 'LFS {:.3E} W'.format(label_ad))
                axs[0].axvline(x= 1, color='black', lw=3, ls='--')
                axs[1].axvline(x= 1, color='black', lw=3, ls='--')
                axs[0].add_artist(anchored_text_1)
                axs[1].add_artist(anchored_text_2)
                axs[0].set_title('Heat flux scan with particle flux = {:.3E} (1/s)'.format(title_ap))
                axs[0].set_xlabel('$\psi_N$')
                axs[1].set_xlabel('$\psi_N$')
                axs[0].legend(loc= 'lower right')
                    
                    
            else:
                print('neteTSplot_structure, please check the scan parameter')




                
    def twinscan_targetNT(self, scan_style, log_flag, match):
        
        
        if self.DF.withshift == False and self.DF.withseries == True:
            print('Opacity_study_radial_plot is not there yet, to be continue...')  
            
            
            if self.DF.series_flag == 'twin_scan':
                
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
                    
                    iter_key, color_dic, scan_title, label_dic = self.twa.twinscan_prep(ta = ta, 
                    keylist_b = keylist_b, scan_style = scan_style)
                    
                    
                    side_list = ['inner target', 'outer target']
                    
                    if match:
                        
                        self.twscan_targetNT_combine_method(iterlist = iter_key, cl_dic = color_dic,
                                    A_dic = label_dic, scan_style = scan_style, log_flag = log_flag,
                                    scandetail = scan_title)
                        
                        self.twscan_targetnd_combine_method(iterlist = iter_key, cl_dic = color_dic,
                                    A_dic = label_dic, scan_style = scan_style, log_flag = log_flag,
                                    scandetail = scan_title)
                        
                    else:
                        
                        for sk in side_list:
                            
                            self.twscan_targetNT_method(iterlist = iter_key, cl_dic = color_dic,
                                        A_dic = label_dic, scan_style = scan_style, log_flag = log_flag,
                                        scandetail = scan_title, side = sk)
                        
            
        else:
            print('Opacity_study_radial_plot has a bug')
        

                
        
        
"""


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


"""   
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        