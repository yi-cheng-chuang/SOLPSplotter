# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 13:00:30 2024

@author: ychuang
"""



from SOLPSplotter_iout_dataprocess import iout_process
from matplotlib.offsetbox import AnchoredText
import load_B2_data_method as lBdm
import numpy as np
import math
import matplotlib.pyplot as plt


class Diff_R_calc(iout_process):
    def __init__(self, DefaultSettings, loadDS):
        iout_process.__init__(self, DefaultSettings, loadDS)
    
    
    def set_plot(self):
        
        plt.rcParams.update({'font.weight': 'normal'})
        plt.rc('lines', linewidth= 5, markersize= 9)
        plt.rcParams.update({'font.size': 16})
        plt.rcParams.update({'figure.facecolor':'w'})
        plt.rcParams.update({'mathtext.default': 'regular'})
  

      

    def diff_quant_y(self, iout_dat):


        # del_x_Line = 
        del_y = np.zeros(np.shape(iout_dat))
        
        # print(np.shape(iout_dat))
        
        for ix in range(np.shape(iout_dat)[0]):
            for iy in range(np.shape(iout_dat)[1]):
                if (iy +1 == np.shape(iout_dat)[1]):
                    del_y[ix][iy] = iout_dat[ix][iy]
                
                else:
                    del_y[ix][iy] = iout_dat[ix][iy] - iout_dat[ix][iy + 1]
        
        return del_y


    def diff_quant_x(self, iout_dat):


        # del_x_Line = 
        del_x = np.zeros(np.shape(iout_dat))
        
        # print(np.shape(iout_dat)[1])
        
        for ix in range(np.shape(iout_dat)[0]):
            for iy in range(np.shape(iout_dat)[1]):
                if (ix + 1 == np.shape(iout_dat)[0]):
                    del_x[ix][iy] = iout_dat[ix][iy]
                    
                else:
                    del_x[ix][iy] = iout_dat[ix + 1][iy] - iout_dat[ix][iy]
        
        return del_x


    def sum_quant_x(self, iout_dat):


        # del_x_Line = 
        sum_x = np.zeros(np.shape(iout_dat))
        
        # print(np.shape(iout_dat))
        
        for ix in range(np.shape(iout_dat)[1]):
            for iy in range(np.shape(iout_dat)[0]):
                if (ix +1 == np.shape(iout_dat)[1]):
                    sum_x[iy][ix] = iout_dat[iy][ix]
                else:
                    sum_x[iy][ix] = iout_dat[iy][ix] + iout_dat[iy][ix+1]
        
        return sum_x

    
    def polNT(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                          'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            fig, axs = plt.subplots(2, 1)
            
            ne_text = AnchoredText('{}'.format('(a) Electron density [$m^{-3}$]'), 
                                         loc='upper center')
            
            te_text = AnchoredText('{}'.format('(b) Electron temperature [eV]'), 
                                         loc='upper center')
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                b2fstate = self.data['b2fstate'][aa]
                
                ne = b2fstate['ne'][1:97, 1:37]
                Te_J = b2fstate['te'][1:97, 1:37]
                
                ev = 1.6021766339999999 * pow(10, -19)
                te = Te_J / ev
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs[0].plot(ang_list, ne[st:ed, 17], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                axs[1].plot(ang_list, te[st:ed, 17], color = color_dic[aa])
            
            axs[0].legend(loc= 'upper left')
            axs[0].add_artist(ne_text)
            axs[1].add_artist(te_text)
            axs[1].set_xlabel('poloidal angle')
            plt.suptitle('$n_e$, $T_e$ at last cell before separatrix')
            plt.subplots_adjust(hspace=.0)
    
    
    
    def parallel_flow(self, pol_list, side):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                          'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            fig, axs = plt.subplots(1, 2)
            
            pol_text = AnchoredText('{}'.format('(a) Poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                         loc='upper center')
            
            S_text = AnchoredText('{}'.format('(b) source [$m^{-3}*s^{-1}$]'), 
                                         loc='upper center')
            
            neuden_text = AnchoredText('{}'.format('(c) neutral density [$m^{-3}$]'), 
                                         loc='upper center')
            
            # side = ['inner target', 'outer target']
            
            # side = 'inner target'
            # for kk in side:
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                
                fnnx = np.divide(fnaxs, vol)
                fnax = np.multiply(fnnx, hx)
                
                
                fnays = self.data['b2wdat'][aa]['b2npc_fnays'][0][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                
                fnny = np.divide(fnays, vol)
                fnay = np.multiply(fnny, hy)
                
                
                
                source = self.data['b2wdat'][aa]['b2npc_sna'][0][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                
                sterm = np.divide(source, vol)
                
                
                # ang_list = self.data['angle']['angle_list'][aa]
                # st = int(pol_list[0])
                # ed = int(pol_list[-1]) + 1
                
                
                if side == 'inner target':
                    psiN = self.data['psi']['psival'][aa][1:37, 1]
                    fnax_t = fnax[1, :]
                    fnay_t = fnay[1, :]
                    sterm_t = sterm[1, :]
                    
                    
                elif side == 'outer target':
                    
                    psiN = self.data['psi']['psival'][aa][1:37, -1]
                    fnax_t = fnax[-1, :]
                    fnay_t = fnay[-1, :]
                    sterm_t = sterm[-1, :]

                
                    
                axs[0].plot(psiN, fnax_t, '-', color = color_dic[aa], 
                              label = 'A = {}'.format(A_dic[aa]))
                axs[1].plot(psiN, sterm_t, '-', color = color_dic[aa])
                # axs[2].plot(psiN, fnay_t, '-', color = color_dic[aa])
                    
            
            axs[0].legend(loc= 'upper left')
            axs[1].legend(loc= 'upper left')
            axs[0].set_title('parallel flux {}'.format(side))
            axs[1].set_title('source {}'.format(side))
            # axs[2].set_title('radial flux')
            axs[0].set_xlabel('$\psi_N$')
            axs[1].set_xlabel('$\psi_N$')
            # axs[2].set_xlabel('$\psi_N$')
            # axs[0].set_yscale('log')
            # axs[1].set_yscale('log')
            plt.subplots_adjust(hspace=.0)
    
    
    def pllflow_evolve(self, pol_list, side):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                          'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            fig, axs = plt.subplots(1, 4, sharey= True)
            
            pol_text = AnchoredText('{}'.format('(a) Poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                         loc='upper center')
            
            S_text = AnchoredText('{}'.format('(b) source [$m^{-3}*s^{-1}$]'), 
                                         loc='upper center')
            
            neuden_text = AnchoredText('{}'.format('(c) neutral density [$m^{-3}$]'), 
                                         loc='upper center')
            
            # side = ['inner target', 'outer target']
            
            # side = 'inner target'
            # for kk in side:
                      
            for ii, aa in enumerate(self.data['dircomp']['multi_shift']):
                
                
                fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                
                fnnx = np.divide(fnaxs, vol)
                fnax = np.multiply(fnnx, hx)
                
                
                
                source = self.data['b2wdat'][aa]['b2npc_sna'][0][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                
                sterm = np.divide(source, vol)
                
                
                if side == 'HFS':
                    psiN = self.data['psi']['psival'][aa][1:37, 1]
                    fnax_t = fnax[0, :]
                    fnax_s = fnax[5, :]
                    fnax_m = fnax[15, :]
                    fnax_px = fnax[23, :]
                    fnax_p = fnax[25, :]
                    
                    axs[ii].plot(psiN, np.abs(fnax_t), '-', color = 'red', 
                                  label = 'inner target')
                    axs[ii].axhline(y= max(np.abs(fnax_t)), ls = '--', color = 'red')
                    axs[ii].plot(psiN, np.abs(fnax_s), '-', color = 'orange', 
                                  label = 'index 5')
                    
                    # axs[ii].axhline(y= max(fnax_s),ls = '--', color = 'orange')
                    axs[ii].plot(psiN, np.abs(fnax_m), '-', color = 'green', 
                                  label = 'index 15')
                    
                    # axs[ii].axhline(y= max(fnax_m), ls = '--', color = 'green')
                    axs[ii].plot(psiN, np.abs(fnax_px), '-', color = 'blue', 
                                  label = 'index 23, before left cut 24')
                    
                    # axs[ii].axhline(y= max(fnax_px), ls = '--', color = 'blue')
                    
                    axs[ii].plot(psiN, np.abs(fnax_p), '-', color = 'purple', 
                                  label = 'poloidal angle 250, index 25')
                    axs[ii].axhline(y= max(np.abs(fnax_p)), ls = '--', color = 'purple')
                    axs[ii].set_xlabel('$\psi_N$')
                    axs[ii].legend(loc= 'lower right')
                    axs[ii].set_title('A = {} {} poloidal flux evolution'.format(A_dic[aa], side))
                    axs[ii].set_yscale('log')
                    
                    
                elif side == 'LFS':
                    
                    psiN = self.data['psi']['psival'][aa][18:37, -1]
                    fnax_t = fnax[95, 17:]
                    fnax_s = fnax[90, 17:]
                    fnax_m = fnax[80, 17:]
                    fnax_px = fnax[73, 17:]
                    fnax_p = fnax[59, 17:]
                    
                    axs[ii].plot(psiN, fnax_t, '-', color = 'red', 
                                  label = 'outer target')
                    axs[ii].axhline(y= max(fnax_t), ls = '--', color = 'red')
                    axs[ii].plot(psiN, fnax_s, '-', color = 'orange', 
                                  label = 'index 90')
                    
                    # axs[ii].axhline(y= max(fnax_s), ls = '--', color = 'orange')
                    
                    axs[ii].plot(psiN, fnax_m, '-', color = 'green', 
                                  label = 'index 80')
                    
                    # axs[ii].axhline(y= max(fnax_m), ls = '--', color = 'green')
                    
                    axs[ii].plot(psiN, fnax_px, '-', color = 'blue', 
                                  label = 'index 73, after right cut 72')
                    
                    # axs[ii].axhline(y= max(fnax_px), ls = '--', color = 'blue')
                    axs[ii].plot(psiN, fnax_p, '-', color = 'purple', 
                                  label = 'inner midplane, index 59')
                    
                    axs[ii].axhline(y= max(fnax_p), ls = '--', color = 'purple')
                    axs[ii].set_xlabel('$\psi_N$')
                    axs[ii].legend(loc= 'lower right')
                    axs[ii].set_title('A = {} {} poloidal flux evolution'.format(A_dic[aa], side))
                    axs[ii].set_yscale('log')
                
                    
            plt.suptitle('poloidal flux {} evolution'.format(side))
            # axs[0].set_yscale('log')
            # axs[1].set_yscale('log')
            plt.subplots_adjust(wspace=.0)
    
    
    
    def tarNT_evolve(self, pol_list, side):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                          'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            fig, axs = plt.subplots(1, 4, sharey= True)
            
            pol_text = AnchoredText('{}'.format('(a) Poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                         loc='upper center')
            
            S_text = AnchoredText('{}'.format('(b) source [$m^{-3}*s^{-1}$]'), 
                                         loc='upper center')
            
            neuden_text = AnchoredText('{}'.format('(c) neutral density [$m^{-3}$]'), 
                                         loc='upper center')
            
            # side = ['inner target', 'outer target']
            
            # side = 'inner target'
            # for kk in side:
                      
            for ii, aa in enumerate(self.data['dircomp']['multi_shift']):
                
                
                "I will use the data size as small version"
                    
                b2fstate = self.data['b2fstate'][aa]
                
                ne_pro = b2fstate['ne'][1:97, 1:37]
                Te_J = b2fstate['te'][1:97, 1:37]
                    
                ev = 1.6021766339999999 * pow(10, -19)
                tepro = Te_J / ev
                
                
                if side == 'HFS':
                    psiN = self.data['psi']['psival'][aa][1:37, 1]
                    te_t = tepro[0, :]
                    te_s = tepro[5, :]
                    te_m = tepro[15, :]
                    te_px = tepro[23, :]
                    te_p = tepro[25, :]
                    
                    axs[ii].plot(psiN, np.abs(te_t), '-', color = 'red', 
                                  label = 'inner target')
                    axs[ii].axhline(y= max(np.abs(te_t)), ls = '--', color = 'red')
                    axs[ii].plot(psiN, np.abs(te_s), '-', color = 'orange', 
                                  label = 'index 5')
                    
                    # axs[ii].axhline(y= max(fnax_s),ls = '--', color = 'orange')
                    axs[ii].plot(psiN, np.abs(te_m), '-', color = 'green', 
                                  label = 'index 15')
                    
                    # axs[ii].axhline(y= max(fnax_m), ls = '--', color = 'green')
                    axs[ii].plot(psiN, np.abs(te_px), '-', color = 'blue', 
                                  label = 'index 23, before left cut 24')
                    
                    # axs[ii].axhline(y= max(fnax_px), ls = '--', color = 'blue')
                    
                    axs[ii].plot(psiN, np.abs(te_p), '-', color = 'purple', 
                                  label = 'poloidal angle 250, index 25')
                    axs[ii].axhline(y= max(np.abs(te_p)), ls = '--', color = 'purple')
                    axs[ii].set_xlabel('$\psi_N$')
                    axs[ii].legend(loc= 'lower right')
                    axs[ii].set_title('A = {} {} poloidal flux evolution'.format(A_dic[aa], side))
                    axs[ii].set_yscale('log')
                    
                    
                elif side == 'LFS':
                    
                    psiN = self.data['psi']['psival'][aa][18:37, -1]
                    te_t = tepro[95, 17:]
                    te_s = tepro[90, 17:]
                    te_m = tepro[80, 17:]
                    te_px = tepro[73, 17:]
                    te_p = tepro[59, 17:]
                    
                    axs[ii].plot(psiN, te_t, '-', color = 'red', 
                                  label = 'outer target')
                    axs[ii].axhline(y= max(te_t), ls = '--', color = 'red')
                    axs[ii].plot(psiN, te_s, '-', color = 'orange', 
                                  label = 'index 90')
                    
                    # axs[ii].axhline(y= max(te_s), ls = '--', color = 'orange')
                    
                    axs[ii].plot(psiN, te_m, '-', color = 'green', 
                                  label = 'index 80')
                    
                    # axs[ii].axhline(y= max(te_m), ls = '--', color = 'green')
                    
                    axs[ii].plot(psiN, te_px, '-', color = 'blue', 
                                  label = 'index 73, after right cut 72')
                    
                    # axs[ii].axhline(y= max(te_px), ls = '--', color = 'blue')
                    axs[ii].plot(psiN, te_p, '-', color = 'purple', 
                                  label = 'inner midplane, index 59')
                    
                    axs[ii].axhline(y= max(te_p), ls = '--', color = 'purple')
                    axs[ii].set_xlabel('$\psi_N$')
                    axs[ii].legend(loc= 'lower right')
                    axs[ii].set_title('A = {} {} poloidal flux evolution'.format(A_dic[aa], side))
                    axs[ii].set_yscale('log')
                
                    
            plt.suptitle('poloidal flux {} evolution'.format(side))
            # axs[0].set_yscale('log')
            # axs[1].set_yscale('log')
            plt.subplots_adjust(wspace=.0)
    
    

    

    def test_case(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
                     
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
                       

            fig, axs = plt.subplots()
            
            
            for aa in self.data['dircomp']['multi_shift']:
                
                
                na = self.data['b2wdat'][aa]['na'][0][1:97, 1:37]
                nadiff = self.diff_quant_y(iout_dat = na)
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                dny = np.divide(nadiff, hy)
                
                
                test_n = self.data['b2wdat'][aa]['b2tfnb_dPat_mdf_gradnay'][0][1:97, 1:37]
                g_term = self.data['b2wdat'][aa]['b2trno_cdnay'][0][1:97, 1:37]
                
                
                gca = self.data['b2wdat'][aa]['b2tqna_dna0'][0][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hy2 = np.multiply(hy, hy)
                geo = np.divide(vol, hy2)
                cgt = np.multiply(gca, geo)
                
                
                dndy = np.divide(test_n, cgt)
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs.plot(ang_list, nadiff[st:ed, 18],'--', color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                axs.plot(ang_list, dndy[st:ed, 18],'-o', color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                
            
            axs.legend(loc= 'best')
            axs.set_title('source')
        
        


    def polsource(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
                     
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
                       

            fig, axs = plt.subplots()
                
            for aa in self.data['dircomp']['multi_shift']:
                
                Sn = self.data['b2wdat'][aa]['b2npc_sna'][0]
                vol = self.data['b2wdat'][aa]['vol']
                
                sng = np.divide(Sn, vol)
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs.plot(ang_list, sng[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
            
            axs.legend(loc= 'best')
            axs.set_title('source term [$s^{-1}$]')
            axs.set_xlabel('poloidal angle')
                
    
    def hxhy(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            fig, axs = plt.subplots(2, 1)
            
            hx_text = AnchoredText('{}'.format('(a) $h_{\Theta}$'), 
                                         loc='upper center')
            
            hy_text = AnchoredText('{}'.format('(b) $h_r$'), 
                                         loc='upper center')
                      
            for aa in self.data['dircomp']['multi_shift']:
                

                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs[0].plot(ang_list, hx[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                axs[1].plot(ang_list, hy[st:ed, 18], color = color_dic[aa])
            
            
            axs[0].legend(loc= 'upper right')
            axs[0].add_artist(hx_text)
            axs[1].add_artist(hy_text)
            axs[1].set_xlabel('poloidal angle')

            plt.subplots_adjust(hspace=.0)
            
    
    
    def term_one(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                if aa == 'org':
                    pass
                else:
                    
                    "org case"
                    

                    fnaxs_org = self.data['b2wdat']['org']['b2npc_fnaxs'][0][1:97, 1:37]
                    vol_org = self.data['b2wdat']['org']['vol'][1:97, 1:37]
                    hx_org = self.data['b2wdat']['org']['hx'][1:97, 1:37]
                    hy_org = self.data['b2wdat']['org']['hy'][1:97, 1:37]
                    
                    
                    g_coe_org = np.divide(vol_org, hx_org)
                    
                    fnax_org = np.divide(fnaxs_org, g_coe_org)
                    cup_fnax_org = np.multiply(fnax_org, hy_org)
                    
                    
                    cpfnaxdiff_org = self.diff_quant_x(iout_dat = cup_fnax_org)
                    acp = np.divide(cpfnaxdiff_org, hx_org)
                    agcoe = np.multiply(hx_org, hy_org)
                    t1_org = np.divide(acp, agcoe)

                    
                    
                    
                    fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                    vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                    hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                    hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                    
                    
                    g_coe = np.divide(vol, hx)
                    
                    fnax = np.divide(fnaxs, g_coe)
                    cup_fnax_org = np.multiply(fnax, hy)
                    
                    
                    cpfnaxdiff_org = self.diff_quant_x(iout_dat = cup_fnax_org)
                    acp = np.divide(cpfnaxdiff_org, hx)
                    agcoe = np.multiply(hx_org, hy_org)
                    t1 = np.divide(acp, agcoe)
                    
                    del_t1 = t1 - t1_org
                    
                    
                    ang_list = self.data['angle']['angle_list'][aa]
                    # print(np.shape(nadiff))
                    st = int(pol_list[0])
                    ed = int(pol_list[-1]) + 1
                    
                    
                    # axs.plot(ang_list, dndx[st:ed, 1])
                    axs.plot(ang_list, del_t1[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                    
                
                
            
            axs.legend(loc= 'best')
            axs.set_title('term 1')
            axs.set_xlabel('poloidal angle')
    
    
           
    def term_three(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                if aa == 'org':
                    pass
                else:
                    
                    "org case"
                    
                    Radloc_org = np.transpose(self.data['grid']['RadLoc']['org'])
                    
                    Rdiff_org = self.diff_quant_y(iout_dat = Radloc_org)[1:97, 1:37]
                    hy_org = self.data['b2wdat']['org']['hy'][1:97, 1:37]
                    
                    
                    fnays_org = self.data['b2wdat']['org']['b2npc_fnays'][0][1:97, 1:37]
                    vol_org = self.data['b2wdat']['org']['vol'][1:97, 1:37]
                                       
                    fnay_coe_org = np.divide(fnays_org, vol_org)                   
                    
                    Rdy_org = np.divide(Rdiff_org, hy_org)
                    Rloc_org = Radloc_org[1:97, 1:37]
                    RRdy_org = np.divide(Rdy_org, Rloc_org)
                    
                    t3_org = np.multiply(fnay_coe_org, RRdy_org)
                    
                    
                    
                    Radloc = np.transpose(self.data['grid']['RadLoc'][aa])
                    
                    Rdiff = self.diff_quant_y(iout_dat = Radloc)[1:97, 1:37]
                    hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                    
                    Rdy = np.divide(Rdiff, hy)
                    Rloc = Radloc[1:97, 1:37]
                    RRdy = np.divide(Rdy, Rloc)
                    
                    
                    fnays = self.data['b2wdat'][aa]['b2npc_fnays'][0][1:97, 1:37]
                    vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                                       
                    fnay_coe = np.divide(fnays, vol)
                    
                    t3 = np.multiply(fnay_coe, RRdy)
                    
                    del_t3 = t3 - t3_org
                    
                    
                    ang_list = self.data['angle']['angle_list'][aa]
                    # print(np.shape(nadiff))
                    st = int(pol_list[0])
                    ed = int(pol_list[-1]) + 1
                    
                    
                    # axs.plot(ang_list, dndx[st:ed, 1])
                    axs.plot(ang_list, del_t3[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                    
                
                
            
            axs.legend(loc= 'best')
            axs.set_title('term 3')
            axs.set_xlabel('poloidal angle')
    
    
    
    
    def term_two(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                if aa == 'org':
                    pass
                else:
                    
                    "org case"
                    
                    Radloc_org = np.transpose(self.data['grid']['RadLoc']['org'])
                    
                    Rdiff_org = self.diff_quant_x(iout_dat = Radloc_org)[1:97, 1:37]
                    hx_org = self.data['b2wdat']['org']['hx'][1:97, 1:37]
                    
                    
                    fnaxs_org = self.data['b2wdat']['org']['b2npc_fnaxs'][0][1:97, 1:37]
                    vol_org = self.data['b2wdat']['org']['vol'][1:97, 1:37]
                                       
                    fnax_coe_org = np.divide(fnaxs_org, vol_org)                   
                    
                    Rdx_org = np.divide(Rdiff_org, hx_org)
                    Rloc_org = Radloc_org[1:97, 1:37]
                    RRdx_org = np.divide(Rdx_org, Rloc_org)
                    
                    t2_org = np.multiply(fnax_coe_org, RRdx_org)
                    
                    
                    
                    Radloc = np.transpose(self.data['grid']['RadLoc'][aa])
                    
                    Rdiff = self.diff_quant_x(iout_dat = Radloc)[1:97, 1:37]
                    hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                    
                    Rdx = np.divide(Rdiff, hx)
                    Rloc = Radloc[1:97, 1:37]
                    RRdx = np.divide(Rdx, Rloc)
                    
                    
                    fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                    vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                                       
                    fnax_coe = np.divide(fnaxs, vol)
                    
                    t2 = np.multiply(fnax_coe, RRdx)
                    
                    del_t2 = t2 - t2_org
                    
                    
                    ang_list = self.data['angle']['angle_list'][aa]
                    # print(np.shape(nadiff))
                    st = int(pol_list[0])
                    ed = int(pol_list[-1]) + 1
                    
                    
                    # axs.plot(ang_list, dndx[st:ed, 1])
                    axs.plot(ang_list, del_t2[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                    
                
                
            
            axs.legend(loc= 'best')
            axs.set_title('term 2')
            axs.set_xlabel('poloidal angle')
    

    
    def R_change(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                Radloc = np.transpose(self.data['grid']['RadLoc'][aa])
                
                Rdiff = self.diff_quant_x(iout_dat = Radloc)[1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                Rdx = np.divide(Rdiff, hx)
                Rloc = Radloc[1:97, 1:37]
                RRdx = np.divide(Rdx, Rloc)
                
                RnRdx = np.multiply(Rdx, Rloc)
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs.plot(ang_list, Rdiff[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
            
            axs.legend(loc= 'best')
            axs.set_title('(1/R)(dR/d$\Theta$)')
            axs.set_xlabel('poloidal angle')
            
            
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                Radloc = np.transpose(self.data['grid']['RadLoc'][aa])
                
                Rdiff = self.diff_quant_y(iout_dat = Radloc)[1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                Rdy = np.divide(Rdiff, hy)
                Rloc = Radloc[1:97, 1:37]
                RRdy = np.divide(Rdy, Rloc)
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs.plot(ang_list, Rdiff[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
            
            axs.legend(loc= 'best')
            axs.set_title('(1/R)(dR/dr)')
            axs.set_xlabel('poloidal angle')
    
    
    
    def sqrtg_change(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                vol = self.data['b2wdat'][aa]['vol']
                
                gdiff = self.diff_quant_x(iout_dat = vol)[1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                gdx = np.divide(gdiff, hx)
                vloc = vol[1:97, 1:37]
                Gdx = np.divide(gdx, vloc)
                
                # RnRdx = np.multiply(Rdx, Rloc)
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs.plot(ang_list, Gdx[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
            
            axs.legend(loc= 'best')
            axs.set_title('(1/$\sqrt{g}$)(d $\sqrt{g}$/dx)')
            
            
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                
                vol = self.data['b2wdat'][aa]['vol']
                
                gdiff = self.diff_quant_y(iout_dat = vol)[1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                vloc = vol[1:97, 1:37]
                gdy = np.divide(gdiff, hy)
                Gdy = np.divide(gdy, vloc)
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs.plot(ang_list, Gdy[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
            
            axs.legend(loc= 'best')
            axs.set_title('(1/$\sqrt{g}$)(d $\sqrt{g}$/dy)')
    
    
    
    

    def mag_constant(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                Radloc = np.transpose(self.data['grid']['RadLoc'][aa])
                
                bx = self.data['b2wdat'][aa]['bbx'][1:97, 1:37]
                bxb = self.data['b2wdat'][aa]['bx'][1:97, 1:37]
                Rloc = Radloc[1:97, 1:37]
                
                try_c = np.multiply(Rloc, bx)
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                axs.plot(ang_list, try_c[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                # axs.plot(ang_list, bxb[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                
            
            axs.set_title('$B_{pol}* R$')
            axs.legend(loc= 'best')




  
    def pol_psiN(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            fig, axs = plt.subplots()
            
            
            for aa in self.data['dircomp']['multi_shift']:
                

                    
                psiN = self.data['psi']['psival'][aa][1:37, 1:97]
                

                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                    
                axs.plot(ang_list, psiN[18, st:ed], '-', color = color_dic[aa], 
                             label = '{}'.format(A_dic[aa]))
            
            axs.legend(loc= 'best')
            axs.set_title('psiN')
            
            
            fig, axs = plt.subplots()
            
            
            for aa in self.data['dircomp']['multi_shift']:
                                
                psiN = self.data['psi']['psival'][aa][1:37, 1:97]
                
    
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                weight = self.data['polpsi_weight'][aa]
                
                psi_list = []
                
                for ii, kk in enumerate(pol_list):
                    tt = psiN[20, int(kk)]*weight[ii] + (1 - weight[ii])*psiN[18, int(kk)]
                    psi_list.append(tt)
                    
                    
                axs.plot(ang_list, psi_list, '-', color = color_dic[aa], 
                         label = '{}'.format(A_dic[aa]))
            
            
            axs.set_title('psiN is 1')
            
            
            fig, axs = plt.subplots()
            
            for aa in self.data['dircomp']['multi_shift']:
                                
                psiN = self.data['psi']['psival'][aa][1:37, 1:97]
                
    
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                # weight = self.data['polpsi_weight'][aa]
                
                # psi_list = []
                
                # for ii, kk in enumerate(pol_list):
                #     tt = psiN[20, int(kk)]*weight[ii] + (1 - weight[ii])*psiN[17, int(kk)]
                #     psi_list.append(tt)
                    
                    
                axs.plot(ang_list, psiN[20, st:ed], '-', color = color_dic[aa], 
                         label = '{}'.format(A_dic[aa]))
                
                axs.plot(ang_list, psiN[18, st:ed], '-', color = color_dic[aa], 
                         label = '{}'.format(A_dic[aa]))
            
            
            axs.set_title('compare psiN')
    

    def paper_fluxes(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            
            fig, axs = plt.subplots(2, 1)
            
            pol_text = AnchoredText('{}'.format('(a) Poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                         loc='upper center')
            
            rad_text = AnchoredText('{}'.format('(b) Radial flux $\Gamma_r$ [$m^{-2} s^{-1}$]'), 
                                         loc='upper center')
            
            
            axs[0].axhline(y=0, color = 'black', linestyle = '--', label= '$\Gamma_{\Theta}$ = 0')
            
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                
                fnn = np.divide(fnaxs, vol)
                fnax = np.multiply(fnn, hx)
                
                fnays = self.data['b2wdat'][aa]['b2npc_fnays'][0][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                
                fny = np.divide(fnays, vol)
                fnay = np.multiply(fny, hy)
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                psi_avg = False
                
                if psi_avg == True:
                
                    weight = self.data['polpsi_weight'][aa]
                    
                    fx_list = []
                    fy_list = []
                    
                    for ii, kk in enumerate(pol_list):
                        fx = fnax[int(kk), 20]*weight[ii] + (1 - weight[ii])*fnax[int(kk), 18]
                        fx_list.append(fx)
                        
                        fy = fnay[int(kk), 20]*weight[ii] + (1 - weight[ii])*fnay[int(kk), 18]
                        fy_list.append(fy)
                    
                    
                    axs[0].plot(ang_list, fx_list, color = color_dic[aa])
                    axs[1].plot(ang_list, fy_list, color = color_dic[aa], 
                                 label = '{}'.format(A_dic[aa]))
                
                else:
                    
                    axs[0].plot(ang_list, fnax[st:ed, 18], color = color_dic[aa])
                    axs[1].plot(ang_list, fnay[st:ed, 18], color = color_dic[aa], 
                                  label = 'A = {}'.format(A_dic[aa]))
                    
                

                

            
            axs[1].legend(loc= 'upper right')
            axs[0].add_artist(pol_text)
            axs[1].add_artist(rad_text)
            axs[0].legend(loc = 'upper right')
            axs[1].set_xlabel('poloidal angle')

            plt.subplots_adjust(hspace=.0)
            
            
            
            
            fig, axs = plt.subplots(2, 1)
            
            pol_text = AnchoredText('{}'.format('(a) Poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                         loc='upper center')
            
            rad_text = AnchoredText('{}'.format('(b) Radial flux $\Gamma_r$ [$m^{-2} s^{-1}$]'), 
                                         loc='upper center')
            
            for aa in self.data['dircomp']['multi_shift']:
                
                fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                
                fnn = np.divide(fnaxs, vol)
                fnax = np.multiply(fnn, hx)
                
                fnays = self.data['b2wdat'][aa]['b2npc_fnays'][0][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                
                fny = np.divide(fnays, vol)
                fnay = np.multiply(fny, hy)
                
                
                index_list = np.linspace(1, 24, 24)
                
                
                psiN = self.data['psi']['psival'][aa][1:37, 1:97]
                
                
                high_psi = psiN[20, :24]
                low_psi = psiN[18, :24]
                
                list_len = len(high_psi)
            
                weight = np.zeros(list_len)
                for x in range(list_len):
                    weight[x] = (1 - low_psi[x])/ (high_psi[x] - low_psi[x])
                
                
                fx_list = []
                fy_list = []
                
                for ii, kk in enumerate(index_list):
                    fx = fnax[int(kk) -1, 20]*weight[ii] + (1 - weight[ii])*fnax[int(kk) -1, 18]
                    fx_list.append(fx)
                    
                    fy = fnay[int(kk) -1, 20]*weight[ii] + (1 - weight[ii])*fnay[int(kk) - 1, 18]
                    fy_list.append(fy)
                
                
                axs[0].plot(index_list, fx_list, color = color_dic[aa])
                axs[1].plot(index_list, fy_list, color = color_dic[aa], 
                             label = '{}'.format(A_dic[aa]))
                
                

                # axs[0].plot(index_list, fnax[:24, 18], color = color_dic[aa])
                # axs[1].plot(index_list, fnay[:24, 18], color = color_dic[aa], 
                #              label = '{}'.format(A_dic[aa]))

            
            axs[1].legend(loc= 'upper right')
            axs[0].add_artist(pol_text)
            axs[1].add_artist(rad_text)
            axs[0].legend(loc = 'upper right')
            axs[1].set_xlabel('poloidal angle')
            plt.suptitle('inner leg')

            plt.subplots_adjust(hspace=.0)
            
            
            
            fig, axs = plt.subplots(2, 1)
            
            pol_text = AnchoredText('{}'.format('(a) Poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                         loc='upper center')
            
            rad_text = AnchoredText('{}'.format('(b) Radial flux $\Gamma_r$ [$m^{-2} s^{-1}$]'), 
                                         loc='upper center')
            
            for aa in self.data['dircomp']['multi_shift']:
                
                fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                
                fnn = np.divide(fnaxs, vol)
                fnax = np.multiply(fnn, hx)
                
                fnays = self.data['b2wdat'][aa]['b2npc_fnays'][0][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                
                fny = np.divide(fnays, vol)
                fnay = np.multiply(fny, hy)
                
                
                index_list = np.linspace(73, 96, 24)
                
                
                psiN = self.data['psi']['psival'][aa][1:37, 1:97]
                
                
                high_psi = psiN[20, 72:]
                low_psi = psiN[18, 72:]
                
                list_len = len(high_psi)
            
                weight = np.zeros(list_len)
                for x in range(list_len):
                    weight[x] = (1 - low_psi[x])/ (high_psi[x] - low_psi[x])
                
                
                fx_list = []
                fy_list = []
                
                for ii, kk in enumerate(index_list):
                    fx = fnax[int(kk) - 1, 20]*weight[ii] + (1 - weight[ii])*fnax[int(kk) -1, 18]
                    fx_list.append(fx)
                    
                    fy = fnay[int(kk) - 1, 20]*weight[ii] + (1 - weight[ii])*fnay[int(kk) -1, 18]
                    fy_list.append(fy)
                
                
                axs[0].plot(index_list, fx_list, color = color_dic[aa])
                axs[1].plot(index_list, fy_list, color = color_dic[aa], 
                             label = '{}'.format(A_dic[aa]))
                

                # axs[0].plot(index_list, fnax[72:, 18], color = color_dic[aa])
                # axs[1].plot(index_list, fnay[72:, 18], color = color_dic[aa], 
                #              label = '{}'.format(A_dic[aa]))

            
            axs[1].legend(loc= 'upper right')
            axs[0].add_artist(pol_text)
            axs[1].add_artist(rad_text)
            axs[0].legend(loc = 'upper right')
            axs[1].set_xlabel('poloidal angle')
            plt.suptitle('inner leg')

            plt.subplots_adjust(hspace=.0)
            
            
    
    
    def fluxes_no_geo(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]

                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                axs.plot(ang_list, fnaxs[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                axs.set_title('parallel flux')
            
            axs.legend(loc= 'best')
            
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fnays = self.data['b2wdat'][aa]['b2npc_fnays'][0][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                
                fny = np.divide(fnays, vol)
                fnay = np.multiply(fny, hy)
                
                share_c = np.multiply(fnay, hy)
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs.plot(ang_list, share_c[st:ed, 0], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
            
            axs.legend(loc= 'best')
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                ua = self.data['b2wdat'][aa]['ua'][0][1:97, 1:37]
                
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs.plot(ang_list, ua[st:ed, 18], color = color_dic[aa], 
                             label = '{}'.format(A_dic[aa]))
                axs.set_title('parallel velocity')
            
            axs.legend(loc= 'best')
    
    
    
    def continuity_terms(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                fxdiff = self.diff_quant_x(iout_dat = fnaxs)[1:97, 1:37]
                
                fnx = np.divide(fxdiff, hx)
                fnax = np.divide(fnx, vol)
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs.plot(ang_list, fnax[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                axs.set_title('continity equation parallel particle flux term')
            
            axs.legend(loc= 'best')
            # axs.set_yscale('log')
        
        
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fnays = self.data['b2wdat'][aa]['b2npc_fnays'][0]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                fydiff = self.diff_quant_y(iout_dat = fnays)[1:97, 1:37]
                
                fny = np.divide(fydiff, hy)
                fnay = np.divide(fny, vol)
                
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs.plot(ang_list, fnay[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                axs.set_title('continity equation radial particle flux term')
                # axs.set_yscale('log')
            
            axs.legend(loc= 'best')
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fnays = self.data['b2wdat'][aa]['b2npc_fnays'][0]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                fydiff = self.diff_quant_y(iout_dat = fnays)[1:97, 1:37]
                
                fny = np.divide(fydiff, hy)
                fnay = np.divide(fny, vol)
                
                fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                fxdiff = self.diff_quant_x(iout_dat = fnaxs)[1:97, 1:37]
                
                fnx = np.divide(fxdiff, hx)
                fnax = np.divide(fnx, vol)
                
                totcon = fnay + fnax
                
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs.plot(ang_list, totcon[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                axs.set_title('continity equation total')
                # axs.set_yscale('log')
            
            axs.legend(loc= 'best')

    
    def poloidal_transcoe(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                d_withcoe = self.data['b2wdat'][aa]['b2trno_cdnay'][0][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                
                hy2 = np.multiply(hy, hy)
                geo = np.divide(vol, hy2)
                Dcoe = np.divide(d_withcoe, geo)
                
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                # weight = self.data['polpsi_weight'][aa]
                
                # D_list = []
                
                # for ii, kk in enumerate(pol_list):
                #     tt = Dcoe[int(kk), 20]*weight[ii] + (1 - weight[ii])*Dcoe[int(kk), 18]
                #     D_list.append(tt)
                
                
                
                # axs.plot(ang_list, D_list, color = color_dic[aa], 
                #               label = '{}'.format(A_dic[aa]))
                
                
                axs.plot(ang_list, Dcoe[st:ed, 17], color = color_dic[aa], 
                              label = 'A = {}'.format(A_dic[aa]))
            
            axs.set_title('Particle diffusion coefficient y at last cell before separatrix')
            axs.set_xlabel('poloidal angle')
            axs.legend(loc= 'best')
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                d_withcoe = self.data['b2wdat'][aa]['b2trno_cdnax'][0][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                
                hx2 = np.multiply(hx, hx)
                geo = np.divide(vol, hx2)
                Dcoe = np.divide(d_withcoe, geo)
                
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                # weight = self.data['polpsi_weight'][aa]
                
                # D_list = []
                
                # for ii, kk in enumerate(pol_list):
                #     tt = Dcoe[int(kk), 20]*weight[ii] + (1 - weight[ii])*Dcoe[int(kk), 18]
                #     D_list.append(tt)
                
                
                
                # axs.plot(ang_list, D_list, color = color_dic[aa], 
                #               label = '{}'.format(A_dic[aa]))
                
                
                axs.plot(ang_list, Dcoe[st:ed, 17], color = color_dic[aa], 
                              label = 'A = {}'.format(A_dic[aa]))
            
            axs.set_title('Particle diffusion coefficient x at last cell before separatrix')
            axs.set_xlabel('poloidal angle')
            axs.legend(loc= 'best')
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                Dcoe = self.data['b2wdat'][aa]['b2tqna_dna0'][0][1:97, 1:37]
                
                # vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                # hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                
                # hx2 = np.multiply(hx, hx)
                # geo = np.divide(vol, hx2)
                # Dcoe = np.divide(d_withcoe, geo)
                
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                # weight = self.data['polpsi_weight'][aa]
                
                # D_list = []
                
                # for ii, kk in enumerate(pol_list):
                #     tt = Dcoe[int(kk), 20]*weight[ii] + (1 - weight[ii])*Dcoe[int(kk), 18]
                #     D_list.append(tt)
                
                
                
                # axs.plot(ang_list, D_list, color = color_dic[aa], 
                #               label = '{}'.format(A_dic[aa]))
                
                
                axs.plot(ang_list, Dcoe[st:ed, 17], color = color_dic[aa], 
                              label = 'A = {}'.format(A_dic[aa]))
            
            axs.set_title('Particle diffusion coefficient [$m^{2}/s$] at last cell before separatrix')
            axs.set_xlabel('poloidal angle')
            axs.legend(loc= 'best')



    def heat_fluxes(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8'}
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fheys = self.data['b2wdat'][aa]['b2nph9_fhey'][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                # fheydiff = self.diff_quant_x(iout_dat = fheys)
                
                fhy = np.multiply(fheys, hy)
                fhey = np.divide(fhy, vol)
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                weight = self.data['polpsi_weight'][aa]
                
                
                fy_list = []
                
                for ii, kk in enumerate(pol_list):
                    
                    fy = fhey[int(kk), 20]*weight[ii] + (1 - weight[ii])*fhey[int(kk), 18]
                    fy_list.append(fy)
                
                
                axs.plot(ang_list, fy_list, color = color_dic[aa], 
                             label = 'A = {}'.format(A_dic[aa]))
                
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                # axs.plot(ang_list, fhey[st:ed, 18], color = color_dic[aa], 
                #              label = '{}'.format(A_dic[aa]))
                axs.set_title('radial electron heat flux at separatrix[W]')
            
            axs.legend(loc= 'best')
            axs.set_xlabel('poloidal angle')
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fhexs = self.data['b2wdat'][aa]['b2nph9_fhex'][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                # fheydiff = self.diff_quant_x(iout_dat = fheys)
                
                fhx = np.multiply(fhexs, hx)
                fhex = np.divide(fhx, vol)
                
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                
                weight = self.data['polpsi_weight'][aa]
                
                
                fx_list = []
                
                for ii, kk in enumerate(pol_list):
                    
                    fx = fhex[int(kk), 20]*weight[ii] + (1 - weight[ii])*fhex[int(kk), 18]
                    fx_list.append(fx)
                
                
                axs.plot(ang_list, fx_list, color = color_dic[aa], 
                              label = 'A = {}'.format(A_dic[aa]))
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                # axs.plot(ang_list, fhex[st:ed, 18], color = color_dic[aa], 
                #               label = '{}'.format(A_dic[aa]))
                
                
                axs.set_title('poloidal electron heat flux at separatrix [W]')
            
            axs.legend(loc= 'best')
            axs.set_xlabel('poloidal angle')
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fhixs = self.data['b2wdat'][aa]['b2nph9_fhix'][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                # fheydiff = self.diff_quant_x(iout_dat = fheys)
                
                fhx = np.multiply(fhixs, hx)
                fhix = np.divide(fhx, vol)
                
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                
                weight = self.data['polpsi_weight'][aa]
                
                
                fx_list = []
                
                for ii, kk in enumerate(pol_list):
                    
                    fx = fhix[int(kk), 20]*weight[ii] + (1 - weight[ii])*fhix[int(kk), 18]
                    fx_list.append(fx)
                
                
                axs.plot(ang_list, fx_list, color = color_dic[aa], 
                              label = 'A = {}'.format(A_dic[aa]))
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                # axs.plot(ang_list, fhex[st:ed, 18], color = color_dic[aa], 
                #               label = '{}'.format(A_dic[aa]))
                
                
                axs.set_title('poloidal ion heat flux at separatrix [W]')
            
            axs.legend(loc= 'best')
            axs.set_xlabel('poloidal angle')
            
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fheys = self.data['b2wdat'][aa]['b2nph9_fhey'][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                # fheydiff = self.diff_quant_x(iout_dat = fheys)
                
                fhy = np.multiply(fheys, hy)
                fhey = np.divide(fhy, vol)
                
                
                index_list = np.linspace(1, 23, 23)
                
                
                psiN = self.data['psi']['psival'][aa][1:37, 1:97]
                
                
                high_psi = psiN[20, :24]
                low_psi = psiN[18, :24]
                
                list_len = len(high_psi)
            
                weight = np.zeros(list_len)
                for x in range(list_len):
                    weight[x] = (1 - low_psi[x])/ (high_psi[x] - low_psi[x])
                
                
                fy_list = []
                
                for ii, kk in enumerate(index_list):
                    
                    fy = fhey[int(kk) -1, 20]*weight[ii] + (1 - weight[ii])*fhey[int(kk) - 1, 18]
                    fy_list.append(fy)
                
                
                axs.plot(index_list, fy_list, color = color_dic[aa], 
                             label = '{}'.format(A_dic[aa]))
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                # axs.plot(index_list, fhey[1:24, 18], color = color_dic[aa], 
                #              label = '{}'.format(A_dic[aa]))
                
                axs.set_title('radial electron heat flux inner leg')
            
            axs.legend(loc= 'best')
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fhexs = self.data['b2wdat'][aa]['b2nph9_fhex'][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                # fheydiff = self.diff_quant_x(iout_dat = fheys)
                
                fhx = np.multiply(fhexs, hx)
                fhex = np.divide(fhx, vol)
                
                
                
                index_list = np.linspace(1, 23, 23)
                
                psiN = self.data['psi']['psival'][aa][1:37, 1:97]
                
                
                high_psi = psiN[20, :24]
                low_psi = psiN[18, :24]
                
                list_len = len(high_psi)
            
                weight = np.zeros(list_len)
                for x in range(list_len):
                    weight[x] = (1 - low_psi[x])/ (high_psi[x] - low_psi[x])
                
                
                fx_list = []
                
                for ii, kk in enumerate(index_list):
                    fx = fhex[int(kk) -1, 20]*weight[ii] + (1 - weight[ii])*fhex[int(kk) -1, 18]
                    fx_list.append(fx)
                
                
                axs.plot(index_list, fx_list, color = color_dic[aa], 
                             label = '{}'.format(A_dic[aa]))
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                # axs.plot(index_list, fhex[1:24, 18], color = color_dic[aa], 
                #              label = '{}'.format(A_dic[aa]))
                axs.set_title('poloidal electron heat flux inner leg')
            
            axs.legend(loc= 'best')
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fheys = self.data['b2wdat'][aa]['b2nph9_fhey'][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                # fheydiff = self.diff_quant_x(iout_dat = fheys)
                
                fhy = np.multiply(fheys, hy)
                fhey = np.divide(fhy, vol)
                
                
                index_list = np.linspace(73, 96, 24)
                
                
                psiN = self.data['psi']['psival'][aa][1:37, 1:97]
                
                
                high_psi = psiN[20, 72:]
                low_psi = psiN[18, 72:]
                
                list_len = len(high_psi)
            
                weight = np.zeros(list_len)
                for x in range(list_len):
                    weight[x] = (1 - low_psi[x])/ (high_psi[x] - low_psi[x])
                
                
                fy_list = []
                
                for ii, kk in enumerate(index_list):
                    
                    fy = fhey[int(kk) -1, 20]*weight[ii] + (1 - weight[ii])*fhey[int(kk) - 1, 18]
                    fy_list.append(fy)
                
                
                axs.plot(index_list, fy_list, color = color_dic[aa], 
                             label = '{}'.format(A_dic[aa]))
                
                
                
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                # axs.plot(index_list, fhey[73:, 18], color = color_dic[aa], 
                #              label = '{}'.format(A_dic[aa]))
                axs.set_title('radial electron heat flux outer leg')
            
            axs.legend(loc= 'best')
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                fhexs = self.data['b2wdat'][aa]['b2nph9_fhex'][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                # fheydiff = self.diff_quant_x(iout_dat = fheys)
                
                fhx = np.multiply(fhexs, hx)
                fhex = np.divide(fhx, vol)
                
                
                index_list = np.linspace(73, 96, 24)
                
                
                psiN = self.data['psi']['psival'][aa][1:37, 1:97]
                
                
                high_psi = psiN[20, 72:]
                low_psi = psiN[18, 72:]
                
                list_len = len(high_psi)
            
                weight = np.zeros(list_len)
                for x in range(list_len):
                    weight[x] = (1 - low_psi[x])/ (high_psi[x] - low_psi[x])
                
                
                fx_list = []
                
                for ii, kk in enumerate(index_list):
                    
                    fx = fhex[int(kk) -1, 20]*weight[ii] + (1 - weight[ii])*fhex[int(kk) - 1, 18]
                    fx_list.append(fx)
                
                
                axs.plot(index_list, fx_list, color = color_dic[aa], 
                             label = '{}'.format(A_dic[aa]))
                
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                # axs.plot(index_list, fhex[73:, 18], color = color_dic[aa], 
                #              label = '{}'.format(A_dic[aa]))
                
                
                axs.set_title('poloidal electron heat flux outer leg')
            
            axs.legend(loc= 'best')
            
            
            
    def heatS_no_geo(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                heat_source = self.data['b2wdat'][aa]['b2nph9_she'][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                
                heatS = np.divide(heat_source, vol)

                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                weight = self.data['polpsi_weight'][aa]
                
                
                hS_list = []
                
                for ii, kk in enumerate(pol_list):
                    
                    hs = heatS[int(kk), 20]*weight[ii] + (1 - weight[ii])*heatS[int(kk), 18]
                    hS_list.append(hs)
                
                
                axs.plot(ang_list, hS_list, color = color_dic[aa], 
                             label = 'A = {}'.format(A_dic[aa]))
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                # axs.plot(ang_list, fnaxs[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                axs.set_title('heat source')
            
            axs.legend(loc= 'best')
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                heat_source = self.data['b2wdat'][aa]['b2nph9_she'][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                
                heatS = np.divide(heat_source, vol)

                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                index_list = np.linspace(73, 96, 24)
                
                
                psiN = self.data['psi']['psival'][aa][1:37, 1:97]
                
                
                high_psi = psiN[20, 72:]
                low_psi = psiN[18, 72:]
                
                list_len = len(high_psi)
            
                weight = np.zeros(list_len)
                for x in range(list_len):
                    weight[x] = (1 - low_psi[x])/ (high_psi[x] - low_psi[x])
                
                
                hS_list = []
                
                for ii, kk in enumerate(index_list):
                    
                    hs = heatS[int(kk) -1, 20]*weight[ii] + (1 - weight[ii])*heatS[int(kk) -1, 18]
                    hS_list.append(hs)
                
                
                axs.plot(index_list, hS_list, color = color_dic[aa], 
                             label = 'A = {}'.format(A_dic[aa]))
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                # axs.plot(ang_list, fnaxs[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                axs.set_title('heat source outer leg')
            
            axs.legend(loc= 'best')
            
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                heat_source = self.data['b2wdat'][aa]['b2nph9_she'][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                
                heatS = np.divide(heat_source, vol)

                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                index_list = np.linspace(1, 23, 23)
                
                
                psiN = self.data['psi']['psival'][aa][1:37, 1:97]
                
                
                high_psi = psiN[20, :24]
                low_psi = psiN[18, :24]
                
                list_len = len(high_psi)
            
                weight = np.zeros(list_len)
                for x in range(list_len):
                    weight[x] = (1 - low_psi[x])/ (high_psi[x] - low_psi[x])
                
                
                hS_list = []
                
                for ii, kk in enumerate(index_list):
                    
                    hs = heatS[int(kk) -1, 20]*weight[ii] + (1 - weight[ii])*heatS[int(kk) -1, 18]
                    hS_list.append(hs)
                
                
                axs.plot(index_list, hS_list, color = color_dic[aa], 
                             label = 'A = {}'.format(A_dic[aa]))
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                # axs.plot(ang_list, fnaxs[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                axs.set_title('heat source inner leg')
            
            axs.legend(loc= 'best')
    
    
    
    
    def check_targetNT(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
                     
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
                       
            div_side_list = ['inner target', 'outer target', 
                             'inner x point boundary', 'outer x point boundary']
            # qu = self.load_iout(filename = dname, simple_quant = quant)
            
            side = ['inner target', 'outer target']
            
            
            for kk in side:
                
                fig, axs = plt.subplots()
                
                
                for aa in self.data['dircomp']['multi_shift']:
                    
                    b2fstate = self.data['b2fstate'][aa]
                    
                    ne_pro = b2fstate['ne'].transpose()
                    
                    if kk == 'inner target':
                        psiN = self.data['psi']['psival'][aa][1:37, 1]
                        plot_data = ne_pro[1:37, 1]
                        
                    elif kk == 'outer target':
                        
                        psiN = self.data['psi']['psival'][aa][1:37, -2]
                        plot_data = ne_pro[1:37, -1]

                    
                        
                    axs.plot(psiN, plot_data, '-', color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))
                
                axs.set_title('{} density [1/(m^3)]'.format(kk))
                axs.set_xlabel('$\psi_N$')
                axs.legend(loc= 'best')
                
                
                fig, axs = plt.subplots()
                
                
                for aa in self.data['dircomp']['multi_shift']:
                    
                    b2fstate = self.data['b2fstate'][aa]
                    
                    ne_pro = b2fstate['ne'].transpose()
                    Te_J = b2fstate['te'].transpose()
                    ev = 1.6021766339999999 * pow(10, -19)
                    te_pro = Te_J / ev
                    
                    if kk == 'inner target':
                        psiN = self.data['psi']['psival'][aa][1:37, 1]
                        plot_data = te_pro[1:37, 1]
                        
                    elif kk == 'outer target':
                        
                        psiN = self.data['psi']['psival'][aa][1:37, -2]
                        plot_data = te_pro[1:37, -1]

                    
                        
                    axs.plot(psiN, plot_data, '-', color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))
                
                axs.set_title('{} temperature [eV]'.format(kk))
                axs.set_xlabel('$\psi_N$')
                axs.legend(loc= 'best')
    
    
    
    def rebu_targetNT(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
                     
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
                       
            inne_text = AnchoredText('{}'.format('(a) Inner target electron density [$m^{-3}$]'), 
                                         loc= 'upper right')
            
            inte_text = AnchoredText('{}'.format('(b) Inner target electron temperature [eV]'), 
                                         loc= 'lower right')
            
            infnax_text = AnchoredText('{}'.format('(c) Inner target poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                         loc= 'lower right')
            
            otne_text = AnchoredText('{}'.format('(d) Outer target electron density [$m^{-3}$]'), 
                                         loc= 'upper right')
            
            otte_text = AnchoredText('{}'.format('(e) Outer target electron temperature [eV]'), 
                                         loc= 'lower right')
            
            otfnax_text = AnchoredText('{}'.format('(f) Outer target poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                         loc= 'upper right')
            
            side = ['inner target', 'outer target']
            
            
            fig, axs = plt.subplots(3, 2)
            
            
            for kk in side:
                               
                for aa in self.data['dircomp']['multi_shift']:
                    
                    b2fstate = self.data['b2fstate'][aa]
                    
                    ne_pro = self.data['b2fstate'][aa]['ne'][1:97, 1:37]
                    
                    Te_J = self.data['b2fstate'][aa]['te'][1:97, 1:37]
                    ev = 1.6021766339999999 * pow(10, -19)
                    te_pro = Te_J / ev
                    
                    fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                    vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                    hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                    
                    
                    fnnx = np.divide(fnaxs, vol)
                    fnax = np.multiply(fnnx, hx)
                    
                    data_list = [ne_pro, te_pro, fnax]
                    
                    intitle_list = [inne_text, inte_text, infnax_text]
                    ottitle_list = [otne_text, otte_text, otfnax_text] 
                    
                    if kk == 'inner target':
                        psiN = self.data['psi']['psival'][aa][1:37, 0]
                        
                        for ii, dl in enumerate(data_list):
                            plot_data = dl[0, :]
                            
                            if ii == 0:
                                axs[ii, 0].plot(psiN, plot_data, '-', color = color_dic[aa], 
                                             label = 'A = {}'.format(A_dic[aa]))
                            else:
                                axs[ii, 0].plot(psiN, plot_data, '-', color = color_dic[aa])
                                
                            
                            axs[ii, 0].add_artist(intitle_list[ii])
                            axs[2, 0].set_xlabel('$\psi_N$')
                            axs[0, 0].legend(loc= 'upper left')
                            
                            ii = ii + 1
                            
                    
                    
                    elif kk == 'outer target':
                        
                        psiN = self.data['psi']['psival'][aa][1:37, -1]
                        
                        for ii, dl in enumerate(data_list):
                            plot_data = dl[-1, :]
                            
                            if ii == 0:
                                axs[ii, 1].plot(psiN, plot_data, '-', color = color_dic[aa], 
                                             label = 'A = {}'.format(A_dic[aa]))
                            
                            else:
                                axs[ii, 1].plot(psiN, plot_data, '-', color = color_dic[aa])
                                

                            axs[ii, 1].add_artist(ottitle_list[ii])
                            axs[2, 1].set_xlabel('$\psi_N$')
                            axs[0, 1].legend(loc= 'upper left')
                            
                            ii = ii + 1

                plt.subplots_adjust(hspace=.0)
                plt.subplots_adjust(wspace=0.07)
                # axs.set_title('{} density [1/(m^3)]'.format(kk))
                
    
    def rebu_targetNTi(self):
        
        if self.withshift == True and self.withseries == False:
            
                     
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
                       
            inne_text = AnchoredText('{}'.format('(a) Inner target electron density [$m^{-3}$]'), 
                                         loc= 'upper right')
            
            inte_text = AnchoredText('{}'.format('(b) Inner target electron temperature [eV]'), 
                                         loc= 'lower right')
            
            infnax_text = AnchoredText('{}'.format('(c) Inner target poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                         loc= 'lower right')
            
            otne_text = AnchoredText('{}'.format('(d) Outer target electron density [$m^{-3}$]'), 
                                         loc= 'upper right')
            
            otte_text = AnchoredText('{}'.format('(e) Outer target electron temperature [eV]'), 
                                         loc= 'lower right')
            
            otfnax_text = AnchoredText('{}'.format('(f) Outer target poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                         loc= 'upper right')
            
            side = ['inner target', 'outer target']
            
            
            fig, axs = plt.subplots(3, 2)
            
            
            for kk in side:
                               
                for aa in self.data['dircomp']['multi_shift']:
                    
                    b2fstate = self.data['b2fstate'][aa]
                    
                    ne_pro = self.data['b2fstate'][aa]['na'][1:97, 1:37, 1]
                    
                    Te_J = self.data['b2fstate'][aa]['ti'][1:97, 1:37]
                    ev = 1.6021766339999999 * pow(10, -19)
                    te_pro = Te_J / ev
                    
                    fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                    vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                    hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                    
                    
                    fnnx = np.divide(fnaxs, vol)
                    fnax = np.multiply(fnnx, hx)
                    
                    data_list = [ne_pro, te_pro, fnax]
                    
                    intitle_list = [inne_text, inte_text, infnax_text]
                    ottitle_list = [otne_text, otte_text, otfnax_text] 
                    
                    if kk == 'inner target':
                        psiN = self.data['psi']['psival'][aa][1:37, 0]
                        
                        for ii, dl in enumerate(data_list):
                            plot_data = dl[0, :]
                            
                            if ii == 0:
                                axs[ii, 0].plot(psiN, plot_data, '-', color = color_dic[aa], 
                                             label = 'A = {}'.format(A_dic[aa]))
                            else:
                                axs[ii, 0].plot(psiN, plot_data, '-', color = color_dic[aa])
                                
                            
                            axs[ii, 0].add_artist(intitle_list[ii])
                            axs[2, 0].set_xlabel('$\psi_N$')
                            axs[0, 0].legend(loc= 'upper left')
                            
                            ii = ii + 1
                            
                    
                    
                    elif kk == 'outer target':
                        
                        psiN = self.data['psi']['psival'][aa][1:37, -1]
                        
                        for ii, dl in enumerate(data_list):
                            plot_data = dl[-1, :]
                            
                            if ii == 0:
                                axs[ii, 1].plot(psiN, plot_data, '-', color = color_dic[aa], 
                                             label = 'A = {}'.format(A_dic[aa]))
                            
                            else:
                                axs[ii, 1].plot(psiN, plot_data, '-', color = color_dic[aa])
                                

                            axs[ii, 1].add_artist(ottitle_list[ii])
                            axs[2, 1].set_xlabel('$\psi_N$')
                            axs[0, 1].legend(loc= 'upper left')
                            
                            ii = ii + 1

                plt.subplots_adjust(hspace=.0)
                plt.subplots_adjust(wspace=0.07)
                # axs.set_title('{} density [1/(m^3)]'.format(kk))        
    
                

    
    def check_source(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
                     
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
                       
            div_side_list = ['inner target', 'outer target', 
                             'inner x point boundary', 'outer x point boundary']
            # qu = self.load_iout(filename = dname, simple_quant = quant)
            
            side = ['inner target', 'outer target']
            
            
            for kk in side:
                
                fig, axs = plt.subplots()
                
                
                for aa in self.data['dircomp']['multi_shift']:
                    
                    source = self.data['b2wdat'][aa]['b2npc_sna'][0]
                    # vol = self.data['b2wdat'][aa]['vol']
                    
                    # source_NG = np.divide(source, vol)
                    
                    # S_pro = source_NG.transpose()
                    
                    
                    S_pro = source.transpose()
                    
                    
                    if kk == 'inner target':
                        psiN = self.data['psi']['psival'][aa][1:37, 1]
                        plot_data = S_pro[1:37, 1]
                        
                    elif kk == 'outer target':
                        
                        psiN = self.data['psi']['psival'][aa][1:37, -2]
                        plot_data = S_pro[1:37, -2]

                    
                        
                    axs.plot(psiN, plot_data, '-', color = color_dic[aa], 
                                 label = '{}'.format(A_dic[aa]))
                
                axs.set_title('{} source'.format(kk))
                
                
    
    
    def check_heatsource(self):
        
        if self.withshift == True and self.withseries == False:
            
                     
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            side = ['inner target', 'outer target']
            
            
            for kk in side:
                
                fig, axs = plt.subplots()
                
                
                for aa in self.data['dircomp']['multi_shift']:
                    
                    ehs = self.data['b2wdat'][aa]['b2nph9_she'][1:97, 1:37]
                    evp = self.data['b2wdat'][aa]['b2sihs_shedu'][1:97, 1:37]
                    vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                    ehsa = ehs + evp
                    nehs = np.divide(ehsa, vol)
                    
                    
                    
                    if kk == 'inner target':
                        psiN = self.data['psi']['psival'][aa][1:37, 0]
                        plot_data = nehs[0, :]
                        
                    elif kk == 'outer target':
                        
                        psiN = self.data['psi']['psival'][aa][1:37, -1]
                        plot_data = nehs[-1, :]

                    
                        
                    axs.plot(psiN, plot_data, '-', color = color_dic[aa], 
                                 label = '{}'.format(A_dic[aa]))
                
                axs.set_title('{} heat source'.format(kk))
    
    
    def check_heatflux(self):
        
        if self.withshift == True and self.withseries == False:
            
                     
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            side = ['inner target', 'outer target']
            
            
            for kk in side:
                
                fig, axs = plt.subplots()
                
                
                for aa in self.data['dircomp']['multi_shift']:
                    
                    fhix = self.data['b2wdat'][aa]['b2nph9_fhix'][1:97, 1:37]
                    fhiy = self.data['b2wdat'][aa]['b2nph9_fhiy'][1:97, 1:37]
                    fhex = self.data['b2wdat'][aa]['b2nph9_fhex'][1:97, 1:37]
                    fhey = self.data['b2wdat'][aa]['b2nph9_fhey'][1:97, 1:37]
                    hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
                    hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                    hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                    tor_area = np.multiply(hz, hy)
                    pol_area = np.multiply(hz, hx)
                    fhixn = np.divide(fhix, tor_area)
                    fhiyn = np.divide(fhiy, pol_area)
                    fhexn = np.divide(fhex, tor_area)
                    fheyn = np.divide(fhey, pol_area)
                    
                    
                    if kk == 'inner target':
                        psiN = self.data['psi']['psival'][aa][1:37, 0]
                        plot_data = fhexn[0, :]
                        
                    elif kk == 'outer target':
                        
                        psiN = self.data['psi']['psival'][aa][1:37, -1]
                        plot_data = fhexn[-1, :]

                    
                        
                    axs.plot(psiN, plot_data, '-', color = color_dic[aa], 
                                 label = '{}'.format(A_dic[aa]))
                
                axs.set_title('Electron {} heat flux (W/m^2)'.format(kk))
                # axs.set_yscale('log')
        
        
    
    
    def rebu_tartri(self):
        
        if self.withshift == True and self.withseries == False:
            
                     
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8'}
                       
            inne_text = AnchoredText('{}'.format('(a) $n_e$ [$m^{-3}$]'), 
                                         loc= 'upper right')
            
            inte_text = AnchoredText('{}'.format('(b) $\sqrt{T_e}$ [eV]'), 
                                         loc= 'lower right')
            
            infnax_text = AnchoredText('{}'.format('(c) Poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                         loc= 'lower right')
            
            side = ['inner target']
            
            
            fig, axs = plt.subplots(3, 1)
            
            
            for aa in self.data['dircomp']['multi_shift']:
                
                b2fstate = self.data['b2fstate'][aa]
                
                ne_pro = self.data['b2fstate'][aa]['na'][1:97, 1:37, 1]
                
                Te_J = self.data['b2fstate'][aa]['te'][1:97, 1:37]
                ev = 1.6021766339999999 * pow(10, -19)
                te_pro = Te_J / ev
                
                
                sqrt_te = np.sqrt(te_pro)
                             
                
                fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                
                
                fnnx = np.divide(fnaxs, vol)
                fnax = np.multiply(fnnx, hx)
                
                data_list = [ne_pro, sqrt_te, fnax]
                
                intitle_list = [inne_text, inte_text, infnax_text]
                
                psiN = self.data['psi']['psival'][aa][1:37, 0]
                
                for ii, dl in enumerate(data_list):
                    plot_data = abs(dl[0, :])
                    
                    if ii == 0:
                        axs[ii].plot(psiN, plot_data, '-', color = color_dic[aa], 
                                     label = 'A = {}'.format(A_dic[aa]))
                    else:
                        axs[ii].plot(psiN, plot_data, '-', color = color_dic[aa])
                        
                    
                    axs[ii].add_artist(intitle_list[ii])
                    axs[2].set_xlabel('$\psi_N$')
                    axs[0].legend(loc= 'upper left')
                    axs[0].set_title('Inner target')
                    axs[0].set_yscale('log')
                    axs[1].set_yscale('log')
                    axs[2].set_yscale('log')
                    
                    ii = ii + 1
                
            plt.subplots_adjust(hspace=.0)
            plt.subplots_adjust(wspace=0.07)
                               
                

    
    
    
    
    
    
    
    
    
    

"""

Not working test plan:
    
    
    
gradnax = self.data['b2wdat'][aa]['b2tfnb_dPat_mdf_gradnax'][0][1:97, 1:37]
cdnax = self.data['b2wdat'][aa]['b2trno_cdnax'][0][1:97, 1:37]

print(np.shape(gradnax))

dndx = np.divide(gradnax, cdnax)

na = self.data['b2wdat'][aa]['na'][0]

nadiff = self.diff_quant_x(iout_dat = na)[1:97, 1:37]
hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]

nadx = np.divide(nadiff, hx)

ang_list = self.data['angle']['angle_list'][aa]
# print(np.shape(nadiff))
st = int(pol_list[0])
ed = int(pol_list[-1]) + 1


# axs.plot(ang_list, dndx[st:ed, 1])
axs.plot(ang_list, nadx[st:ed, 1])



def load_b2wdat(self):
    
    if self.withshift == True and self.withseries == False:
        
        b2wdat_dic = {}

        for aa in self.data['dircomp']['multi_shift']:
            
            file_loc = '{}/'.format(self.data['dirdata']['simudir'][aa])
            na_dat = self.data['b2fstate']['org']['na']
            
            
            b2wdat = lBdm.read_b2wdat(b2wdatLoc = file_loc, 
                                      nSpec = np.shape(na_dat)[2])
            b2wdat_dic[aa] = vars(b2wdat)   

        self.data['b2wdat'] = b2wdat_dic

"""






    

