# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 16:22:49 2024

@author: ychuang
"""

from show_flux import show_flow
from matplotlib.offsetbox import AnchoredText
import load_B2_data_method as lBdm
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


class withparticles(show_flow):
    
    
    def __init__(self, DefaultSettings, loadDS):
        show_flow.__init__(self, DefaultSettings, loadDS)
        
        
    
    def allcover_three(self):
        
        
        fig, axs = plt.subplots(2, 1)
            
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        

        
        tpf_text = AnchoredText('{}'.format('(a) Poloidal particle flux $\Gamma_{\Theta}$ [#/s]'), 
                                     loc='upper center')
        
        and_text = AnchoredText('{}'.format('(b) Number of atomic neutrals $n_0$ [#]'), 
                                     loc='upper center')
        
    
        for aa in self.data['dircomp']['multi_shift']:
            
            
            fnax = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
            hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
            hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
            tor_area = np.multiply(hz, hy)
            fnaxs = np.divide(fnax, tor_area)
            
            
            s_term = self.data['b2wdat'][aa]['b2npc_sna'][0][1:97, 1:37]
            vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
            source = np.divide(s_term, vol)
            
            neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
            pn = np.multiply(neuden, vol)
            
            rad_grid = self.data['grid']['RadLoc'][aa]
            vert_grid = self.data['grid']['VertLoc'][aa]
            
            R_loc = rad_grid[18, 1:97][::-1]
            
            Z_loc = vert_grid[18, 1:97][::-1]

            
                      
            arclength, interpolated_points = self.calc_sepx_length(cr = R_loc, 
                                                    cz = Z_loc, plot = False)
            
            fnaxs_list = fnaxs[:, 18]
            source_list = source[:, 18]
            neuden_list = neuden[:, 18]
            fnax_list = abs(fnax[:, 18][::-1])
            fnaxs_list = abs(fnaxs[:, 18][::-1])
            pn_list = pn[:, 18][::-1]
            


            kk = len(fnax_list)
            int_list = list(range(0, kk))
            
            
            axs[0].plot(arclength, fnaxs_list, color = color_dic[aa], 
                         label = 'A = {}'.format(A_dic[aa]))
            axs[1].plot(arclength, fnax_list, color = color_dic[aa])
            # axs[2].plot(arclength, pn_list, color = color_dic[aa])
            
            
            
            

                
        axs[0].legend(loc= 'lower left')
        # axs[0].add_artist(tpf_text)
        # axs[1].add_artist(and_text)
        axs[1].set_xlabel('outer to inner target arc length [m]')
        axs[0].set_yscale('log')
        axs[1].set_yscale('log')
        # axs[2].set_yscale('log')
        axs[0].axvline(x= arclength[0], color='brown', linestyle='--')
        axs[0].axvline(x= arclength[-1], color='gray', linestyle='--')
        axs[0].axvline(x= arclength[23], color='black', linestyle='--')
        axs[0].axvline(x= arclength[71], color='black', linestyle='--')
        axs[1].axvline(x= arclength[0], color='brown', linestyle='--', 
                       label= 'outer target')
        axs[1].axvline(x= arclength[-1], color='gray', linestyle='--', 
                       label= 'inner target')
        axs[1].axvline(x= arclength[23], color='black', linestyle='--', 
                       label= 'x point')
        axs[1].axvline(x= arclength[71], color='black', linestyle='--')
        
        axs[1].legend(loc= 'lower left')
        # axs[2].axvline(x= arclength[0], color='brown', linestyle='--')
        # axs[2].axvline(x= arclength[-1], color='gray', linestyle='--')
        # axs[2].axvline(x= arclength[23], color='black', linestyle='--')
        # axs[2].axvline(x= arclength[71], color='black', linestyle='--')
        
            
            
        plt.subplots_adjust(hspace=.0)
    
    
    
    def allcover_twosum(self):
        
        
        fig, axs = plt.subplots(2, 1)
            
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        

        
        tpf_text = AnchoredText('{}'.format('(a) Poloidal particle flux $\Gamma_{\Theta}$ [#/s]'), 
                                     loc='upper center')
        
        and_text = AnchoredText('{}'.format('(b) Number of atomic neutrals $n_0$ [#]'), 
                                     loc='upper center')
        
    
        for aa in self.data['dircomp']['multi_shift']:
            
            
            fnax = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
            hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
            hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
            hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
            tor_area = np.multiply(hz, hy)
            pol_area = np.multiply(hx, hz)
            fnaxs = np.divide(fnax, tor_area)
            
            
            s_term = self.data['b2wdat'][aa]['b2npc_sna'][0][1:97, 1:37]
            vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
            source = np.divide(s_term, vol)
            
            neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
            pn = np.multiply(neuden, vol)
            sn = np.multiply(source, vol)
            
            rad_grid = self.data['grid']['RadLoc'][aa]
            vert_grid = self.data['grid']['VertLoc'][aa]
            
            R_loc = rad_grid[18, 1:97][::-1]
            
            Z_loc = vert_grid[18, 1:97][::-1]

            
                      
            arclength, interpolated_points = self.calc_sepx_length(cr = R_loc, 
                                                    cz = Z_loc, plot = False)
            
            fnaxs_list = abs(fnaxs[:, 18][::-1])
            source_list = source[:, 18][::-1]
            neuden_list = neuden[:, 18][::-1]
            fnax_list = abs(fnax[:, 18][::-1])
            pn_list = pn[:, 18][::-1]
            sn_list = sn[:, 18][::-1]
            
            
            sfnax_list = abs(fnax[:, 18:][::-1])
            res_fnax_list = np.sum(sfnax_list, axis=1)
            sfnaxs_list = fnaxs[:, 18:][::-1]
            rfnaxs_list = np.sum(sfnaxs_list, axis=1)
            # print("The shape is:")
            # print(np.shape(res_fnax_list))
            
            spn_list = pn[:, 18:][::-1]
            rpn_list = np.sum(spn_list, axis=1)
            snd_list = neuden[:, 18:][::-1]
            rnd_list = np.sum(snd_list, axis=1)
            # print("The shape is:")
            # print(np.shape(rpn_list))
            
            kk = len(fnax_list)
            int_list = list(range(0, kk))
            
            
            axs[0].plot(arclength, sn_list, color = color_dic[aa], 
                         label = 'A = {}'.format(A_dic[aa]))
            axs[1].plot(arclength, neuden_list, color = color_dic[aa])

                
        axs[0].legend(loc= 'lower left')
        # axs[0].add_artist(tpf_text)
        # axs[1].add_artist(and_text)
        axs[1].set_xlabel('outer to inner target arc length')
        axs[0].set_yscale('log')
        axs[1].set_yscale('log')
            
            
        plt.subplots_adjust(hspace=.0)
        plt.title('twosum')
    
    
    
    
    
    
    # def term_one(self, pol_list):
        
    #     if self.withshift == True and self.withseries == False:
            
    #         color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
    #                      'dot7': 'blue', 'one': 'purple'}
            
    #         A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
    #                   'dot7': '2.8', 'one': '3.4'}
            
    #         fig, axs = plt.subplots()
                      
    #         for aa in self.data['dircomp']['multi_shift']:
                
                    
                
    #             fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
    #             vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
    #             hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
    #             hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                
                
    #             g_coe = np.divide(vol, hx)
                
    #             fnax = np.divide(fnaxs, g_coe)
    #             cup_fnax_org = np.multiply(fnax, hy)
                
                
    #             cpfnaxdiff_org = self.diff_quant_x(iout_dat = cup_fnax_org)
    #             acp = np.divide(cpfnaxdiff_org, hx)
    #             agcoe = np.multiply(hx_org, hy_org)
    #             t1 = np.divide(acp, agcoe)
                
    #             del_t1 = t1 - t1_org
                
                
    #             ang_list = self.data['angle']['angle_list'][aa]
    #             # print(np.shape(nadiff))
    #             st = int(pol_list[0])
    #             ed = int(pol_list[-1]) + 1
                
                
    #             # axs.plot(ang_list, dndx[st:ed, 1])
    #             axs.plot(ang_list, del_t1[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                    
                
                
            
    #         axs.legend(loc= 'best')
    #         axs.set_title('term 1')
    #         axs.set_xlabel('poloidal angle')
        
        
        
        