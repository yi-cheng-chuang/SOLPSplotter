# -*- coding: utf-8 -*-
"""
Created on Sun Oct 20 01:49:17 2024

@author: ychuang
"""


from R_diff_calc import Diff_R_calc
from matplotlib.offsetbox import AnchoredText
import load_B2_data_method as lBdm
import numpy as np
import matplotlib.pyplot as plt


class show_flow(Diff_R_calc):
    
    
    def __init__(self, DefaultSettings, loadDS):
        Diff_R_calc.__init__(self, DefaultSettings, loadDS)
    
    
    
    def totpolflux(self):
        
            
        fig, axs = plt.subplots(2, 1)
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        HFS_text = AnchoredText('{}'.format('(a) Total poloidal flux $\Gamma_{\Theta}$ [$s^{-1}$] at HFS'), 
                                     loc='upper center')
        
        LFS_text = AnchoredText('{}'.format('(b) Total poloidal flux $\Gamma_{\Theta}$ [$s^{-1}$] at LFS'), 
                                     loc='upper center')
        
        for side in ['HFS', 'LFS']:
        
        
            for aa in self.data['dircomp']['multi_shift']:
                
                fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                
                # hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
                # fnaxs = np.divide(fnax, hz)
                       
                if side == 'HFS':
                    
                    index_list = np.linspace(0, 25, 26)
                    
                    fx_list = []
                    
                    for ii in index_list:
                        
                        
                        fx_list.append(abs(sum(fnaxs[int(ii), 18:])))
                    
                    axs[0].plot(index_list, fx_list, color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))
                    
                    
                elif side == 'LFS':
                    
                    index_list = np.linspace(60, 95, 36)
                    
                    fx_list = []
                    
                    for ii in index_list:
                        
                        fx_list.append(sum(fnaxs[int(ii), 18:]))
                    
                    axs[1].plot(index_list, fx_list, color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))

                

        axs[0].axvline(x= 0,ls = '--', color = 'black', label = 'inner target')
        axs[0].axvline(x= 25,ls = '--', color = 'brown', label = 'poloidal angle 250')
        axs[1].axvline(x= 95, ls= '--', color = 'black', label = 'outer target')
        axs[1].axvline(x= 60,ls = '--', color = 'brown', label = 'outer midplane')
        axs[0].axvline(x= 24, ls= '--', color = 'gray', label = 'left cut')
        axs[1].axvline(x= 72,ls = '--', color = 'gray', label = 'right cut')
        axs[0].set_yscale('log')
        axs[1].set_yscale('log')
        axs[1].legend(loc= 'upper right')
        axs[0].legend(loc = 'upper right')
        axs[0].add_artist(HFS_text)
        axs[1].add_artist(LFS_text)
        axs[0].set_xlabel('poloidal index')
        axs[1].set_xlabel('poloidal index')
    
    
    def totpolfluxNR(self):
        
            
        fig, axs = plt.subplots(2, 1)
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        HFS_text = AnchoredText('{}'.format('(a) Poloidal flux $\Gamma_{\Theta}$ [$m^{-1} s^{-1}$] at HFS'), 
                                     loc='upper center')
        
        LFS_text = AnchoredText('{}'.format('(b) Poloidal flux $\Gamma_{\Theta}$ [$m^{-1} s^{-1}$] at LFS'), 
                                     loc='upper center')
        
        for side in ['HFS', 'LFS']:
        
        
            for aa in self.data['dircomp']['multi_shift']:
                
                fnax = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                
                hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
                fnaxs = np.divide(fnax, hz)
                       
                if side == 'HFS':
                    
                    index_list = np.linspace(0, 25, 26)
                    
                    fx_list = []
                    
                    for ii in index_list:
                        
                        
                        fx_list.append(abs(sum(fnaxs[int(ii), 18:])))
                    
                    axs[0].plot(index_list, fx_list, color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))
                    
                    
                elif side == 'LFS':
                    
                    index_list = np.linspace(60, 95, 36)
                    
                    fx_list = []
                    
                    for ii in index_list:
                        
                        fx_list.append(sum(fnaxs[int(ii), 18:]))
                    
                    axs[1].plot(index_list, fx_list, color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))

                

        axs[0].axvline(x= 0,ls = '--', color = 'black', label = 'inner target')
        axs[0].axvline(x= 25,ls = '--', color = 'brown', label = 'poloidal angle 250')
        axs[1].axvline(x= 95, ls= '--', color = 'black', label = 'outer target')
        axs[1].axvline(x= 60,ls = '--', color = 'brown', label = 'outer midplane')
        axs[0].axvline(x= 24, ls= '--', color = 'gray', label = 'left cut')
        axs[1].axvline(x= 72,ls = '--', color = 'gray', label = 'right cut')
        axs[0].set_yscale('log')
        axs[1].set_yscale('log')
        axs[1].legend(loc= 'upper right')
        axs[0].legend(loc = 'upper right')
        axs[0].add_artist(HFS_text)
        axs[1].add_artist(LFS_text)
        axs[0].set_xlabel('poloidal index')
        axs[1].set_xlabel('poloidal index')
    
    
        
    def totfluxes(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            
            fig, axs = plt.subplots(3, 1)
            
            for aa in self.data['dircomp']['multi_shift']:
                
                fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                source = self.data['b2wdat'][aa]['b2npc_sna'][0][1:97, 1:37]
                
                neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                nd = np.multiply(neuden, vol)
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                fx_list = []
                s_list = []
                nd_list = []
                
                for ii in pol_list:
                    
                    fx_list.append(sum(fnaxs[int(ii), 18:]))
                    s_list.append(sum(source[int(ii), 18:]))
                    nd_list.append(sum(nd[int(ii), 18:]))
                
                
                    
                axs[0].plot(ang_list, fx_list, color = color_dic[aa], 
                            label = '{}'.format(A_dic[aa]))
                axs[1].plot(ang_list, s_list, color = color_dic[aa])
                axs[2].plot(ang_list, nd_list, color = color_dic[aa])
                    
            
            # axs[0].set_xlabel('poloidal angle')
            # axs[1].set_xlabel('poloidal angle')
            axs[2].set_xlabel('poloidal angle')
            axs[1].set_yscale('log')
            axs[2].set_yscale('log')
            plt.subplots_adjust(hspace=.0)
    
    
    def triNR_HFS(self, side):
        
        if side == 'HFS':
            fig, axs = plt.subplots(4, 1)
        
        elif side == 'LFS':
            
            fig, axs = plt.subplots(3, 1)
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        cspf_text = AnchoredText('{}'.format('(a) Cross section poloidal flux $\Gamma_{\Theta}$ [$m^{-1} s^{-1}$]'), 
                                     loc='upper center')
        
        as_text = AnchoredText('{}'.format('(b) Area source [$m^{-1} s^{-1}$]'), 
                                     loc='upper center')
        
        and_text = AnchoredText('{}'.format('(c) Area atomic neutral density $n_0$ [$m^{-1}$]'), 
                                     loc='upper center')
        
        tpf_text = AnchoredText('{}'.format('(d) Poloidal flux $\Gamma_{\Theta}$ [$s^{-1}$]'), 
                                     loc='upper center')
        

        for aa in self.data['dircomp']['multi_shift']:
            
            fnax = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
            hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
            fnaxs = np.divide(fnax, hz)
            
            
            s_term = self.data['b2wdat'][aa]['b2npc_sna'][0][1:97, 1:37]
            source = np.divide(s_term, hz)
            
            neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
            vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
            hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
            hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
            area = np.multiply(hx, hy)
            nd = np.multiply(neuden, area)
            
            if side == 'HFS':
                
                index_list = np.linspace(0, 25, 26)
            
            elif side == 'LFS':
                
                index_list = np.linspace(60, 95, 36)
                
                
            fx_list = []
            s_list = []
            nd_list = []
            fxt_list = []
            
            for ii in index_list:
                
                
                fx_list.append(abs(sum(fnaxs[int(ii), 18:])))
                fxt_list.append(abs(sum(fnax[int(ii), 18:])))
                s_list.append(sum(source[int(ii), 18:]))
                nd_list.append(sum(nd[int(ii), 18:]))
            
            axs[0].plot(index_list, fx_list, color = color_dic[aa], 
                         label = 'A = {}'.format(A_dic[aa]))
            axs[1].plot(index_list, s_list, color = color_dic[aa])
            axs[2].plot(index_list, nd_list, color = color_dic[aa])
                
            if side == 'HFS':
                
                axs[3].plot(index_list, fxt_list, color = color_dic[aa])
                
                

                
        if side == 'HFS':
            
            axs[0].axvline(x= 0,ls = '--', color = 'black')
            axs[0].axvline(x= 24,ls = '--', color = 'gray')
            axs[0].axvline(x= 25,ls = '--', color = 'brown')
            axs[1].axvline(x= 0, ls= '--', color = 'black', label = 'inner target')
            axs[1].axvline(x= 24,ls = '--', color = 'gray', label = 'x point')
            axs[1].axvline(x= 25,ls = '--', color = 'brown', label = 'poloidal angle 250')
            axs[2].axvline(x= 0, ls= '--', color = 'black')
            axs[2].axvline(x= 24,ls = '--', color = 'gray')
            axs[2].axvline(x= 25,ls = '--', color = 'brown')
            axs[3].axvline(x= 0, ls= '--', color = 'black')
            axs[3].axvline(x= 24,ls = '--', color = 'gray')
            axs[3].axvline(x= 25,ls = '--', color = 'brown')
        
        elif side == 'LFS':
            
            axs[0].axvline(x= 95,ls = '--', color = 'black')
            axs[0].axvline(x= 72,ls = '--', color = 'gray')
            axs[0].axvline(x= 60,ls = '--', color = 'brown')
            axs[1].axvline(x= 95, ls= '--', color = 'black', label = 'outer target')
            axs[1].axvline(x= 72,ls = '--', color = 'gray', label = 'x point')
            axs[1].axvline(x= 60,ls = '--', color = 'brown', label = 'outer midplane')
            axs[2].axvline(x= 95, ls= '--', color = 'black')
            axs[2].axvline(x= 72,ls = '--', color = 'gray')
            axs[2].axvline(x= 60,ls = '--', color = 'brown')
            
        
        axs[1].set_yscale('log')
        axs[2].set_yscale('log')
        axs[0].legend(loc= 'upper right')
        axs[1].legend(loc = 'upper right')
        axs[0].add_artist(cspf_text)
        axs[1].add_artist(as_text)
        axs[2].add_artist(and_text)
        
        
        if side == 'HFS':
            axs[3].set_xlabel('poloidal index')
            axs[3].add_artist(tpf_text)
        
        elif side == 'LFS':
            
            axs[2].set_xlabel('poloidal index')
            
            
        plt.subplots_adjust(hspace=.0)
        plt.suptitle('{}'.format(side))
    
    
    
    def triNR_eq(self, side):
        

            
        fig, axs = plt.subplots(3, 1)
            
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        cspf_text = AnchoredText('{}'.format('(a) Cross section poloidal flux $\Gamma_{\Theta}$ [$m^{-1} s^{-1}$]'), 
                                     loc='upper center')
        
        as_text = AnchoredText('{}'.format('(b) Area source [$m^{-1} s^{-1}$]'), 
                                     loc='upper center')
        
        and_text = AnchoredText('{}'.format('(c) Area atomic neutral density $n_0$ [$m^{-1}$]'), 
                                     loc='upper center')
        
        tpf_text = AnchoredText('{}'.format('(d) Poloidal flux $\Gamma_{\Theta}$ [$s^{-1}$]'), 
                                     loc='upper center')
        

        for aa in self.data['dircomp']['multi_shift']:
            
            fnax = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
            hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
            fnaxs = np.divide(fnax, hz)
            
            
            s_term = self.data['b2wdat'][aa]['b2npc_sna'][0][1:97, 1:37]
            source = np.divide(s_term, hz)
            
            neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
            vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
            hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
            hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
            area = np.multiply(hx, hy)
            nd = np.multiply(neuden, area)
            
            if side == 'HFS':
                
                index_list = np.linspace(0, 25, 26)
            
            elif side == 'LFS':
                
                index_list = np.linspace(60, 95, 36)
                
                
            fx_list = []
            s_list = []
            nd_list = []
            fxt_list = []
            
            for ii in index_list:
                
                
                fx_list.append(abs(sum(fnaxs[int(ii), 18:])))
                fxt_list.append(abs(sum(fnax[int(ii), 18:])))
                s_list.append(sum(source[int(ii), 18:]))
                nd_list.append(sum(nd[int(ii), 18:]))
            
            axs[0].plot(index_list, fx_list, color = color_dic[aa], 
                         label = 'A = {}'.format(A_dic[aa]))
            axs[1].plot(index_list, s_list, color = color_dic[aa])
            axs[2].plot(index_list, nd_list, color = color_dic[aa])
                
            
        if side == 'HFS':
            
            axs[0].axvline(x= 0,ls = '--', color = 'black')
            axs[0].axvline(x= 24,ls = '--', color = 'gray')
            axs[0].axvline(x= 25,ls = '--', color = 'brown')
            axs[1].axvline(x= 0, ls= '--', color = 'black', label = 'inner target')
            axs[1].axvline(x= 24,ls = '--', color = 'gray', label = 'x point')
            axs[1].axvline(x= 25,ls = '--', color = 'brown', label = 'poloidal angle 250')
            axs[2].axvline(x= 0, ls= '--', color = 'black')
            axs[2].axvline(x= 24,ls = '--', color = 'gray')
            axs[2].axvline(x= 25,ls = '--', color = 'brown')
        
        elif side == 'LFS':
            
            axs[0].axvline(x= 95,ls = '--', color = 'black')
            axs[0].axvline(x= 72,ls = '--', color = 'gray')
            axs[0].axvline(x= 60,ls = '--', color = 'brown')
            axs[1].axvline(x= 95, ls= '--', color = 'black', label = 'outer target')
            axs[1].axvline(x= 72,ls = '--', color = 'gray', label = 'x point')
            axs[1].axvline(x= 60,ls = '--', color = 'brown', label = 'outer midplane')
            axs[2].axvline(x= 95, ls= '--', color = 'black')
            axs[2].axvline(x= 72,ls = '--', color = 'gray')
            axs[2].axvline(x= 60,ls = '--', color = 'brown')
            
        
        axs[1].set_yscale('log')
        axs[2].set_yscale('log')
        axs[0].legend(loc= 'upper right')
        axs[1].legend(loc = 'upper right')
        axs[0].add_artist(cspf_text)
        axs[1].add_artist(as_text)
        axs[2].add_artist(and_text)
        axs[2].set_xlabel('poloidal index')
            
            
        plt.subplots_adjust(hspace=.0)
        plt.suptitle('{}'.format(side))
    
    
    
    def triNG(self, side):
        

            
        fig, axs = plt.subplots(3, 1)
            
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        cspf_text = AnchoredText('{}'.format('(a) Poloidal flux $\Gamma_{\Theta}$ [$m^{-3} s^{-1}$]'), 
                                     loc='upper center')
        
        as_text = AnchoredText('{}'.format('(b) Source [$m^{-3} s^{-1}$]'), 
                                     loc='upper center')
        
        and_text = AnchoredText('{}'.format('(c) Atomic neutral density $n_0$ [$m^{-3}$]'), 
                                     loc='upper center')
        
        tpf_text = AnchoredText('{}'.format('(d) Poloidal flux $\Gamma_{\Theta}$ [$s^{-1}$]'), 
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
            
            # hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
            # hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
            # area = np.multiply(hx, hy)
            # nd = np.multiply(neuden, area)
            
            if side == 'HFS':
                
                index_list = np.linspace(0, 25, 26)
            
            elif side == 'LFS':
                
                index_list = np.linspace(60, 95, 36)
                
                
            fx_list = []
            s_list = []
            nd_list = []
            
            for ii in index_list:
                
                
                fx_list.append(abs(fnaxs[int(ii), 18]))
                s_list.append(source[int(ii), 18])
                nd_list.append(neuden[int(ii), 18])
            
            axs[0].plot(index_list, fx_list, color = color_dic[aa], 
                         label = 'A = {}'.format(A_dic[aa]))
            axs[1].plot(index_list, s_list, color = color_dic[aa])
            axs[2].plot(index_list, nd_list, color = color_dic[aa])
                
            
        if side == 'HFS':
            
            axs[0].axvline(x= 0,ls = '--', color = 'black')
            axs[0].axvline(x= 24,ls = '--', color = 'gray')
            axs[0].axvline(x= 25,ls = '--', color = 'brown')
            axs[1].axvline(x= 0, ls= '--', color = 'black', label = 'inner target')
            axs[1].axvline(x= 24,ls = '--', color = 'gray', label = 'x point')
            axs[1].axvline(x= 25,ls = '--', color = 'brown', label = 'poloidal angle 250')
            axs[2].axvline(x= 0, ls= '--', color = 'black')
            axs[2].axvline(x= 24,ls = '--', color = 'gray')
            axs[2].axvline(x= 25,ls = '--', color = 'brown')
        
        elif side == 'LFS':
            
            axs[0].axvline(x= 95,ls = '--', color = 'black')
            axs[0].axvline(x= 72,ls = '--', color = 'gray')
            axs[0].axvline(x= 60,ls = '--', color = 'brown')
            axs[1].axvline(x= 95, ls= '--', color = 'black', label = 'outer target')
            axs[1].axvline(x= 72,ls = '--', color = 'gray', label = 'x point')
            axs[1].axvline(x= 60,ls = '--', color = 'brown', label = 'outer midplane')
            axs[2].axvline(x= 95, ls= '--', color = 'black')
            axs[2].axvline(x= 72,ls = '--', color = 'gray')
            axs[2].axvline(x= 60,ls = '--', color = 'brown')
            
        
        axs[1].set_yscale('log')
        axs[2].set_yscale('log')
        axs[0].legend(loc= 'upper right')
        axs[1].legend(loc = 'upper right')
        axs[0].add_artist(cspf_text)
        axs[1].add_artist(as_text)
        axs[2].add_artist(and_text)
        axs[2].set_xlabel('poloidal index')
            
            
        plt.subplots_adjust(hspace=.0)
        plt.suptitle('{}'.format(side))
    
     
    
    def allcover_fhix(self):
        
        
        fig, axs = plt.subplots(4, 1)
            
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        fhix_text = AnchoredText('{}'.format('(a) Ion poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
                                     loc='upper center')
        
        fhiy_text = AnchoredText('{}'.format('(b) Ion radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
                                     loc='upper center')
        
        fhex_text = AnchoredText('{}'.format('(c) Electron poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
                                     loc='upper center')
        
        fhey_text = AnchoredText('{}'.format('(d) Electron radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
                                     loc='upper center')
        
    
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
            
            
            kk = len(fhixn)
            int_list = list(range(0, kk))
                
                
            fx_list = fhixn[:, 18]
            fy_list = fhiyn[:, -1]
            fex_list = fhixn[:, 18]
            fey_list = fhiyn[:, -1]
            
            
            axs[0].plot(int_list, fx_list, color = color_dic[aa], 
                         label = 'A = {}'.format(A_dic[aa]))
            axs[1].plot(int_list, fy_list, color = color_dic[aa])
            axs[2].plot(int_list, fex_list, color = color_dic[aa])
            axs[3].plot(int_list, fey_list, color = color_dic[aa])
                
        axs[0].legend(loc= 'upper right')
        axs[0].add_artist(fhix_text)
        axs[1].add_artist(fhiy_text)
        axs[2].add_artist(fhex_text)
        axs[3].add_artist(fhey_text)
        axs[1].set_xlabel('poloidal index')
            
            
        plt.subplots_adjust(hspace=.0)
    
    
    
    def allcover_tri(self):
        
        
        fig, axs = plt.subplots(4, 1)
            
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        cspf_text = AnchoredText('{}'.format('(a) Poloidal flux $\Gamma_{\Theta}$ [$m^{-3} s^{-1}$]'), 
                                     loc='upper center')
        
        as_text = AnchoredText('{}'.format('(b) Source [$m^{-3} s^{-1}$]'), 
                                     loc='upper center')
        
        and_text = AnchoredText('{}'.format('(c) Atomic neutral density $n_0$ [$m^{-3}$]'), 
                                     loc='upper center')
        
        tpf_text = AnchoredText('{}'.format('(d) Poloidal particle flux $\Gamma_{\Theta}$ [$s^{-1}$]'), 
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
            
            
                
            fnaxs_list = fnaxs[:, 18]
            source_list = source[:, 18]
            neuden_list = neuden[:, 18]
            fnax_list = fnax[:, 18]
            
            kk = len(fnax_list)
            int_list = list(range(0, kk))
            
            
            axs[0].plot(int_list, fnaxs_list, color = color_dic[aa], 
                         label = 'A = {}'.format(A_dic[aa]))
            axs[1].plot(int_list, source_list, color = color_dic[aa])
            axs[2].plot(int_list, neuden_list, color = color_dic[aa])
            axs[3].plot(int_list, fnax_list, color = color_dic[aa])
                
        axs[0].legend(loc= 'upper right')
        axs[0].add_artist(cspf_text)
        axs[1].add_artist(as_text)
        axs[2].add_artist(and_text)
        axs[3].add_artist(tpf_text)
        axs[1].set_xlabel('poloidal index')
            
            
        plt.subplots_adjust(hspace=.0)
    
    
    
    def allcover_neuden(self):
        
        
        fig, axs = plt.subplots()
            
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        

        
        tpf_text = AnchoredText('{}'.format('Poloidal particle flux $\Gamma_{\Theta}$ [$s^{-1}$]'), 
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
            
            
                
            fnaxs_list = fnaxs[:, 18]
            source_list = source[:, 18]
            neuden_list = neuden[:, 18]
            fnax_list = fnax[:, 18]
            pn_list = pn[:, 18]
            
            kk = len(fnax_list)
            int_list = list(range(0, kk))
            
            

            axs.plot(int_list, pn_list, color = color_dic[aa])
                
        axs.legend(loc= 'upper right')
        axs.add_artist(tpf_text)
        axs.set_xlabel('poloidal index')
        axs.set_yscale('log')
            
            
        plt.subplots_adjust(hspace=.0)
    
    

    
    
    
    
    
    def spac_fhix(self):
        
        
        fig, axs = plt.subplots()
            
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        # fhix_text = AnchoredText('{}'.format('(a) Ion poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        fhiy_text = AnchoredText('{}'.format('(a) Ion radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
                                     loc='upper center')
        
        # fhex_text = AnchoredText('{}'.format('(c) Electron poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        # fhey_text = AnchoredText('{}'.format('(d) Electron radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
    
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
            
            
            kk = len(fhixn)
            int_list = list(range(0, kk))
                
                
            # fx_list = fhixn[:, 18]
            fy_list = fhiyn[:, 0]
            # fex_list = fhixn[:, 18]
            # fey_list = fhiyn[:, -1]
            
            
            # axs[0].plot(int_list, fx_list, color = color_dic[aa], 
            #              label = 'A = {}'.format(A_dic[aa]))
            axs.plot(int_list, fy_list, color = color_dic[aa], 
                        label = 'A = {}'.format(A_dic[aa]))
            # axs[2].plot(int_list, fex_list, color = color_dic[aa])
            # axs[3].plot(int_list, fey_list, color = color_dic[aa])
                
        axs.legend(loc= 'upper right')
        # axs[0].add_artist(fhix_text)
        axs.add_artist(fhiy_text)
        # axs[2].add_artist(fhex_text)
        # axs[3].add_artist(fhey_text)
        axs.set_xlabel('poloidal index')
            
            
        plt.subplots_adjust(hspace=.0)
    
    
    
    
    def hflux_bar(self):
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        # fhix_text = AnchoredText('{}'.format('(a) Ion poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        fhiy_text = AnchoredText('{}'.format('(a) Ion radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
                                     loc='upper center')
        
        # fhex_text = AnchoredText('{}'.format('(c) Electron poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        # fhey_text = AnchoredText('{}'.format('(d) Electron radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        

        def return_value():
            vlist = []
            gl_list = []
            cl_list = []
            for ii, aa in enumerate(self.data['dircomp']['multi_shift']):
                
                
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
                
                
                kk = len(fhixn)
                int_list = list(range(0, kk))
                    
                    
                # fx_list = fhixn[:, 18]
                if aa == 'org':
                    fyi = sum(fhiy[24:72, 3])
                elif aa == 'dot3':
                    fyi = sum(fhiy[24:72, 0])
                else:
                    fyi = sum(fhiy[24:72, 2])
                fyb = sum(fhiy[:, -1])
                fxti = sum(abs(fhix[0, :]))
                fxto = sum(fhix[-1, :])
                fypfr1 = sum(fhiy[:24, 0])
                fypfr2 = sum(fhiy[72:, 0])
                fypfr = fypfr1 + fypfr2
                fsum = fyb + fxti + fxto


                vlist.append([fyb/fyi, fxti/fyi, fxto/fyi])
                gl_list.append('A = {}'.format(A_dic[aa]))
                cl_list.append('{}'.format(color_dic[aa]))
            
            return vlist, gl_list, cl_list
            
            
            
            
        # Data
        
        vlist, gl_list, cl_list = return_value()
        
        categories = ['SOL boundary', 'inner target', 'outer target']  # Categories
        values = vlist
        group_labels = gl_list
        colors = cl_list

        # Number of categories and groups
        n_categories = len(categories)
        n_groups = len(values)

        # X positions for the categories
        x = np.arange(n_categories)

        # Bar width
        bar_width = 0.2

        # Create the plot
        plt.figure()

        # Loop through each group to plot
        for i, group_values in enumerate(values):
            plt.bar(
                x + i * bar_width,  # Shift each group by bar_width
                group_values,       # Heights of the bars
                width=bar_width,    # Width of the bars
                label=group_labels[i],  # Label for the group
                color=colors[i]     # Color for the group
            )

        # Add labels, title, and legend
        plt.xlabel('Categories')
        # plt.ylabel('Values')
        plt.title('Ion Power (W) fraction Bar Chart')
        plt.xticks(x + bar_width * (n_groups - 1) / 2, categories)  # Center the category labels
        plt.legend()

        # Show the plot
        plt.tight_layout()
        plt.show()
        
    
        
        
    def hfluxe_bar(self):
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        # fhix_text = AnchoredText('{}'.format('(a) Ion poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        fhiy_text = AnchoredText('{}'.format('(a) Ion radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
                                     loc='upper center')
        
        # fhex_text = AnchoredText('{}'.format('(c) Electron poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        # fhey_text = AnchoredText('{}'.format('(d) Electron radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        def return_value():
            vlist = []
            gl_list = []
            cl_list = []
            for ii, aa in enumerate(self.data['dircomp']['multi_shift']):
                
                
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
                
                
                kk = len(fhixn)
                int_list = list(range(0, kk))
                    
                    
                # fx_list = fhixn[:, 18]
                if aa == 'org':
                    fyi = sum(fhey[24:72, 3])
                elif aa == 'dot3':
                    fyi = sum(fhey[24:72, 0])
                else:
                    fyi = sum(fhey[24:72, 2])
                fyb = sum(fhey[:, -1])
                fxti = sum(abs(fhex[0, :]))
                fxto = sum(fhex[-1, :])
                fypfr1 = sum(fhey[:24, 0])
                fypfr2 = sum(fhey[72:, 0])
                fypfr = fypfr1 + fypfr2
                
                vlist.append([fyb/fyi, fxti/fyi, fxto/fyi])
                gl_list.append('A = {}'.format(A_dic[aa]))
                cl_list.append('{}'.format(color_dic[aa]))
            
            return vlist, gl_list, cl_list
                
                
        # Data
        
        vlist, gl_list, cl_list = return_value()
        
        categories = ['SOL boundary', 'inner target', 'outer target']  # Categories
        values = vlist
        group_labels = gl_list
        colors = cl_list

        # Number of categories and groups
        n_categories = len(categories)
        n_groups = len(values)

        # X positions for the categories
        x = np.arange(n_categories)

        # Bar width
        bar_width = 0.2

        # Create the plot
        plt.figure()

        # Loop through each group to plot
        for i, group_values in enumerate(values):
            plt.bar(
                x + i * bar_width,  # Shift each group by bar_width
                group_values,       # Heights of the bars
                width=bar_width,    # Width of the bars
                label=group_labels[i],  # Label for the group
                color=colors[i]     # Color for the group
            )

        # Add labels, title, and legend
        plt.xlabel('Categories')
        # plt.ylabel('Values')
        plt.title('Electron Power (W) fraction Bar Chart')
        plt.xticks(x + bar_width * (n_groups - 1) / 2, categories)  # Center the category labels
        plt.legend()

        # Show the plot
        plt.tight_layout()
        plt.show()
        
    
    
    
    
    def hflux_tot_bar(self):
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        # fhix_text = AnchoredText('{}'.format('(a) Ion poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        fhiy_text = AnchoredText('{}'.format('(a) Ion radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
                                     loc='upper center')
        
        # fhex_text = AnchoredText('{}'.format('(c) Electron poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        # fhey_text = AnchoredText('{}'.format('(d) Electron radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        def return_value():
            vlist = []
            gl_list = []
            cl_list = []
            for ii, aa in enumerate(self.data['dircomp']['multi_shift']):
                
                
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
                
                
                kk = len(fhixn)
                int_list = list(range(0, kk))
                    
                    
                # fx_list = fhixn[:, 18]
                if aa == 'org':
                    fyi = sum(fhey[24:72, 3])/sum(pol_area[24:72, 3])
                elif aa == 'dot3':
                    fyi = sum(fhey[24:72, 0])/sum(pol_area[24:72, 0])
                else:
                    fyi = sum(fhey[24:72, 2])/sum(pol_area[24:72, 2])
                fyb = sum(fhey[:, -1])/sum(pol_area[:, -1])
                fxti = sum(abs(fhex[0, 18:]))/sum(tor_area[0, 18:])
                fxto = sum(fhex[-1, 18:])/sum(tor_area[-1, 18:])
                fypfr1 = sum(fhiy[:24, 0])
                fypfr2 = sum(fhiy[72:, 0])
                fypfr = fypfr1 + fypfr2
                
                
                vlist.append([fyi, fyb, fxti, fxto])
                gl_list.append('A = {}'.format(A_dic[aa]))
                cl_list.append('{}'.format(color_dic[aa]))
            
            return vlist, gl_list, cl_list

        # Data
        
        vlist, gl_list, cl_list = return_value()
        
        categories = ['core input', 'SOL boundary', 'inner target', 'outer target']  # Categories
        values = vlist
        group_labels = gl_list
        colors = cl_list

        # Number of categories and groups
        n_categories = len(categories)
        n_groups = len(values)

        # X positions for the categories
        x = np.arange(n_categories)

        # Bar width
        bar_width = 0.2

        # Create the plot
        plt.figure()

        # Loop through each group to plot
        for i, group_values in enumerate(values):
            plt.bar(
                x + i * bar_width,  # Shift each group by bar_width
                group_values,       # Heights of the bars
                width=bar_width,    # Width of the bars
                label=group_labels[i],  # Label for the group
                color=colors[i]     # Color for the group
            )

        # Add labels, title, and legend
        # plt.xlabel('Categories')
        # plt.ylabel('Values')
        plt.yscale('log')
        plt.title('Average electron heat flux (W/$m^2$) value Bar Chart')
        plt.xticks(x + bar_width * (n_groups - 1) / 2, categories)  # Center the category labels
        plt.legend()

        # Show the plot
        plt.tight_layout()
        plt.show()
    
    
    
    def hflux_value_bar(self):
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        # fhix_text = AnchoredText('{}'.format('(a) Ion poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        fhiy_text = AnchoredText('{}'.format('(a) Ion radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
                                     loc='upper center')
        
        # fhex_text = AnchoredText('{}'.format('(c) Electron poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        # fhey_text = AnchoredText('{}'.format('(d) Electron radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        def return_value():
            vlist = []
            gl_list = []
            cl_list = []
            for ii, aa in enumerate(self.data['dircomp']['multi_shift']):
                
                
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
                
                
                kk = len(fhixn)
                int_list = list(range(0, kk))
                    
                    
                # fx_list = fhixn[:, 18]
                if aa == 'org':
                    fyi = sum(fhey[24:72, 3])
                elif aa == 'dot3':
                    fyi = sum(fhey[24:72, 0])
                else:
                    fyi = sum(fhey[24:72, 2])
                fyb = sum(fhey[:, -1])
                fxti = sum(abs(fhex[0, :]))
                fxto = sum(fhex[-1, :])
                fypfr1 = sum(fhiy[:24, 0])
                fypfr2 = sum(fhiy[72:, 0])
                fypfr = fypfr1 + fypfr2
                
                
                vlist.append([fyi, fyb, fxti, fxto])
                gl_list.append('A = {}'.format(A_dic[aa]))
                cl_list.append('{}'.format(color_dic[aa]))
            
            return vlist, gl_list, cl_list

        # Data
        
        vlist, gl_list, cl_list = return_value()
        
        categories = ['core input', 'north boundary', 'inner target', 'outer target']  # Categories
        values = vlist
        group_labels = gl_list
        colors = cl_list

        # Number of categories and groups
        n_categories = len(categories)
        n_groups = len(values)

        # X positions for the categories
        x = np.arange(n_categories)

        # Bar width
        bar_width = 0.2

        # Create the plot
        plt.figure()

        # Loop through each group to plot
        for i, group_values in enumerate(values):
            plt.bar(
                x + i * bar_width,  # Shift each group by bar_width
                group_values,       # Heights of the bars
                width=bar_width,    # Width of the bars
                label=group_labels[i],  # Label for the group
                color=colors[i]     # Color for the group
            )

        # Add labels, title, and legend
        plt.xlabel('Categories')
        plt.ylabel('Values')
        plt.title('Electron Power (W) value Bar Chart')
        plt.xticks(x + bar_width * (n_groups - 1) / 2, categories)  # Center the category labels
        plt.legend()

        # Show the plot
        plt.tight_layout()
        plt.show()
    
    
    
    def pflux_tot_bar(self):
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        # fhix_text = AnchoredText('{}'.format('(a) Ion poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        fhiy_text = AnchoredText('{}'.format('(a) Ion radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
                                     loc='upper center')
        
        # fhex_text = AnchoredText('{}'.format('(c) Electron poloidal heat flux $q_{i\Theta}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        # fhey_text = AnchoredText('{}'.format('(d) Electron radial heat flux $q_{ir}$ [$m^{-2} s^{-1}$]'), 
        #                              loc='upper center')
        
        def return_value():
            vlist = []
            gl_list = []
            cl_list = []
            for ii, aa in enumerate(self.data['dircomp']['multi_shift']):
                
                fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                fnays = self.data['b2wdat'][aa]['b2npc_fnays'][0][1:97, 1:37]
                hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                tor_area = np.multiply(hz, hy)
                pol_area = np.multiply(hz, hx)
                fnax = np.divide(fnaxs, tor_area)
                fnay = np.divide(fnays, tor_area)

                    
                    
                # fx_list = fhixn[:, 18]
                if aa == 'org':
                    fyi = sum(fnays[24:72, 3])/sum(pol_area[24:72, 3])
                elif aa == 'dot3':
                    fyi = sum(fnays[24:72, 0])/sum(pol_area[24:72, 0])
                else:
                    fyi = sum(fnays[24:72, 2])/sum(pol_area[24:72, 2])
                fyb = sum(fnays[:, -1])/sum(pol_area[:, -1])
                fxti = sum(abs(fnaxs[0, 18:]))/sum(tor_area[0, 18:])
                fxto = sum(fnaxs[-1, 18:])/sum(tor_area[-1, 18:])
                
                
                vlist.append([fyi, fyb, fxti, fxto])
                gl_list.append('A = {}'.format(A_dic[aa]))
                cl_list.append('{}'.format(color_dic[aa]))
            
            return vlist, gl_list, cl_list

        # Data
        
        vlist, gl_list, cl_list = return_value()
        
        categories = ['core input', 'SOL boundary', 'inner target', 'outer target']  # Categories
        values = vlist
        group_labels = gl_list
        colors = cl_list

        # Number of categories and groups
        n_categories = len(categories)
        n_groups = len(values)

        # X positions for the categories
        x = np.arange(n_categories)

        # Bar width
        bar_width = 0.2

        # Create the plot
        plt.figure()

        # Loop through each group to plot
        for i, group_values in enumerate(values):
            plt.bar(
                x + i * bar_width,  # Shift each group by bar_width
                group_values,       # Heights of the bars
                width=bar_width,    # Width of the bars
                label=group_labels[i],  # Label for the group
                color=colors[i]     # Color for the group
            )

        # Add labels, title, and legend
        plt.xlabel('Categories')
        # plt.ylabel('Values')
        plt.yscale('log')
        plt.title('Averaged ion particle flux (W/$m^{2}$) value Bar Chart')
        plt.xticks(x + bar_width * (n_groups - 1) / 2, categories)  # Center the category labels
        plt.legend()

        # Show the plot
        plt.tight_layout()
        plt.show()
        

    
    
    
    def polneteNG(self, side):
        

            
        fig, axs = plt.subplots(3, 1)
            
            
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                      'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        cspf_text = AnchoredText('{}'.format('(a) Poloidal flux $\Gamma_{\Theta}$ [$m^{-2} s^{-1}$]'), 
                                      loc='upper center')
        
        ane_text = AnchoredText('{}'.format('(b) Electron density [$m^{-3}$]'), 
                                      loc='upper center')
        
        te_text = AnchoredText('{}'.format('(c) Electron temperature [eV]'), 
                                      loc='upper center')
        

        for aa in self.data['dircomp']['multi_shift']:
            
            fnax = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
            hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
            hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
            tor_area = np.multiply(hz, hy)
            fnaxs = np.divide(fnax, tor_area)
            
            
            ne = self.data['b2fstate'][aa]['na'][1:97, 1:37, 1]
            
            Te_J = self.data['b2fstate'][aa]['te'][1:97, 1:37]
            ev = 1.6021766339999999 * pow(10, -19)
            te_pro = Te_J / ev
            
            
            
            
            if side == 'HFS':
                
                index_list = np.linspace(0, 25, 26)
            
            elif side == 'LFS':
                
                index_list = np.linspace(60, 95, 36)
                
                
            fx_list = []
            te_list = []
            ne_list = []
            fxt_list = []
            
            for ii in index_list:
                
                
                fx_list.append(abs(fnaxs[int(ii), 18]))
                ne_list.append(ne[int(ii), 18])
                te_list.append(te_pro[int(ii), 18])
                # nd_list.append(sum(nd[int(ii), 18:]))
            
            axs[0].plot(index_list, fx_list, color = color_dic[aa], 
                          label = 'A = {}'.format(A_dic[aa]))
            axs[1].plot(index_list, ne_list, color = color_dic[aa])
            axs[2].plot(index_list, te_list, color = color_dic[aa])
                
            
                
                

                
        if side == 'HFS':
            
            axs[0].axvline(x= 0,ls = '--', color = 'black')
            axs[0].axvline(x= 24,ls = '--', color = 'gray')
            axs[0].axvline(x= 25,ls = '--', color = 'brown')
            axs[1].axvline(x= 0, ls= '--', color = 'black', label = 'inner target')
            axs[1].axvline(x= 24,ls = '--', color = 'gray', label = 'x point')
            axs[1].axvline(x= 25,ls = '--', color = 'brown', label = 'poloidal angle 250')
            axs[2].axvline(x= 0, ls= '--', color = 'black')
            axs[2].axvline(x= 24,ls = '--', color = 'gray')
            axs[2].axvline(x= 25,ls = '--', color = 'brown')
        
        elif side == 'LFS':
            
            axs[0].axvline(x= 95,ls = '--', color = 'black')
            axs[0].axvline(x= 72,ls = '--', color = 'gray')
            axs[0].axvline(x= 60,ls = '--', color = 'brown')
            axs[1].axvline(x= 95, ls= '--', color = 'black', label = 'outer target')
            axs[1].axvline(x= 72,ls = '--', color = 'gray', label = 'x point')
            axs[1].axvline(x= 60,ls = '--', color = 'brown', label = 'outer midplane')
            axs[2].axvline(x= 95, ls= '--', color = 'black')
            axs[2].axvline(x= 72,ls = '--', color = 'gray')
            axs[2].axvline(x= 60,ls = '--', color = 'brown')
            
        
        axs[1].set_yscale('log')
        # axs[2].set_yscale('log')
        axs[0].legend(loc= 'upper right')
        axs[1].legend(loc = 'upper right')
        axs[0].add_artist(cspf_text)
        axs[1].add_artist(ane_text)
        axs[2].add_artist(te_text)
        axs[2].set_xlabel('poloidal index')
            
            
        plt.subplots_adjust(hspace=.0)
        plt.suptitle('{}'.format(side))
    
    
    
    def totfluxesNG(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            
            fig, axs = plt.subplots(3, 1)
            
            
            
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
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                fx_list = []
                s_list = []
                nd_list = []
                
                for ii in pol_list:
                    
                    fx_list.append(fnaxs[int(ii), 18])
                    s_list.append(source[int(ii), 18])
                    nd_list.append(neuden[int(ii), 18])
                
                
                    
                axs[0].plot(ang_list, fx_list, color = color_dic[aa])
                axs[1].plot(ang_list, s_list, color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                axs[2].plot(ang_list, nd_list, color = color_dic[aa])
                    
            
            # axs[0].set_xlabel('poloidal angle')
            # axs[1].set_xlabel('poloidal angle')
            axs[0].axvline(x= 0,ls = '--', color = 'black', label = 'LFS')
            axs[0].axvline(x= 180,ls = '--', color = 'brown', label = 'HFS')
            axs[1].axvline(x= 0, ls= '--', color = 'black')
            axs[1].axvline(x= 180,ls = '--', color = 'brown')
            axs[2].axvline(x= 0, ls= '--', color = 'black')
            axs[2].axvline(x= 180,ls = '--', color = 'brown')
            axs[0].axhline(y= 0, ls = '--', color = 'gray', label= '$\Gamma_{\Theta}$ = 0')
            axs[0].legend(loc = 'upper right')
            axs[1].legend(loc= 'upper right')
            axs[2].set_xlabel('poloidal angle')
            axs[1].set_yscale('log')
            axs[2].set_yscale('log')
            plt.subplots_adjust(hspace=.0)
    
    
    
    def rebu_fluxesNG(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            
            fig, axs = plt.subplots()
            
            
            for aa in self.data['dircomp']['multi_shift']:
                
                fnax = self.data['b2wdat'][aa]['b2nph9_fhix'][1:97, 1:37]
                hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                tor_area = np.multiply(hz, hy)
                fnaxs = np.divide(fnax, tor_area)
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                psiN = self.data['psi']['psival'][aa][1:37, 1:97]
                std_psiN = psiN[18, st:ed]
                h_psiN = max(std_psiN)
                l_psiN = min(std_psiN)
                
                print('{} case first grid outside of separatrix psiN value max is: {}, and min is: {}'.format(aa, h_psiN, l_psiN))
                
                
                fx_list = []
                
                for ii in pol_list:
                    
                    fx_list.append(fnaxs[int(ii), 18])
                
                
                    
                axs.plot(ang_list, fx_list, color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                    
            
            # axs[0].set_xlabel('poloidal angle')
            # axs[1].set_xlabel('poloidal angle')
            axs.axvline(x= 0,ls = '--', color = 'black', label = 'LFS')
            axs.axvline(x= 180,ls = '--', color = 'brown', label = 'HFS')
            axs.axhline(y= 0, ls = '--', color = 'gray', label= '$q_{\Theta}$ = 0')
            axs.legend(loc = 'upper right')
            axs.set_xlabel('poloidal angle')
    
        
    
    def totfluxesNR(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            
            fig, axs = plt.subplots(3, 1)
            
            
            
            for aa in self.data['dircomp']['multi_shift']:
                
                fnax = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
                fnaxs = np.divide(fnax, hz)
                
                
                s_term = self.data['b2wdat'][aa]['b2npc_sna'][0][1:97, 1:37]
                source = np.divide(s_term, hz)
                
                neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                area = np.multiply(hx, hy)
                nd = np.multiply(neuden, area)
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                fx_list = []
                s_list = []
                nd_list = []
                
                for ii in pol_list:
                    
                    fx_list.append(sum(fnaxs[int(ii), 18:]))
                    s_list.append(sum(source[int(ii), 18:]))
                    nd_list.append(sum(nd[int(ii), 18:]))
                
                
                    
                axs[0].plot(ang_list, fx_list, color = color_dic[aa])
                axs[1].plot(ang_list, s_list, color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                axs[2].plot(ang_list, nd_list, color = color_dic[aa])
                    
            
            # axs[0].set_xlabel('poloidal angle')
            # axs[1].set_xlabel('poloidal angle')
            axs[0].axvline(x= 0,ls = '--', color = 'black', label = 'LFS')
            axs[0].axvline(x= 180,ls = '--', color = 'brown', label = 'HFS')
            axs[1].axvline(x= 0, ls= '--', color = 'black')
            axs[1].axvline(x= 180,ls = '--', color = 'brown')
            axs[2].axvline(x= 0, ls= '--', color = 'black')
            axs[2].axvline(x= 180,ls = '--', color = 'brown')
            axs[0].axhline(y= 0, ls = '--', color = 'gray', label= '$\Gamma_{\Theta}$ = 0')
            axs[0].legend(loc = 'upper right')
            axs[1].legend(loc= 'upper right')
            axs[2].set_xlabel('poloidal angle')
            axs[1].set_yscale('log')
            axs[2].set_yscale('log')
            plt.subplots_adjust(hspace=.0)
    
    
    def avgfluxesNG(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            
            fig, axs = plt.subplots(2, 1)
            
            
            
            for aa in self.data['dircomp']['multi_shift']:
                
                fnax = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
                hy = self.data['b2wdat'][aa]['hy'][1:97, 1:37]
                g_coe = np.multiply(hz, hy)
                fnaxs = np.divide(fnax, g_coe)
                
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                fx_list = []
                hy_list = []
                
                for ii in pol_list:
                    
                    fx_list.append(sum(fnaxs[int(ii), 18:])/len(fnaxs[int(ii), 18:]))
                    hy_list.append(sum(hy[int(ii), 18:])/len(hy[int(ii), 18:]))
                
                
                    
                axs[0].plot(ang_list, fx_list, color = color_dic[aa])
                axs[1].plot(ang_list, hy_list, color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
                    
            
            # axs[0].set_xlabel('poloidal angle')
            # axs[1].set_xlabel('poloidal angle')
            axs[0].axvline(x= 0,ls = '--', color = 'black', label = 'LFS')
            axs[0].axvline(x= 180,ls = '--', color = 'brown', label = 'HFS')
            axs[1].axvline(x= 0, ls= '--', color = 'black')
            axs[1].axvline(x= 180,ls = '--', color = 'brown')
            axs[0].axhline(y= 0, ls = '--', color = 'gray', label= '$\Gamma_{\Theta}$ = 0')
            axs[0].legend(loc = 'upper right')
            axs[1].legend(loc= 'upper right')
            axs[1].set_xlabel('poloidal angle')
            plt.subplots_adjust(hspace=.0)
    
    
    
    
    
    def totsource(self):
        
            
        fig, axs = plt.subplots(2, 1)
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        HFS_text = AnchoredText('{}'.format('(a) Total source [$s^{-1}$] at HFS'), 
                                     loc='upper center')
        
        LFS_text = AnchoredText('{}'.format('(b) Total source [$s^{-1}$] at LFS'), 
                                     loc='upper center')
        
        for side in ['HFS', 'LFS']:
        
        
            for aa in self.data['dircomp']['multi_shift']:
                
                source = self.data['b2wdat'][aa]['b2npc_sna'][0][1:97, 1:37]
                
                # hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
                # fnaxs = np.divide(fnax, hz)
                       
                if side == 'HFS':
                    
                    index_list = np.linspace(0, 25, 26)
                    
                    fx_list = []
                    
                    for ii in index_list:
                        
                        
                        fx_list.append(abs(sum(source[int(ii), 18:])))
                    
                    axs[0].plot(index_list, fx_list, color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))
                    
                    
                elif side == 'LFS':
                    
                    index_list = np.linspace(60, 95, 36)
                    
                    fx_list = []
                    
                    for ii in index_list:
                        
                        fx_list.append(sum(source[int(ii), 18:]))
                    
                    axs[1].plot(index_list, fx_list, color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))

                

        axs[0].axvline(x= 0,ls = '--', color = 'black', label = 'inner target')
        axs[0].axvline(x= 25,ls = '--', color = 'brown', label = 'poloidal angle 250')
        axs[1].axvline(x= 95, ls= '--', color = 'black', label = 'outer target')
        axs[1].axvline(x= 60,ls = '--', color = 'brown', label = 'outer midplane')
        axs[0].axvline(x= 24, ls= '--', color = 'gray', label = 'left cut')
        axs[1].axvline(x= 72,ls = '--', color = 'gray', label = 'right cut')
        axs[0].set_yscale('log')
        axs[1].set_yscale('log')
        axs[1].legend(loc= 'upper right')
        axs[0].legend(loc = 'upper right')
        axs[0].add_artist(HFS_text)
        axs[1].add_artist(LFS_text)
        axs[0].set_xlabel('poloidal index')
        axs[1].set_xlabel('poloidal index')
    
    
    def totnd(self):
        
            
        fig, axs = plt.subplots(2, 1)
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        HFS_text = AnchoredText('{}'.format('(a) Total neuden at HFS'), 
                                     loc='upper center')
        
        LFS_text = AnchoredText('{}'.format('(b) Total neuden at LFS'), 
                                     loc='upper center')
        
        for side in ['HFS', 'LFS']:
        
        
            for aa in self.data['dircomp']['multi_shift']:
                
                neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                nd = np.multiply(neuden, vol)
                       
                if side == 'HFS':
                    
                    index_list = np.linspace(0, 25, 26)
                    
                    nd_list = []
                    
                    for ii in index_list:
                        
                        
                        nd_list.append(abs(sum(nd[int(ii), 18:])))
                    
                    axs[0].plot(index_list, nd_list, color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))
                    
                    
                elif side == 'LFS':
                    
                    index_list = np.linspace(60, 95, 36)
                    
                    nd_list = []
                    
                    for ii in index_list:
                        
                        nd_list.append(sum(nd[int(ii), 18:]))
                    
                    axs[1].plot(index_list, nd_list, color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))

                

        axs[0].axvline(x= 0,ls = '--', color = 'black', label = 'inner target')
        axs[0].axvline(x= 25,ls = '--', color = 'brown', label = 'poloidal angle 250')
        axs[1].axvline(x= 95, ls= '--', color = 'black', label = 'outer target')
        axs[1].axvline(x= 60,ls = '--', color = 'brown', label = 'outer midplane')
        axs[0].axvline(x= 24, ls= '--', color = 'gray', label = 'left cut')
        axs[1].axvline(x= 72,ls = '--', color = 'gray', label = 'right cut')
        axs[0].set_yscale('log')
        axs[1].set_yscale('log')
        axs[1].legend(loc= 'upper right')
        axs[0].legend(loc = 'upper right')
        axs[0].add_artist(HFS_text)
        axs[1].add_artist(LFS_text)
        axs[0].set_xlabel('poloidal index')
        axs[1].set_xlabel('poloidal index')
    
    
    
    def totndNR(self):
        
            
        fig, axs = plt.subplots(2, 1)
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        HFS_text = AnchoredText('{}'.format('(a) Total area neuden [$m^{-1}$] at HFS'), 
                                     loc='upper center')
        
        LFS_text = AnchoredText('{}'.format('(b) Total area neuden [$m^{-1}$] at LFS'), 
                                     loc='upper center')
        
        for side in ['HFS', 'LFS']:
        
        
            for aa in self.data['dircomp']['multi_shift']:
                
                neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
                vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
                area = np.divide(vol, hz)
                nd = np.multiply(neuden, area)
                       
                if side == 'HFS':
                    
                    index_list = np.linspace(0, 25, 26)
                    
                    nd_list = []
                    
                    for ii in index_list:
                        
                        
                        nd_list.append(abs(sum(nd[int(ii), 18:])))
                    
                    axs[0].plot(index_list, nd_list, color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))
                    
                    
                elif side == 'LFS':
                    
                    index_list = np.linspace(60, 95, 36)
                    
                    nd_list = []
                    
                    for ii in index_list:
                        
                        nd_list.append(sum(nd[int(ii), 18:]))
                    
                    axs[1].plot(index_list, nd_list, color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))

                

    
        axs[0].axvline(x= 0,ls = '--', color = 'black', label = 'inner target')
        axs[0].axvline(x= 25,ls = '--', color = 'brown', label = 'poloidal angle 250')
        axs[1].axvline(x= 95, ls= '--', color = 'black', label = 'outer target')
        axs[1].axvline(x= 60,ls = '--', color = 'brown', label = 'outer midplane')
        axs[0].axvline(x= 24, ls= '--', color = 'gray', label = 'left cut')
        axs[1].axvline(x= 72,ls = '--', color = 'gray', label = 'right cut')
        axs[0].set_yscale('log')
        axs[1].set_yscale('log')
        axs[1].legend(loc= 'upper right')
        axs[0].legend(loc = 'upper right')
        axs[0].add_artist(HFS_text)
        axs[1].add_artist(LFS_text)
        axs[0].set_xlabel('poloidal index')
        axs[1].set_xlabel('poloidal index')
    
    
    
    
    
    def totsourceNR(self):
        
            
        fig, axs = plt.subplots(2, 1)
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        HFS_text = AnchoredText('{}'.format('(a) Total area source [$m^{-1}*s^{-1}$] at HFS'), 
                                     loc='upper center')
        
        LFS_text = AnchoredText('{}'.format('(b) Total area source [$m^{-1}*s^{-1}$] at LFS'), 
                                     loc='upper center')
        
        for side in ['HFS', 'LFS']:
        
        
            for aa in self.data['dircomp']['multi_shift']:
                
                source = self.data['b2wdat'][aa]['b2npc_sna'][0][1:97, 1:37]
                
                hz = self.data['b2wdat'][aa]['hz'][1:97, 1:37]
                sx = np.divide(source, hz)
                       
                if side == 'HFS':
                    
                    index_list = np.linspace(0, 25, 26)
                    
                    sx_list = []
                    
                    for ii in index_list:
                        
                        
                        sx_list.append(abs(sum(sx[int(ii), 18:])))
                    
                    axs[0].plot(index_list, sx_list, color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))
                    
                    
                elif side == 'LFS':
                    
                    index_list = np.linspace(60, 95, 36)
                    
                    sx_list = []
                    
                    for ii in index_list:
                        
                        sx_list.append(sum(sx[int(ii), 18:]))
                    
                    axs[1].plot(index_list, sx_list, color = color_dic[aa], 
                                 label = 'A = {}'.format(A_dic[aa]))

                

        axs[0].axvline(x= 0,ls = '--', color = 'black', label = 'inner target')
        axs[0].axvline(x= 25,ls = '--', color = 'brown', label = 'poloidal angle 250')
        axs[1].axvline(x= 95, ls= '--', color = 'black', label = 'outer target')
        axs[1].axvline(x= 60,ls = '--', color = 'brown', label = 'outer midplane')
        axs[0].axvline(x= 24, ls= '--', color = 'gray', label = 'left cut')
        axs[1].axvline(x= 72,ls = '--', color = 'gray', label = 'right cut')
        axs[0].set_yscale('log')
        axs[1].set_yscale('log')
        axs[1].legend(loc= 'upper right')
        axs[0].legend(loc = 'upper right')
        axs[0].add_artist(HFS_text)
        axs[1].add_artist(LFS_text)
        axs[0].set_xlabel('poloidal index')
        axs[1].set_xlabel('poloidal index')
        
        
    
    
