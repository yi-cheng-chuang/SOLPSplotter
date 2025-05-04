# -*- coding: utf-8 -*-
"""
Created on Thu May  1 16:21:58 2025

@author: ychuang
"""



import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnchoredText


class result_explain:
    
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    def plot_addition_ndop(self):
        
        
        
        # fig, axs = plt.subplots(3, 1)
        
        fig, axs = plt.subplots(2, 1, sharex = True)  # Default: no shared x
        

        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        nd_text = AnchoredText('{}'.format('(c) Atomic neutral density $n_0^{pf}$ [$m^{-3}$]'), 
                                     loc='upper center')
        
        op_text = AnchoredText('{}'.format('(b) Neutral opaqueness'), loc='upper center')
        
        
        for aa in self.data['dircomp']['multi_shift']:
            
            nx = self.data['b2fgeo'][aa]['nx']
            ny = self.data['b2fgeo'][aa]['ny']
            
            fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:nx+1, 1:ny+1]
            
            kk = len(fnaxs)
            int_list = list(range(0, kk))
            
            fnaxs_list = abs(fnaxs[:, 18][::-1])
            
            neuden_fit = self.data['opacity_poloidal'][aa]['neutral_density']
            op_fit = self.data['opacity_poloidal'][aa]['dimensionless_opaqueness']
            
            neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
            nd_list = neuden[:, 18]
            
            ang_list = self.data['angle']['angle_list'][aa]
            arc_length = self.data['sepl'][aa]['sepl']
            
            
            axs[0].plot(ang_list, op_fit, color = color_dic[aa], 
                        label = 'A = {}'.format(A_dic[aa]))
            axs[1].plot(ang_list, neuden_fit, color = color_dic[aa])
            
            
        
            
        
        
        axs[1].set_yscale('log')
        axs[1].set_xlabel('poloidal angle')
        axs[0].axvline(x= 0, color='darkorange', linestyle='--', label = 'outer midplane')
        axs[0].axvline(x= 180, color='darkslategray', linestyle='--', label = 'inner midplane')
        axs[1].axvline(x= 0, color='darkorange', linestyle='--')
        axs[1].axvline(x= 180, color='darkslategray', linestyle='--')        
        axs[0].legend(loc= 'best')
        axs[1].add_artist(nd_text)
        axs[0].add_artist(op_text)
        
        plt.subplots_adjust(hspace=.0)
        

    
    
    def plot_addition_ndflux(self):
        
        
        
        fig, axs = plt.subplots(2, 1)
        

        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        nd_text = AnchoredText('{}'.format('(c) Atomic neutral density $n_0^{pf}$ [$m^{-3}$]'), 
                                     loc='upper center')
        
        gamma_text = AnchoredText('{}'.format('(d) Absolute value of surface integrated poloidal flux $\Gamma_{\Theta}$ [$s^{-1}$]'), 
                                     loc='upper center')
        
        
        for aa in self.data['dircomp']['multi_shift']:
            
            nx = self.data['b2fgeo'][aa]['nx']
            ny = self.data['b2fgeo'][aa]['ny']
            
            fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:nx+1, 1:ny+1]
            
            kk = len(fnaxs)
            int_list = list(range(0, kk))
            
            fnaxs_list = abs(fnaxs[:, 18][::-1])
            
            neuden_fit = self.data['opacity_poloidal'][aa]['neutral_density']
            op_fit = self.data['opacity_poloidal'][aa]['dimensionless_opaqueness']
            
            neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
            nd_list = neuden[:, 18]
            
            ang_list = self.data['angle']['angle_list'][aa]
            arc_length = self.data['sepl'][aa]['sepl']
            
            
            axs[0].plot(ang_list, neuden_fit, color = color_dic[aa], 
                        label = 'A = {}'.format(A_dic[aa]))
            axs[1].plot(arc_length, fnaxs_list, color = color_dic[aa])
            
            
        
            
        
        
        axs[0].set_yscale('log')
        axs[1].set_yscale('log')
        axs[0].set_xlabel('poloidal angle')
        axs[1].set_xlabel('outer to inner target arc length [m]')
        axs[0].axvline(x= 0, color='darkorange', linestyle='--', label = 'outer midplane')
        axs[0].axvline(x= 180, color='darkslategray', linestyle='--', label = 'inner midplane')        
        axs[1].axvline(x= arc_length[0], color='brown', linestyle='--', 
                       label= 'outer target')
        axs[1].axvline(x= arc_length[-1], color='gray', linestyle='--', 
                       label= 'inner target')
        axs[0].legend(loc= 'best')
        axs[1].legend(loc= 'lower left')
        axs[0].add_artist(nd_text)
        axs[1].add_artist(gamma_text)
    
    
    
    def shape_plot_addition(self, pol_list):
        
        aa = 'org'
        
        
        RadLoc = self.data['grid']['RadLoc'][aa]
        VertLoc = self.data['grid']['VertLoc'][aa]
        
        
        "shell"
        
        north_R = RadLoc[-1, :]
        north_Z = VertLoc[-1, :]
        
        south_R_core = RadLoc[0, 25:72]
        south_Z_core = VertLoc[0, 25:72]
        
        south_R_pfr1 = RadLoc[0, 73:]
        south_Z_pfr1 = VertLoc[0, 73:]
        
        south_R_pfr2 = RadLoc[0, :24]
        south_Z_pfr2 = VertLoc[0, :24]
        
        sep_R_up = RadLoc[18, :]
        sep_Z_up = VertLoc[18, :]
        
        sep_R_down = RadLoc[17, :]
        sep_Z_down = VertLoc[17, :]
        
        sep_R = 0.5*(sep_R_up + sep_R_down)
        sep_Z = 0.5*(sep_Z_up + sep_Z_down)
        
        gfile_data = self.data['gfile']['g']
        st_R = gfile_data['rmaxis']
        st_Z = gfile_data['zmaxis']
        
        ed_R = self.data['midplane_calc'][aa]['mid_R'][-1]
        ed_Z = self.data['midplane_calc'][aa]['mid_Z'][-1]
        
        outl = int(pol_list[-1])
        innl = int(pol_list[0])
        
        outl_R = RadLoc[:, outl]
        outl_Z = VertLoc[:, outl]
        
        
        innl_R = RadLoc[:, innl]
        innl_Z = VertLoc[:, innl]
        
        
        color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                     'dot7': 'blue'}
        
        A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                  'dot7': '2.8'}
        
        
        fig, axs = plt.subplots()
        
        "plot the shell"
        
        
        for ad in self.data['dircomp']['multi_shift']:
            
            vessel = self.data['vessel'][ad]
            if ad == 'org':
                
                axs.plot(vessel[:,0]/1000, vessel[:,1]/1000, color = color_dic[ad], 
                            label = 'A = {}'.format(A_dic[ad]))
            
            else:
                
                axs.plot(vessel[:,0]/1000, vessel[:,1]/1000, alpha = 0.4, color = color_dic[ad], 
                            label = 'A = {}'.format(A_dic[ad]))
                
            
        
        
        
        
        axs.plot(north_R, north_Z,'-', color = 'purple', label= 'grid boundary', linewidth = 5)
        axs.plot(south_R_pfr1, south_Z_pfr1,'-', color = 'purple',linewidth = 5)
        axs.plot(south_R_pfr2, south_Z_pfr2,'-', color = 'purple',linewidth = 5)
        axs.plot(south_R_core, south_Z_core,'-', color = 'purple',linewidth = 5)
        axs.plot([st_R, ed_R], [st_Z, ed_Z],'-', color = 'goldenrod', label = 'outer midplane')
        axs.plot(sep_R, sep_Z,'--', color = 'darkcyan', label = 'separatrix', linewidth = 2)
        axs.plot(outl_R, outl_Z,'-', color = 'black', label = 'angle limits')
        axs.plot(innl_R, innl_Z,'-', color = 'black')
        axs.legend(loc= 'upper right')

        
        axs.set_xlabel("R: [m]")
        axs.set_ylabel("Z: [m]")
        
        plt.gca().set_aspect('equal')
        
            
            
    
    
    
