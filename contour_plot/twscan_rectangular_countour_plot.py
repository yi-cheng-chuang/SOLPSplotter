# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 10:58:12 2025

@author: ychuang
"""

from matplotlib.offsetbox import AnchoredText
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker, colors
import matplotlib.tri as tri
from matplotlib.colors import LogNorm
from numpy import ma
from contour_plot.contourplot_toolbox import contour_plot_method_collect


class Rectangular_contour:
    
    
    
    def __init__(self, DF, data, cpmc: contour_plot_method_collect):
        
        self.DF = DF
        self.data = data
        self.cpmc = cpmc

    

    def twrectangular_test(self):
        
            

        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        
        psi_map = np.transpose(self.data['psi']['psival'])[1:nx+1, 1:ny+1]


        Ra = np.arange(1, nx+1)
        extend_Ra = np.tile(Ra, (ny, 1))
        Ra_map = np.transpose(extend_Ra)
        
        self.data['save_rectangular'] = Ra_map
        
        # pol_list = self.DF.poloidal_loc_list
        
        # print(pol_list)



    
    def twrectangular_contour_plot(self, scan_style, plot_name, norm_type, label_type, pol_loc_list):
        

        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        
        psi_map = np.transpose(self.data['psi']['psival'])[1:nx+1, 1:ny+1]

        


        ang_label = self.data['angle']['angle_list']
        extend_ang = np.tile(ang_label, (ny, 1))
        ang_map = np.transpose(extend_ang)
        
        
        Ra = np.arange(1, nx+1)
        extend_Ra = np.tile(Ra, (ny, 1))
        Ra_map = np.transpose(extend_Ra)
        
        rec_dic = {'angle': ang_map, 'index': Ra_map}
        
        self.data['save_rectangular'] = rec_dic
        
        
        RadLoc = np.transpose(self.data['grid']['RadLoc'])[1:nx + 1, 1:ny + 1]
        VertLoc = np.transpose(self.data['grid']['VertLoc'])[1:nx + 1, 1:ny + 1]
            
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries

        
        if withshift == False and withseries == True:
            
            
            
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
                    
                    # color_list = ['red', 'orange', 'green', 'blue', 'purple']
                    
                    # color_dic = self.pair_dic(keys = keylist_b, values = color_list)
                    

                    # print('scan_list after initial:')
                    # print(scan_list)
                    iterlist = []
                    
                    
                    for tb in keylist_b:
                        
                        if scan_style == 'tempscan':
                            
                            it_in = (ta, tb)
                        
                        elif scan_style == 'denscan':
                            
                            it_in = (tb, ta)
                        
                        else:
                            print('twinscan_plot_method, please check the scan_style!')
                        
                        iterlist.append(it_in)
                    
                    
                    
                    for aa in iterlist:
                        
                            
                        if scan_style == 'tempscan':
                            
                            ad = aa[1]
                            ap = aa[0]
                        
                        elif scan_style == 'denscan':
                            
                            ad = aa[0]
                            ap = aa[1]
                        
                        else:
                            print('twscan_contour_plot, please check scan_style')
                            
                        
                
                        if scan_style == 'denscan':
                            
                            fig, axs = plt.subplots()
                            
                            title_ap = float(ap)*pow(10, 5)
                            label_ad = float(ad)*pow(10, 20)
                            
                            nf = aa[0]
                            tf = aa[1]
                            
                            

                            b2fstate = self.data['b2fstate'][nf][tf]
                            ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
                            Te_J = b2fstate['te'][1:nx+1, 1:ny+1]
                            
                            ev = 1.6021766339999999 * pow(10, -19)
                            te_dat = Te_J / ev
                            
                            source = self.data['b2wdat'][nf][tf]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                            vol = self.data['b2wdat'][nf][tf]['vol'][1:nx+1, 1:ny+1]
                            sx = np.divide(source, vol)
                            
                            neuden_dat = self.data['ft44'][nf][tf]['dab2'][:, :, 0]
                            
                            
                            
                            if plot_name == 'Te':
                                
                                datname = 'Te'
                                title_name = '{0} $\Gamma_r ={1}$*$10^{{20}}$ 1/s, $q_r = {2}*10^5$ W'.format(datname, ad, ap)
                                
                                
                                if label_type == 'angle':
                                    
                                    st = int(pol_loc_list[0])
                                    ed = int(pol_loc_list[-1]) + 1
                                                                        
                                    te_limit = te_dat[st:ed, 14:]
                                    psi_limit = psi_map[st:ed, 14:]
                                    Ra_limit = ang_map[:, 14:]
                                    
                                    self.cpmc.twcontourp(plot_2dval = te_limit, R_coord = Ra_limit, 
                                            Z_coord = psi_limit, quantity = title_name, axs = axs, 
                                                      norm_type = norm_type, lv = 40)
                                    
                                    axs.set_ylabel("$\psi_N$")
                                    axs.set_xlabel("poloidal angle")
                                
                                elif label_type == 'index':
                                    
                                    te_limit = te_dat[:, 14:]
                                    psi_limit = psi_map[:, 14:]
                                    Ra_limit = Ra_map[:, 14:]
                                    
                                    
                                    self.cpmc.twcontourp(plot_2dval = te_limit, R_coord = Ra_limit, 
                                            Z_coord = psi_limit, quantity = title_name, axs = axs, 
                                                      norm_type = norm_type, lv = 40)
                                    
                                    axs.set_ylabel("$\psi_N$")
                                    axs.set_xlabel("poloidal index")

                                
                                
                                
                            
                            
                            elif plot_name == 'neuden':
                                
                                datname = 'Atomic neutral density'
                                title_name = '{0} $\Gamma_r ={1}$*$10^{{20}}$ 1/s, $q_r = {2}*10^5$ W'.format(datname, ad, ap)
                                
                                
                                if label_type == 'angle':
                                    
                                    st = int(pol_loc_list[0])
                                    ed = int(pol_loc_list[-1]) + 1
                                                                        
                                    sx_limit = neuden_dat[st:ed, 14:]
                                    psi_limit = psi_map[st:ed, 14:]
                                    Ra_limit = ang_map[:, 14:]
                                    
                                    self.cpmc.twcontourp(plot_2dval = sx_limit, R_coord = Ra_limit, 
                                            Z_coord = psi_limit, quantity = title_name, axs = axs, 
                                                      norm_type = norm_type, lv = 40)
                                    
                                    axs.set_ylabel("$\psi_N$")
                                    axs.set_xlabel("poloidal angle")
                                    
                                    
                                    
                                
                                elif label_type == 'index':
                                    
                                    sx_limit = neuden_dat[:, 14:]
                                    psi_limit = psi_map[:, 14:]
                                    Ra_limit = Ra_map[:, 14:]
                                    
                                    
                                    self.cpmc.twcontourp(plot_2dval = sx_limit, R_coord = Ra_limit, 
                                            Z_coord = psi_limit, quantity = title_name, axs = axs, 
                                                      norm_type = norm_type, lv = 40)
                                    
                                    axs.set_ylabel("$\psi_N$")
                                    axs.set_xlabel("poloidal index")
                                    
                                    
                                                                 
                                
                                
                            
                            
                            
                            
                            elif plot_name == 'sx':
                                
                                datname = 'Source'
                                title_name = '{0} $\Gamma_r ={1}$*$10^{{20}}$ 1/s, $q_r = {2}*10^5$ W'.format(datname, ad, ap)
                                
                                
                                if label_type == 'angle':
                                    
                                    st = int(pol_loc_list[0])
                                    ed = int(pol_loc_list[-1]) + 1
                                                                        
                                    sx_limit = sx[st:ed, 14:]
                                    psi_limit = psi_map[st:ed, 14:]
                                    Ra_limit = ang_map[:, 14:]
                                    
                                    self.cpmc.twcontourp(plot_2dval = sx_limit, R_coord = Ra_limit, 
                                            Z_coord = psi_limit, quantity = title_name, axs = axs, 
                                                      norm_type = norm_type, lv = 40)
                                    
                                    axs.set_ylabel("$\psi_N$")
                                    axs.set_xlabel("poloidal angle")
                                
                                elif label_type == 'index':
                                    
                                    sx_limit = sx[:, 14:]
                                    psi_limit = psi_map[:, 14:]
                                    Ra_limit = Ra_map[:, 14:]
                                    
                                    self.cpmc.twcontourp(plot_2dval = sx_limit, R_coord = Ra_limit, 
                                            Z_coord = psi_limit, quantity = title_name, axs = axs, 
                                                      norm_type = norm_type, lv = 40)
                                    
                                    axs.set_ylabel("$\psi_N$")
                                    axs.set_xlabel("poloidal index")
                                                                 
                                
                                
                            



             

                        elif scan_style == 'tempscan':
                            
                            title_ap = float(ap)*pow(10, 20)
                            label_ad = float(ad)*pow(10, 5)
                            # exp_an_fit = fit_dat['exp_fit']
                            # xcoord_cut = fit_dat['x_coord_cut']

                                
                            nf = aa[0]
                            tf = aa[1]
                            


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
                            
                            
                            self.contour_plot(plot_2dval = te_dat, R_coord = RadLoc, Z_coord = VertLoc, 
                                              quantity = f'{aa[0]:.2f/aa[1]:.2f} Electron temperature')
                                
                                
                        else:
                            print('neteTSplot_structure, please check the scan parameter')