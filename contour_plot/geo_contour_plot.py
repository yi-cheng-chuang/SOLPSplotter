# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 19:20:00 2025

@author: ychuang
"""

import numpy as np
from fit_data.fitting_method import fit_method_collection
from contour_plot.contourplot_toolbox import contour_plot_method_collect
from midplane_data.SOLPSplotter_PRmap import RP_mapping


class fluxexpansion_contour:
    
    
    
    def __init__(self, DF, data, fmc: fit_method_collection, cpmc: contour_plot_method_collect, rp: RP_mapping):
        
        self.DF = DF
        self.data = data
        self.fmc = fmc
        self.cpmc = cpmc
        self.rp = rp

    

    def calc_flux_expansion_line_method(self, RR_sep, arcR):
                   
        flux_fit_dic = self.fmc.flux_expand_fit(RRsep = RR_sep, arclength = arcR)
        
        flux_expand = flux_fit_dic['flux_fitcoe'][0]
        a_flux_exp = flux_expand*np.ones(self.data['b2fgeo']['ny'])
        
        return a_flux_exp


    def flux_expansion_contour_plot_method(self, RR_sep, flux_expand_map, itername):
        
        for pol_loc in range(self.data['b2fgeo']['nx']):
            
            if itername == None:
                arcR = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]

            else:
                arcR = self.data['dsa']['dsa_{}'.format(pol_loc)][itername]['dsa_{}_val'.format(pol_loc)]
            
            
            flux_fit_dic = self.fmc.flux_expand_fit(RRsep = RR_sep, arclength = arcR)
            
            flux_expand = flux_fit_dic['flux_fitcoe'][0]
            a_flux_exp = flux_expand*np.ones(self.data['b2fgeo']['ny'])
            
            flux_expand_map[:, pol_loc] = a_flux_exp
        
        
        RadLoc = self.data['grid']['RadLoc']
        VertLoc = self.data['grid']['VertLoc']
        
        R_con = RadLoc[1:37, 1:97]
        Z_con = VertLoc[1:37, 1:97]
        
        # contour_dic = {'R_coord': R_con, 'Z_coord': Z_con, 
        #                'flux_map': flux_expand_map}
        
        contour_dic = {'flux_map': flux_expand_map}
        
        # map_flat = flux_expand_map.flatten()
        
        
        self.cpmc.contour_plot(plot_2dval = flux_expand_map, R_coord = R_con, 
                         Z_coord = Z_con, quantity = 'flux expansion')
        
        
        # self.plot_vessel(itername = itername, independent = False, meter = True)
        # plt.show()
        
        return contour_dic
        
       
    def flux_expansion_contour_plot(self):
        
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        
        if withshift == False and withseries == False:
            
            for pol_loc in range(self.data['b2fgeo']['nx']):
                self.calc_dsa(pol_loc)
            
                        
            RR_sep = self.data['midplane_calc']['R_Rsep']
            flux_expand_map = np.zeros([self.data['b2fgeo']['ny'], self.data['b2fgeo']['nx']])
            
            contour_dic = self.flux_expansion_contour_plot_method(RR_sep = RR_sep, 
                            flux_expand_map = flux_expand_map, itername = None)
            
            self.data['flux_contour'] = contour_dic
            
            
        
        elif withshift == True and withseries == False:
            
            contour_dic = {}
            
            for pol_loc in range(self.data['b2fgeo']['nx']):
                self.calc_dsa(pol_loc)
            
            for aa in self.data['dircomp']['multi_shift']:
            
                RR_sep = self.data['midplane_calc'][aa]['R_Rsep']
                flux_expand_map = np.zeros([self.data['b2fgeo']['ny'], self.data['b2fgeo']['nx']])
            
            
                contour_dic[aa] = self.flux_expansion_contour_plot_method(RR_sep = RR_sep, 
                                flux_expand_map = flux_expand_map, itername = aa)
                
            
            self.data['flux_contour'] = contour_dic
            
        
        elif withshift == False and withseries == True:
            
            for pol_loc in range(self.data['b2fgeo']['nx']):
                self.rp.calc_dsa(pol_loc)
                                    
            RR_sep = self.data['midplane_calc']['R_Rsep']
            flux_expand_map = np.zeros([self.data['b2fgeo']['ny'], self.data['b2fgeo']['nx']])
            
            contour_dic = self.flux_expansion_contour_plot_method(RR_sep = RR_sep, 
                            flux_expand_map = flux_expand_map, itername = None)
            
            self.data['flux_contour'] = contour_dic
        
        elif withshift == True and withseries == True:
            print('calc_flux_expansion is not there yet, to be continue...')
            
        else:
            print('There is a bug')