# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 18:24:05 2023

@author: user
"""

from B2plotter_plot import Opacity_study
import matplotlib.pyplot as plt
import matplotlib.cm as pcm
import fitting_method as fm
import numpy as np




class PlotContour(Opacity_study):
    def __init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters, Publish):
        Opacity_study.__init__(self, DEV, withshift, withseries, 
                               DefaultSettings, loadDS, Parameters, Publish)

    def flux_expansion_map(self, pol_loc, iter_index):
        
        if self.withshift == False and self.withseries == False:
            
            RR_sep = self.data['midplane_calc']['R_Rsep']
            flux_expand_map = np.zeros([self.data['b2fgeo']['ny'], self.data['b2fgeo']['nx']])
            
            
            for pol_loc in range(self.data['b2fgeo']['nx']):
                self.calc_dsa(pol_loc)
                arcR = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]
                flux_fit_dic = fm.flux_expand_fit(RRsep = RR_sep, arclength = arcR)
                
                flux_expand = flux_fit_dic['flux_fitcoe'][0]
                a_flux_exp = flux_expand*np.ones(self.data['b2fgeo']['ny'])
                flux_expand_map[:, pol_loc] = a_flux_exp
                
            
            Attempt = self.data['dircomp']['Attempt']
            DRT = self.data['dirdata']['outputdir']['Output']
            XDIM = self.data['b2fgeo']['nx'] + 2
            YDIM = self.data['b2fgeo']['ny'] + 2
            
            
            RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                        usecols = (3)).reshape((YDIM, XDIM))
            VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                          usecols = (3)).reshape((YDIM,XDIM))
            
            R_con = RadLoc[1:37, 1:97]
            Z_con = VertLoc[1:37, 1:97]
            
            contour_dic = {'R_coord': R_con, 'Z_coord': Z_con, 
                           'flux_map': flux_expand_map}
            self.data['flux_contour'] = contour_dic
            
            
            map_flat = flux_expand_map.flatten()
            
            
            CMAP = pcm.viridis
            NORM= plt.Normalize(map_flat.min(), map_flat.max())
            
            plt.figure(figsize=(6,12))
            plt.contourf(R_con, Z_con, flux_expand_map, levels= 20, cmap=CMAP,norm=NORM)
            
            SM= pcm.ScalarMappable(NORM,CMAP)    
            plt.colorbar(SM)
            
            
            return flux_expand
        
        elif self.withshift == True and self.withseries == False:
            
            arcR = self.data['dsa']['dsa_{}'.format(pol_loc)][iter_index]['dsa_{}_val'.format(pol_loc)]
            RR_sep = self.data['midplane_calc'][iter_index]['R_Rsep']
            
            arcR_inv = list(reversed(arcR))
            RRsep_inv = list(reversed(RR_sep))
                                              
            flux_fit_dic = fm.flux_expand_fit(RRsep = RR_sep, arclength = arcR)
            
            flux_expand = flux_fit_dic['flux_fitcoe'][0]
            
            return flux_expand
        
        elif self.withshift == False and self.withseries == True:
            
            arcR = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]
            RR_sep = self.data['midplane_calc']['R_Rsep']
            
            arcR_inv = list(reversed(arcR))
            RRsep_inv = list(reversed(RR_sep))
            
            flux_fit_dic = fm.flux_expand_fit(RRsep = RR_sep, arclength = arcR)
            
            flux_expand = flux_fit_dic['flux_fitcoe'][0]
            
            return flux_expand
        
        elif self.withshift == True and self.withseries == True:
            print('calc_flux_expansion is not there yet, to be continue...')
            
        else:
            print('There is a bug')

    def load_vessel(self):
        # try:
        #     WallFile = np.loadtxt('{}/mesh.extra'.format(self.data['dirdata']['tbase']))
        # except:
        #     print('mesh.extra file not found! Using vvfile.ogr instead')
        #     WallFile=None
            
        VVFILE = np.loadtxt('{}/vvfile.ogr'.format(self.data['dirdata']['simutop']))
        
        self.data['vessel'] = VVFILE
        
        # if plot:
        #     plt.plot
    
    def plot_all_radial(self):
        
        if self.data['outputdata']['Ne'].all() == None or self.data['outputdata']['Te'].all() == None:
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
        else:
            pass
        
        ne_pro = self.data['outputdata']['Ne']
        te_pro = self.data['outputdata']['Te']
        
        core_ne_pro = ne_pro[:, 25:71]
        core_te_pro = te_pro[:, 25:71]
        
    
        mean_core_ne = np.mean(core_ne_pro, axis=1)
        std_core_ne = np.std(core_ne_pro, axis=1)
        # print(std_core_ne)
        mean_core_te = np.mean(core_te_pro, axis=1)
        std_core_te = np.std(core_te_pro, axis=1)
        
        
        psiN = self.data['experimental_fit']['psiN']
        ne = self.data['experimental_fit']['ne']*pow(10, 20)
        te = self.data['experimental_fit']['te']*pow(10, 3)
        
        
        plt.figure(figsize=(7,7))
        plt.yscale('log')
        plt.errorbar(psiN, mean_core_ne, yerr= std_core_ne, fmt = 'o', color = 'g', label= 'ne_solps')
        plt.plot(psiN, ne, 'o', color = 'r', label= 'ne_exp_fit')
        plt.xlabel('psiN')
        plt.title('electron density with experimental fit')
        plt.legend()
        
        
        plt.figure(figsize=(7,7))
        plt.yscale('log')
        plt.errorbar(psiN, mean_core_te, yerr= std_core_te, fmt = 'o', color = 'g', label= 'te_solps')
        plt.plot(psiN, te, 'o', color = 'r', label= 'te_exp_fit')
        plt.xlabel('psiN')
        plt.title('electron temperature with experimental fit')
        plt.legend()
        
        
    