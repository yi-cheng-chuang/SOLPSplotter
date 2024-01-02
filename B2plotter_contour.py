# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 18:24:05 2023

@author: user
"""

from B2plotter_plot import Opacity_study
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import matplotlib.tri as tri

import fitting_method as fm
import numpy as np




class PlotContour(Opacity_study):
    def __init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters, Publish):
        Opacity_study.__init__(self, DEV, withshift, withseries, 
                               DefaultSettings, loadDS, Parameters, Publish)
    
    
   
    def contour_plot(self, plot_2dval, R_coord, Z_coord, quantity):
        CMAP = cm.viridis
        NORM= plt.Normalize(plot_2dval.min(), plot_2dval.max())
        
        plt.figure(figsize=(6,12))
        plt.contourf(R_coord, Z_coord, plot_2dval, levels= 20, cmap=CMAP,norm=NORM)
        plt.title('{} contour plot'.format(quantity))
        
        
        SM= cm.ScalarMappable(NORM,CMAP)    
        plt.colorbar(SM)
    
    
   
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
            
            
            CMAP = cm.viridis
            NORM= plt.Normalize(map_flat.min(), map_flat.max())
            
            plt.figure(figsize=(6,12))
            plt.contourf(R_con, Z_con, flux_expand_map, levels= 20, cmap=CMAP,norm=NORM)
            plt.title('flux expansion contour plot')
            
            
            SM= cm.ScalarMappable(NORM,CMAP)    
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
            
        VVFILE = np.loadtxt('{}/baserun/vvfile.ogr'.format(self.data['dirdata']['simutop']))
        
        self.data['vessel'] = VVFILE
        
        # if plot:
        #     plt.plot
    
    def plot_all_radial(self):
        
        
        if self.withshift == False and self.withseries == False:
        
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
            
            core_ne_pro = ne_pro[:, 25:71]
            core_te_pro = te_pro[:, 25:71]
            core_neu_pro = neu_pro[:, 25:71]
            
            innerleg_ne = ne_pro[:, :25]
            innerleg_te = te_pro[:, :25]
            innerleg_neu = neu_pro[:, :25]
            
            outerleg_ne = ne_pro[:, 73:96]
            outerleg_te = te_pro[:, 73:96]
            outerleg_neu = neu_pro[:, 73:96]
            
        
            mean_core_ne = np.mean(core_ne_pro, axis=1)
            std_core_ne = np.std(core_ne_pro, axis=1)
            # print(std_core_ne)
            
            mean_core_te = np.mean(core_te_pro, axis=1)
            std_core_te = np.std(core_te_pro, axis=1)
            
            mean_core_neu = np.mean(core_neu_pro, axis=1)
            std_core_neu = np.std(core_neu_pro, axis=1)
            
            
            
            mean_innerleg_ne = np.mean(innerleg_ne, axis=1)
            std_innerleg_ne = np.std(innerleg_ne, axis=1)
            # print(std_innerleg_ne)
            mean_innerleg_te = np.mean(innerleg_te, axis=1)
            std_innerleg_te = np.std(innerleg_te, axis=1)
            
            mean_innerleg_neu = np.mean(innerleg_neu, axis=1)
            std_innerleg_neu = np.std(innerleg_neu, axis=1)
            
            
            
            mean_outerleg_ne = np.mean(outerleg_ne, axis=1)
            std_outerleg_ne = np.std(outerleg_ne, axis=1)
            # print(std_outerleg_ne)
            mean_outerleg_te = np.mean(outerleg_te, axis=1)
            std_outerleg_te = np.std(outerleg_te, axis=1)
            
            mean_outerleg_neu = np.mean(outerleg_neu, axis=1)
            std_outerleg_neu = np.std(outerleg_neu, axis=1)
            
            
            
            psiN = self.data['experimental_fit']['psiN']
            ne = self.data['experimental_fit']['ne']*pow(10, 20)
            te = self.data['experimental_fit']['te']*pow(10, 3)
            
            'core'
            
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
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN, mean_core_neu, yerr= std_core_neu, fmt = 'o', color = 'g', label= 'Neuden_solps')
            plt.xlabel('psiN')
            plt.title('Neutral density')
            plt.legend()
            
            'inner leg'
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN, mean_innerleg_ne, yerr= std_innerleg_ne, fmt = 'o', color = 'g', label= 'ne_solps')
            plt.xlabel('psiN')
            plt.title('inner leg electron density')
            plt.legend()
            
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN, mean_innerleg_te, yerr= std_innerleg_te, fmt = 'o', color = 'g', label= 'te_solps')
            plt.xlabel('psiN')
            plt.title('inner leg electron temperature')
            plt.legend()
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN, mean_innerleg_neu, yerr= std_innerleg_neu, fmt = 'o', color = 'g', label= 'neuden_solps')
            plt.xlabel('psiN')
            plt.title('inner leg neutral density')
            plt.legend()
            
            
            'outerleg'
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN, mean_outerleg_ne, yerr= std_outerleg_ne, fmt = 'o', color = 'g', label= 'ne_solps')
            plt.xlabel('psiN')
            plt.title('outer leg electron density')
            plt.legend()
            
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN, mean_outerleg_te, yerr= std_outerleg_te, fmt = 'o', color = 'g', label= 'te_solps')
            plt.xlabel('psiN')
            plt.title('outer leg electron temperature')
            plt.legend()
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN, mean_outerleg_neu, yerr= std_outerleg_neu, fmt = 'o', color = 'g', label= 'neuden_solps')
            plt.xlabel('psiN')
            plt.title('outer leg neutral density')
            plt.legend()
            
            
            
            'contour plot for ne & te'
            
            Attempt = self.data['dircomp']['Attempt']
            DRT = self.data['dirdata']['outputdir']['Output']
            Eirout = self.data['dirdata']['outputdir']['EirOutput']
            simudir = self.data['dirdata']['simudir']
            XDIM = self.data['b2fgeo']['nx'] + 2
            YDIM = self.data['b2fgeo']['ny'] + 2
            
            RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                        usecols = (3)).reshape((YDIM, XDIM))
            VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                          usecols = (3)).reshape((YDIM,XDIM))
            
            
            tz=np.loadtxt('{}/TriangVertLoc{}'.format(Eirout, str(Attempt)), 
                          usecols = (2))
                          
            tr=np.loadtxt('{}/TriangRadLoc{}'.format(Eirout, str(Attempt)), 
                          usecols = (2))
            
            Eiratom =np.loadtxt('{}/EirAtom{}'.format(Eirout, str(Attempt)), 
                          usecols = (2))
            
            
            CMAP = cm.viridis
            
            NORM_ne = plt.Normalize(ne_pro.min(), ne_pro.max())
            
            plt.figure(figsize=(6,12))
            plt.contourf(RadLoc, VertLoc, ne_pro, levels= 20, cmap=CMAP,norm=NORM_ne)
            plt.title('electron density contour plot')
            
            SM_ne= cm.ScalarMappable(NORM_ne,CMAP)    
            plt.colorbar(SM_ne)
            
            
            NORM_te = plt.Normalize(te_pro.min(), te_pro.max())
            
            plt.figure(figsize=(6,12))
            plt.contourf(RadLoc, VertLoc, te_pro, levels= 20, cmap=CMAP,norm=NORM_te)
            plt.title('electron temperature contour plot')
            
            SM_te= cm.ScalarMappable(NORM_te,CMAP)    
            plt.colorbar(SM_te)
            
            base_start_core = np.log10(np.nanmin(core_neu_pro))
            base_end_core = np.log10(np.nanmax(core_neu_pro))
            log_level_core = np.logspace(base_start_core, base_end_core, num=20, base= 10)
            print(log_level_core)
            
            NORM_neu_core = colors.LogNorm(np.nanmin(core_neu_pro), np.nanmax(core_neu_pro))
            
            
            plt.figure(figsize=(6,12))
            plt.contourf(RadLoc[:, 25:71], VertLoc[:, 25:71], core_neu_pro,
                         levels= log_level_core, cmap=CMAP, norm=NORM_neu_core)
            plt.title('Neutral density contour plot')
            
            SM_neu_core= cm.ScalarMappable(NORM_neu_core,CMAP)    
            plt.colorbar(SM_neu_core)
            
            base_start_inleg = np.log10(np.nanmin(innerleg_neu))
            base_end_inleg = np.log10(np.nanmax(innerleg_neu))
            log_level_inleg = np.logspace(base_start_inleg, base_end_inleg, num=20, base= 10)
            # print(log_level_inleg)
            
            NORM_neu_inleg = colors.LogNorm(np.nanmin(innerleg_neu), np.nanmax(innerleg_neu))
            
            
            plt.figure(figsize=(6,12))
            plt.contourf(RadLoc[:, :25], VertLoc[:, :25], innerleg_neu,
                         levels= log_level_inleg, cmap=CMAP, norm=NORM_neu_inleg)
            plt.title('Neutral density contour plot innerleg')
            
            SM_neu_inleg= cm.ScalarMappable(NORM_neu_inleg, CMAP)    
            plt.colorbar(SM_neu_inleg)
            
            base_start_outleg = np.log10(np.nanmin(outerleg_neu))
            base_end_outleg = np.log10(np.nanmax(outerleg_neu))
            log_level_outleg = np.logspace(base_start_outleg, base_end_outleg, num=20, base= 10)
            # print(log_level_outleg)
            
            NORM_neu_outleg = colors.LogNorm(np.nanmin(outerleg_neu), np.nanmax(outerleg_neu))
            
            
            plt.figure(figsize=(6,12))
            plt.contourf(RadLoc[:, 73:96], VertLoc[:, 73:96], outerleg_neu,
                         levels= log_level_outleg, cmap=CMAP, norm=NORM_neu_outleg)
            plt.title('Neutral density contour plot outerleg')
            
            SM_neu_outleg= cm.ScalarMappable(NORM_neu_outleg, CMAP)    
            plt.colorbar(SM_neu_outleg)
            
            
            
            base_start_eiratom = np.log10(np.nanmin(Eiratom))
            base_end_eiratom = np.log10(np.nanmax(Eiratom))
            log_level_eiratom = np.logspace(base_start_eiratom, base_end_eiratom, num=20, base= 10)
            # print(log_level_eiratom)
            
            NORM_neu_eiratom = colors.LogNorm(np.nanmin(Eiratom), np.nanmax(Eiratom))
            
            
            Nodes=np.fromfile('{}/fort.33'.format(simudir),sep=' ') #Alternatively use fort.33
            NN=int(Nodes[0])
            XNodes=Nodes[1:NN+1]
            YNodes=Nodes[NN+1:]
            
            
            numberlist = np.zeros(NN)
            for i in range(NN):
                numberlist[i] = i
            
            plt.figure(figsize=(7,7))
            plt.scatter(XNodes[:500], YNodes[:500])
            
            

            Triangles = np.loadtxt('{}/fort.34'.format(simudir), 
                skiprows=1, usecols=(1,2,3)) #Alternatively use fort.34
            # print(Triangles -1)

            TP = tri.Triangulation(XNodes, YNodes, triangles= (Triangles -1))
            
            
            plt.figure(figsize=(6,12))
            plt.tripcolor(TP, Eiratom, shading='flat', cmap= CMAP, norm= NORM_neu_eiratom)
            plt.title('Eiratom')
            # plt.title('Neutral density contour plot outerleg')
            
            SM_neu_eiratom= cm.ScalarMappable(NORM_neu_eiratom, CMAP)    
            plt.colorbar(SM_neu_eiratom)
            
            
            
        else:
            print('plot_all_radial is not there yet...')
        
        
    