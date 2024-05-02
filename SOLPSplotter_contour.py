# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 18:24:05 2023

@author: user
"""

from SOLPSplotter_fit import profile_fit
import matplotlib.pyplot as plt
import SOLPS_set as ss
import Contourplot_method as cpm
from matplotlib import colors, cm, ticker
import matplotlib.tri as tri
from matplotlib.colors import LogNorm
import fitting_method as fm
import numpy as np
from numpy import ma
from matplotlib.offsetbox import AnchoredText




class PlotContour(profile_fit):
    def __init__(self, DefaultSettings, loadDS):
        profile_fit.__init__(self, DefaultSettings, loadDS)
        
        self.Publish = DefaultSettings['Publish']
        self.data['DefaultSettings']['Publish'] = self.Publish
    
    
    def set_plot(self):
        if self.Publish == 'b2plottersetting':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 4, markersize= 7)
            plt.rcParams.update({'font.size': 12})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
  
        else:
            print('Publish setting is incorrect or add another setting')
    
    
    
    def calc_flux_expansion_line_method(self, RR_sep, arcR):
                   
        flux_fit_dic = fm.flux_expand_fit(RRsep = RR_sep, arclength = arcR)
        
        flux_expand = flux_fit_dic['flux_fitcoe'][0]
        a_flux_exp = flux_expand*np.ones(self.data['b2fgeo']['ny'])
        
        return a_flux_exp
    
    
    def contour_plot(self, plot_2dval, R_coord, Z_coord, quantity, itername, 
                     log_bar, ma100, bounds, color_dic, A_dic):
        
        
        vessel = self.data['vessel'][itername]
        
        fig, axs = plt.subplots()
        if log_bar:
            if np.all(plot_2dval == 0):
                print('data_file is an zero matrix')
            elif np.any(plot_2dval == 0):
                plot_2dval = ma.masked_where(plot_2dval <= 0, plot_2dval)
                
                datamap = np.abs(plot_2dval)
                
                
                CPB = cm.viridis
                Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())
                plt.contourf(R_coord, Z_coord, datamap, levels= 20, 
                             cmap = CPB, norm = Lnorm)
                
                smap = cm.ScalarMappable(Lnorm, CPB)    
                plt.colorbar(smap)
            else:
                
                datamap = np.abs(plot_2dval)
                
                if len(bounds) == 0:
                    print('the mask upper and lower bound is not set!')
                
                else:
                    if ma100:
                        datamap = ma.masked_where(datamap >= bounds['max'], datamap)
                        datamap = ma.masked_where(datamap <= bounds['min'], datamap)
                    else:
                        pass
                    
                self.data['mask'] = datamap
                CPB = cm.viridis
                Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())
                plt.contourf(R_coord, Z_coord, datamap,levels= 20, 
                             cmap = CPB, norm = Lnorm)
                
                smap = cm.ScalarMappable(Lnorm, CPB)    
                plt.colorbar(smap)
                
            
        else:
            
            # datamap = np.abs(plot_2dval)
            
            if ma100 == True and bounds != None:
                plot_2dval = ma.masked_where(plot_2dval >= 1500, plot_2dval)
                plot_2dval = ma.masked_where(plot_2dval <= -1500, plot_2dval)
                NORM = plt.Normalize(vmin = bounds['min'], vmax = bounds['max'])
            else:
                NORM = plt.Normalize(plot_2dval.min(), plot_2dval.max())
            
            
            # CMAP = 'RdBu'
            CMAP = cm.viridis
            
            plt.contourf(R_coord, Z_coord, plot_2dval, levels= 40, cmap= CMAP, norm = NORM)
            
            SM= cm.ScalarMappable(NORM,CMAP)    
            plt.colorbar(SM)
        
        
        # cs = ax.contourf(X, Y, z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
        # cbar = fig.colorbar(cs) 
        
        
        
        if itername == None:
            plt.title('{} contour plot'.format(quantity))
        
        else:
            plt.title('{} for aspect ratio {}'.format(quantity, A_dic[itername]), )
            
        plt.plot(vessel[:,0]/1000, vessel[:,1]/1000, color = color_dic[itername])
                
        
        
    
    def flux_expansion_contour_plot_method(self, RR_sep, flux_expand_map, itername):
        
        for pol_loc in range(self.data['b2fgeo']['nx']):
            
            if itername == None:
                arcR = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]

            else:
                arcR = self.data['dsa']['dsa_{}'.format(pol_loc)][itername]['dsa_{}_val'.format(pol_loc)]
            
            
            flux_fit_dic = fm.flux_expand_fit(RRsep = RR_sep, arclength = arcR)
            
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
        
        
        cpm.contour_plot(plot_2dval = flux_expand_map, R_coord = R_con, 
                         Z_coord = Z_con, quantity = 'flux expansion')
        
        
        # self.plot_vessel(itername = itername, independent = False, meter = True)
        # plt.show()
        
        return contour_dic
        
       
    def flux_expansion_contour_plot(self):
        
        if self.withshift == False and self.withseries == False:
            
            for pol_loc in range(self.data['b2fgeo']['nx']):
                self.calc_dsa(pol_loc)
            
                        
            RR_sep = self.data['midplane_calc']['R_Rsep']
            flux_expand_map = np.zeros([self.data['b2fgeo']['ny'], self.data['b2fgeo']['nx']])
            
            contour_dic = self.flux_expansion_contour_plot_method(RR_sep = RR_sep, 
                            flux_expand_map = flux_expand_map, itername = None)
            
            self.data['flux_contour'] = contour_dic
            
            
        
        elif self.withshift == True and self.withseries == False:
            
            contour_dic = {}
            
            for pol_loc in range(self.data['b2fgeo']['nx']):
                self.calc_dsa(pol_loc)
            
            for aa in self.data['dircomp']['multi_shift']:
            
                RR_sep = self.data['midplane_calc'][aa]['R_Rsep']
                flux_expand_map = np.zeros([self.data['b2fgeo']['ny'], self.data['b2fgeo']['nx']])
            
            
                contour_dic[aa] = self.flux_expansion_contour_plot_method(RR_sep = RR_sep, 
                                flux_expand_map = flux_expand_map, itername = aa)
                
            
            self.data['flux_contour'] = contour_dic
            
        
        elif self.withshift == False and self.withseries == True:
            
            for pol_loc in range(self.data['b2fgeo']['nx']):
                self.calc_dsa(pol_loc)
                                    
            RR_sep = self.data['midplane_calc']['R_Rsep']
            flux_expand_map = np.zeros([self.data['b2fgeo']['ny'], self.data['b2fgeo']['nx']])
            
            contour_dic = self.flux_expansion_contour_plot_method(RR_sep = RR_sep, 
                            flux_expand_map = flux_expand_map, itername = None)
            
            self.data['flux_contour'] = contour_dic
        
        elif self.withshift == True and self.withseries == True:
            print('calc_flux_expansion is not there yet, to be continue...')
            
        else:
            print('There is a bug')


        
    def load_vessel(self):
        if self.withshift == False and self.withseries == False:
            filedir = self.data['dirdata']['simutop']
            vessel_file = cpm.load_vessel_method(fdir = filedir)
            self.data['vessel'] = vessel_file
        
        elif self.withshift == True and self.withseries == False:
            vessel_file_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                filedir = self.data['dirdata']['simutop'][aa]
                vessel_file = cpm.load_vessel_method(fdir = filedir)
                vessel_file_dic[aa] = vessel_file
            
            self.data['vessel'] = vessel_file_dic
        
        elif self.withshift == False and self.withseries == True:
            # series_rep = list(self.data['dircomp']['Attempt'].keys())[0]
            filedir = self.data['dirdata']['simutop']
            vessel_file = cpm.load_vessel_method(fdir = filedir)
            self.data['vessel'] = vessel_file
        
        elif self.withshift == True and self.withseries == True:
            print('load_vessel function is not there yet!')
        
        else:
            print('load_vessel function has a bug')
            
    
    
    def plot_vessel_method(self, vessel_data, shift_value, independent, meter, color_dic, itername):
        
        if independent:
            plt.figure(figsize=(7,7))
        else:
            pass
        
        if meter:
            plt.plot(vessel_data[:,0]/1000, vessel_data[:,1]/1000, color = color_dic[itername])
            # plt.xlabel('R')
                    
            tick_label = np.arange(0 + shift_value, 2.1 + shift_value, 0.5)
            
            ytick = np.arange(-1.1, 1.1, 0.5)
            plt.yticks(ytick)
            
            plt.xticks(tick_label)
            
            
        else:
            plt.plot(vessel_data[:,0], vessel_data[:,1], color = color_dic[itername])
            # plt.xlabel('R')
                    
            tick_label = np.arange(0 + shift_value, 2100 + shift_value, 500)
            
            ytick = np.arange(-1100, 1100, 500)
            plt.yticks(ytick)
            
            plt.xticks(tick_label)
            
        
        if independent:
            plt.title('vessel')
            plt.show()
        
        else:
            pass
    
            
    def paper_vessel_method(self, vessel_data, shift_value, meter, 
                           color_dic, A_dic, itername, axs):
        
        if meter:
            axs.plot(vessel_data[:,0]/1000, vessel_data[:,1]/1000, 
        color = color_dic[itername], label= 'aspect ratio = {}'.format(A_dic[itername]))
            
            
        else:
            axs.plot(vessel_data[:,0], vessel_data[:,1], 
    color = color_dic[itername], label= 'aspect ratio = {}'.format(A_dic[itername]))
            
    
            
    def plot_vessel(self, itername, independent, meter):
        
        if self.withshift == False and self.withseries == False:
            
            vessel = self.data['vessel']
            shift = self.data['dircomp']['shift_value']*1000
            
            self.plot_vessel_method(vessel_data = vessel, shift_value = shift, 
                                    independent = independent, meter= meter)
        
        elif self.withshift == True and self.withseries == False:
            
            vessel = self.data['vessel'][itername]
            shift = self.data['dircomp']['shift_dic'][itername]*1000
            
            self.plot_vessel_method(vessel_data = vessel, shift_value = shift, 
                                    independent = independent, meter = meter)
        
        elif self.withshift == False and self.withseries == True:
            
            vessel = self.data['vessel']
            shift = self.data['dircomp']['shift_value']
            
            self.plot_vessel_method(vessel_data = vessel, shift_value = shift, 
                                    independent = independent, meter = meter)
        
        else:
            
            print('plot_vessel function is not there yet!')
    
    
    
    def shaded_area(self):
        
            
        # CMAP = 'RdBu'
        CMAP = cm.viridis
        
        R_coord = self.data['grid']['RadLoc']['org']
        Z_coord = self.data['grid']['VertLoc']['org']
    
        shade = np.ones([38, 98])
        
        plt.contourf(R_coord, Z_coord, shade, levels= 1, cmap= 'Blues')
        
        
        # cs = ax.contourf(X, Y, z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
        # cbar = fig.colorbar(cs)
    
    
    
    
    
    def shift_vessel_in_one(self):
        
        if self.withshift == True and self.withseries == False:
            
            fig, axs = plt.subplots()
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            anchored_text = AnchoredText('{}'.format('vessel cross section'), loc='upper right')
            
            self.shaded_area()
            
            # ylabel_text = AnchoredText('{}'.format('Z [m]'), loc='center left')
            
            for aa in self.data['dircomp']['multi_shift']:
                
                vessel = self.data['vessel'][aa]
                shift = self.data['dircomp']['shift_dic'][aa]*1000
                
                self.paper_vessel_method(vessel_data = vessel, shift_value = shift,
            meter = True, color_dic = color_dic, itername = aa, axs = axs, A_dic = A_dic)
                
                axs.add_artist(anchored_text)
                # axs.add_artist(ylabel_text)
                axs.set_xlabel('R [m]')
                axs.set_ylabel('Z [m]')
                axs.legend(loc= 'center right', fontsize=10)
            
            
            axs.set_aspect('equal')
            
            fig.savefig('vessel.pdf')
        
        
        
    def iout_contour_plot(self, quant, log_bar, ma100, bounds):
        
        if self.withshift == False and self.withseries == False:
            
            data = self.data['iout_data'][quant]
            
            RadLoc = self.data['grid']['RadLoc']
            VertLoc = self.data['grid']['VertLoc']
            
            R_con = RadLoc[1:37, 1:97]
            Z_con = VertLoc[1:37, 1:97]
                        
            self.contour_plot(plot_2dval = data, R_coord = R_con, 
                             Z_coord = Z_con, quantity = quant, ma100 = ma100, 
                             itername = None, log_bar = log_bar, bounds = bounds)
            
            fig_dir  = ss.set_figdir()
            plt.savefig('{}/{}.png'.format(fig_dir, quant), format='png')
        
        elif self.withshift == True and self.withseries == False:
            
            data_dic = self.data['iout_data'][quant]
            
            RadLoc = self.data['grid']['RadLoc']
            VertLoc = self.data['grid']['VertLoc']
            
            for aa in self.data['dircomp']['multi_shift']:
                
                data = self.data['iout_data'][quant][aa]
                
                R_con = RadLoc[aa][1:37, 1:97]
                Z_con = VertLoc[aa][1:37, 1:97]
                            
                self.contour_plot(plot_2dval = data, R_coord = R_con, 
                                 Z_coord = Z_con, quantity = quant, ma100 = ma100,
                                 itername = aa, log_bar = log_bar, bounds = bounds)
                
                fig_dir  = ss.set_figdir()
                plt.savefig('{}/{}_{}.png'.format(fig_dir, quant, aa), format='png')
        
        elif self.withshift == False and self.withseries == True:
            
            data_dic = self.data['iout_data'][quant]
            
            RadLoc = self.data['grid']['RadLoc']
            VertLoc = self.data['grid']['VertLoc']
            
            
            for aa in list(self.data['dircomp']['Attempt'].keys()):
                
                data = self.data['iout_data'][quant][aa]
            
                R_con = RadLoc[1:37, 1:97]
                Z_con = VertLoc[1:37, 1:97]
                            
                self.contour_plot(plot_2dval = data, R_coord = R_con, 
                                 Z_coord = Z_con, quantity = quant, ma100 = ma100, 
                                 itername = aa, log_bar = log_bar, bounds = bounds)
        
        else:
            print('iout_contour_plot function is not there yet!')

    
    
    def plot_change(self, quant, log_bar, itername, ma100):
        
        if self.withshift == True and self.withseries == False:
            
            data = self.data['iout_data'][quant][itername]
            
            RadLoc = self.data['grid']['RadLoc'][itername]
            VertLoc = self.data['grid']['VertLoc'][itername]
            
            R_con = RadLoc[1:37, 1:97]
            Z_con = VertLoc[1:37, 1:97]
                        
            self.contour_plot(plot_2dval = data, R_coord = R_con, 
                             Z_coord = Z_con, quantity = quant, 
                             itername = itername, log_bar = log_bar, ma100= ma100)
            
            fig_dir  = ss.set_figdir()
            plt.savefig('{}/{}.png'.format(fig_dir, quant), format='png')
    
    
    def plot_change_data(self, data, log_bar, itername, quant, ma100, 
                         bounds, color_dic, A_dic):
        
        if self.withshift == True and self.withseries == False:
            
            
            RadLoc = self.data['grid']['RadLoc'][itername]
            VertLoc = self.data['grid']['VertLoc'][itername]
            
            R_con = RadLoc[1:37, 1:97]
            Z_con = VertLoc[1:37, 1:97]
                        
            self.contour_plot(plot_2dval = data, R_coord = R_con, 
                             Z_coord = Z_con, quantity = quant, bounds = bounds, 
                                     itername = itername, log_bar = log_bar, 
                ma100= ma100, color_dic = color_dic, A_dic= A_dic)
            
            # fig_dir  = ss.set_figdir()
            # plt.savefig('{}/{}.png'.format(fig_dir, quant), format='png')
            
            
        
        
    
    


"""    
Need to fix!!
  
    def plot_all_radial(self):
            
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
        
        
"""       
        
        
        
"""
backup:
    
def load_vessel_method(self, itername):
    # try:
    #     WallFile = np.loadtxt('{}/mesh.extra'.format(self.data['dirdata']['tbase']))
    # except:
    #     print('mesh.extra file not found! Using vvfile.ogr instead')
    #     WallFile=None
    
    if self.withshift == False and self.withseries == False:
        VVFILE = np.loadtxt('{}/baserun/vvfile.ogr'.format(self.data['dirdata']['simutop']))
    
    elif self.withshift == True and self.withseries == False:
        VVFILE = np.loadtxt('{}/baserun/vvfile.ogr'.format(self.data['dirdata']['simutop'][itername]))
    
    elif self.withshift == False and self.withseries == True:
        VVFILE = np.loadtxt('{}/baserun/vvfile.ogr'.format(self.data['dirdata']['simutop']))

        
    elif self.withshift == True and self.withseries == True:
        print('load_vessel_method function is not there yet!')
    
    else:
        print('load_vessel_method function has a bug')
    
    # if plot:
    #     plt.plot
    return VVFILE



def contour_plot(self, plot_2dval, R_coord, Z_coord, quantity):
    CMAP = cm.viridis
    NORM= plt.Normalize(plot_2dval.min(), plot_2dval.max())
    
    plt.figure(figsize=(6,12))
    plt.contourf(R_coord, Z_coord, plot_2dval, levels= 20, cmap=CMAP,norm=NORM)
    plt.title('{} contour plot'.format(quantity))
    
    
    SM= cm.ScalarMappable(NORM,CMAP)    
    plt.colorbar(SM)


for pol_loc in range(self.data['b2fgeo']['nx']):
    self.calc_dsa(pol_loc)
    arcR = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]
    
    a_flux_exp = self.calc_flux_expansion_line_method(RR_sep = RR_sep, 
                                                      arcR = arcR)
    
    flux_expand_map[:, pol_loc] = a_flux_exp


RadLoc = self.data['grid']['RadLoc']
VertLoc = self.data['grid']['VertLoc']

R_con = RadLoc[1:37, 1:97]
Z_con = VertLoc[1:37, 1:97]

contour_dic = {'R_coord': R_con, 'Z_coord': Z_con, 
               'flux_map': flux_expand_map}
self.data['flux_contour'] = contour_dic


# map_flat = flux_expand_map.flatten()


cpm.contour_plot(plot_2dval = flux_expand_map, R_coord = RadLoc, 
                 Z_coord = VertLoc, quantity = 'flux expansion')



"""