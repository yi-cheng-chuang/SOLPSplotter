# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 18:24:05 2023

@author: user
"""

from SOLPSplotter_fit import profile_fit
import matplotlib.pyplot as plt
import SOLPS_set as ss
from matplotlib import colors, cm, ticker
import matplotlib.tri as tri
from matplotlib.colors import LogNorm
import fitting_method as fm
import numpy as np
from numpy import ma
from matplotlib.offsetbox import AnchoredText




class PlotContour(profile_fit):
    
       
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    

    
   
    
    def contour_plot(self, plot_2dval, R_coord, Z_coord, quantity, itername, 
                     log_bar, ma100, bounds, color_dic, A_dic):
        
        
        vessel = self.data['vessel'][itername]
        
        
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
                
    
    
    
    def improve_paper_plot(self, plotstyle):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            if plotstyle == 'single':
                
                for aa in self.data['dircomp']['multi_shift']:
                    
                    
                    neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
                    Radloc = np.transpose(self.data['grid']['RadLoc'][aa])[1:97, 1:37]
                    VertLoc = np.transpose(self.data['grid']['VertLoc'][aa])[1:97, 1:37]
                    
                    
                    self.contour_plot(plot_2dval = neuden, R_coord = Radloc , 
                    Z_coord = VertLoc, quantity = 'Neutral density', itername = aa, 
                             log_bar = True, ma100 = False, bounds = [], 
                             color_dic = color_dic, A_dic = A_dic)
            
            
            elif plotstyle == 'paper':
                
                
                comp_list = ['org', 'dot7']
                
                fig, axs = plt.subplots(1, 2, sharey= True)
                
                org_text = AnchoredText('{}'.format('(a) A = 1.4'), 
                                             loc='upper center')
                
                dot7_text = AnchoredText('{}'.format('(b) A = 2.8'), 
                                             loc='upper center')
                
                
                text_list = [org_text, dot7_text]
                
                

                for ii, aa in enumerate(comp_list):
                    
                    neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
                    RadLoc = np.transpose(self.data['grid']['RadLoc'][aa])[1:97, 1:37]
                    VertLoc = np.transpose(self.data['grid']['VertLoc'][aa])[1:97, 1:37]
                    
                    
                    if np.all(neuden == 0):
                        print('data_file is an zero matrix')
                        
                    elif np.any(neuden == 0):
                        plot_2dval = ma.masked_where(neuden <= 0, neuden)
                        
                        datamap = np.abs(plot_2dval)
                    
                    else:
                        
                        datamap = neuden
                    
                    
                    CPB = cm.viridis
                    Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())
                    
                    
                    self.paper_contour(plot_2dval = neuden, R_coord = RadLoc, Z_coord = VertLoc, 
                            quantity = 'Neutral density', itername = aa, 
                log_bar = True, color_dic = color_dic, A_dic = A_dic, axs = axs[ii], 
                cmap = CPB, norm = Lnorm, levels = 20)
                    
                    axs[ii].add_artist(text_list[ii])
                    axs[ii].set_xlabel('R [m]')
                    
                
                fig.supylabel('Z [m]')
                # fig.supxlabel('R [m]')
                fig.suptitle('Atomic neutral density $[m^{-3}]$ contour plot')
                smap = cm.ScalarMappable(Lnorm, CPB)    
                fig.colorbar(smap)
                plt.tight_layout()
    
    
    
    def rebuttal_NTplot(self, plotstyle):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                          'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            

            if plotstyle == 'paper':
                
                
                comp_list = ['org', 'dot3', 'dot5', 'dot7']
                
                fig, axs = plt.subplots(1, 4, sharey= True)
                
                org_text = AnchoredText('{}'.format('(a) A = 1.4'), 
                                              loc='upper center')
                
                dot3_text = AnchoredText('{}'.format('(b) A = 2.0'), 
                                              loc='upper center')
                
                dot5_text = AnchoredText('{}'.format('(c) A = 2.4'), 
                                              loc='upper center')
                
                dot7_text = AnchoredText('{}'.format('(d) A = 2.8'), 
                                              loc='upper center')
                
                
                text_list = [org_text, dot3_text, dot5_text, dot7_text]
                
                

                for ii, aa in enumerate(comp_list):
                    
                    
                    
                    
                    Te_J = self.data['b2fstate'][aa]['te'][1:97, 1:37]
                    
                    ev = 1.6021766339999999 * pow(10, -19)
                    neuden = Te_J / ev
                    
                    
                    RadLoc = np.transpose(self.data['grid']['RadLoc'][aa])[1:97, 1:37]
                    VertLoc = np.transpose(self.data['grid']['VertLoc'][aa])[1:97, 1:37]
                    
                    
                    if np.all(neuden == 0):
                        print('data_file is an zero matrix')
                        
                    elif np.any(neuden == 0):
                        plot_2dval = ma.masked_where(neuden <= 0, neuden)
                        
                        datamap = np.abs(plot_2dval)
                    
                    else:
                        
                        datamap = neuden
                    
                    
                    CPB = cm.viridis
                    Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())
                    
                    
                    self.paper_contour(plot_2dval = neuden, R_coord = RadLoc, Z_coord = VertLoc, 
                            quantity = 'Electron density', itername = aa, 
                log_bar = True, color_dic = color_dic, A_dic = A_dic, axs = axs[ii], 
                cmap = CPB, norm = Lnorm, levels = 20)
                    
                    axs[ii].add_artist(text_list[ii])
                    axs[ii].set_xlabel('R [m]')
                    
                
                axs[0].set_ylabel('Z [m]')
                fig.suptitle('Electron temperature $[eV]$ contour plot')
                smap = cm.ScalarMappable(Lnorm, CPB)    
                fig.colorbar(smap)
                plt.tight_layout(w_pad = 0.05)
    
    
    
    
    def rebuttal_NTchpplot(self, plotstyle):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                          'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            

            if plotstyle == 'paper':
                
                
                comp_list = ['org', 'dot3', 'dot5', 'dot7']
                
                fig, axs = plt.subplots(1, 4, sharey= True)
                
                org_text = AnchoredText('{}'.format('(a) A = 1.4'), 
                                              loc='upper center')
                
                dot3_text = AnchoredText('{}'.format('(b) A = 2.0'), 
                                              loc='upper center')
                
                dot5_text = AnchoredText('{}'.format('(c) A = 2.4'), 
                                              loc='upper center')
                
                dot7_text = AnchoredText('{}'.format('(d) A = 2.8'), 
                                              loc='upper center')
                
                
                text_list = [org_text, dot3_text, dot5_text, dot7_text]
                
                

                for ii, aa in enumerate(comp_list):
                    
                    neuden = self.data['b2fstate'][aa]['na'][1:97, 1:37, 1]
                    
                    
                    RadLoc = np.transpose(self.data['grid']['RadLoc'][aa])[1:97, 1:37]
                    VertLoc = np.transpose(self.data['grid']['VertLoc'][aa])[1:97, 1:37]
                    
                    
                    if np.all(neuden == 0):
                        print('data_file is an zero matrix')
                        
                    elif np.any(neuden == 0):
                        plot_2dval = ma.masked_where(neuden <= 0, neuden)
                        
                        datamap = np.abs(plot_2dval)
                    
                    else:
                        
                        datamap = neuden
                    
                    
                    CPB = cm.viridis
                    NORM= plt.Normalize(datamap.min(), datamap.max())
                    
                    
                    self.paper_contour(plot_2dval = neuden, R_coord = RadLoc, Z_coord = VertLoc, 
                            quantity = 'Electron density', itername = aa, 
                log_bar = True, color_dic = color_dic, A_dic = A_dic, axs = axs[ii], 
                cmap = CPB, norm = NORM, levels = 20)
                    
                    axs[ii].add_artist(text_list[ii])
                    axs[ii].set_xlabel('R [m]')
                    
                
                axs[0].set_ylabel('Z [m]')
                fig.suptitle('Electron temperature $[eV]$ contour plot')
                smap = cm.ScalarMappable(NORM, CPB)    
                fig.colorbar(smap)
                plt.tight_layout(w_pad = 0.05)
    
 
    
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
                
                # fig_dir  = ss.set_figdir()
                # plt.savefig('{}/{}_{}.png'.format(fig_dir, quant, aa), format='png')
        
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
            
            # fig_dir  = ss.set_figdir()
            # plt.savefig('{}/{}.png'.format(fig_dir, quant), format='png')
    
    
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
            
     
            
    def paper_contour_method(self, data, log_bar, itername, quant, color_dic,
                             A_dic, axs, cmap, norm, levels):
        
        if self.withshift == True and self.withseries == False:
            
            
            RadLoc = self.data['grid']['RadLoc'][itername]
            VertLoc = self.data['grid']['VertLoc'][itername]
            
            R_con = RadLoc[1:37, 1:97]
            Z_con = VertLoc[1:37, 1:97]
            
            

                
            self.paper_contour(plot_2dval = data, R_coord = R_con, norm = norm, 
                                 Z_coord = Z_con, quantity = quant, cmap= cmap,  
                                  itername = itername, log_bar = log_bar, 
                    color_dic = color_dic, A_dic= A_dic, axs = axs, levels = levels)
            
            
        
        
       
    
    


"""
Retire the load vessel function at contour file:
    
def load_vessel(self):
    if self.withshift == False and self.withseries == False:
        filedir = self.data['dirdata']['simutop']
        vessel_file = self.load_vessel_method(fdir = filedir)
        self.data['vessel'] = vessel_file
    
    elif self.withshift == True and self.withseries == False:
        vessel_file_dic = {}
        for aa in self.data['dircomp']['multi_shift']:
            filedir = self.data['dirdata']['simutop'][aa]
            vessel_file = self.load_vessel_method(fdir = filedir)
            vessel_file_dic[aa] = vessel_file
        
        self.data['vessel'] = vessel_file_dic
    
    elif self.withshift == False and self.withseries == True:
        # series_rep = list(self.data['dircomp']['Attempt'].keys())[0]
        filedir = self.data['dirdata']['simutop']
        vessel_file = self.load_vessel_method(fdir = filedir)
        self.data['vessel'] = vessel_file
    
    elif self.withshift == True and self.withseries == True:
        print('load_vessel function is not there yet!')
    
    else:
        print('load_vessel function has a bug')


def load_vessel_method(self, fdir):
    # try:
    #     WallFile = np.loadtxt('{}/mesh.extra'.format(self.data['dirdata']['tbase']))
    # except:
    #     print('mesh.extra file not found! Using vvfile.ogr instead')
    #     WallFile=None
    
    try:
        VVFILE = np.loadtxt('{}/baserun/vvfile.ogr'.format(fdir))
    except:
        print('load_vessel_method has a bug!')

    return VVFILE

"""


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
        
"""       