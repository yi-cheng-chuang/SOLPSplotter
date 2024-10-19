# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 10:02:30 2024

@author: ychuang
"""



from R_diff_calc import Diff_R_calc
from SOLPSplotter_contour import PlotContour
from matplotlib.offsetbox import AnchoredText
import load_B2_data_method as lBdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker
import matplotlib.tri as tri
from matplotlib.colors import LogNorm
from numpy import ma


class target_contour(Diff_R_calc, PlotContour):
    def __init__(self, DefaultSettings, loadDS):
        Diff_R_calc.__init__(self, DefaultSettings, loadDS)
        PlotContour.__init__(self, DefaultSettings, loadDS)
    
        
    
    
    def iout_paper_plot(self, plotstyle, dataname, sideswitch):
        
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
                
                
                comp_list = ['org','dot3', 'dot5', 'dot7']
                
                fig, axs = plt.subplots(2, 2, sharey= True)
                
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
                    
                    if dataname == 'neutral density':
                        neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
                        input_dat = neuden
                    
                    elif dataname == 'Poloidal flux':
                        
                        fnaxs = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                        vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                        hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                                               
                        fnnx = np.divide(fnaxs, vol)
                        fnax = np.multiply(fnnx, hx)
                        input_dat = fnax
                    
                    elif dataname == 'Source':
                        
                        source = self.data['b2wdat'][aa]['b2npc_sna'][0][1:97, 1:37]
                        vol = self.data['b2wdat'][aa]['vol'][1:97, 1:37]
                        
                        sterm = np.divide(source, vol)
                        input_dat = sterm
                    
                    elif dataname == 'hx':
                        
                        hx = self.data['b2wdat'][aa]['hx'][1:97, 1:37]
                                               
                        input_dat = hx
                    
                    
                    
                    RadLoc = np.transpose(self.data['grid']['RadLoc'][aa])[1:97, 1:37]
                    VertLoc = np.transpose(self.data['grid']['VertLoc'][aa])[1:97, 1:37]
                    
                    
                    if np.all(input_dat == 0):
                        print('data_file is an zero matrix')
                        
                    elif np.any(input_dat <= 0):
                        
                        if dataname == 'Poloidal flux':
                            
                            if sideswitch == 'negative':
                                
                                plot_2dval = ma.masked_where(input_dat >= 0, input_dat)
                                
                                datamap = np.abs(plot_2dval)
                            
                            elif sideswitch == 'positive':
                                
                                plot_2dval = ma.masked_where(input_dat <= 0, input_dat)
                                
                                datamap = np.abs(plot_2dval)
                            
                            elif sideswitch == 'both':
                                
                                plot_2dval = np.abs(input_dat)
                                
                                datamap = ma.masked_where(plot_2dval <= 0, plot_2dval)
                            
                            else:
                                print('check sideswitch for poloidal flux')
                                
                                
                        elif dataname == 'Source':
                            
                            if sideswitch == 'positive':
                                
                                plot_2dval = ma.masked_where(input_dat <= 0, input_dat)
                                
                                datamap = np.abs(plot_2dval)
                            
                            elif sideswitch == 'negative':
                                
                                plot_2dval = ma.masked_where(input_dat >= 0, input_dat)
                                
                                datamap = np.abs(plot_2dval)
                                
                            elif sideswitch == 'both':
                                
                                plot_2dval = np.abs(input_dat)
                                
                                datamap = plot_2dval
                            
                            else:
                                print('check sideswitch for source')
                        
                                
                        elif dataname == 'hx':
                            
                            if sideswitch == 'negative':
                                
                                plot_2dval = ma.masked_where(input_dat >= 0, input_dat)
                                
                                datamap = np.abs(plot_2dval)
                            
                            elif sideswitch == 'positive':
                                
                                plot_2dval = ma.masked_where(input_dat <= 0, input_dat)
                                
                                datamap = np.abs(plot_2dval)
                            
                            elif sideswitch == 'both':
                                
                                plot_2dval = np.abs(input_dat)
                                
                                datamap = ma.masked_where(plot_2dval <= 0, plot_2dval)
                            
                            else:
                                print('check sideswitch for hx')
                        
                        
                        
                        else:
                            plot_2dval = ma.masked_where(input_dat <= 0, input_dat)
                            
                            datamap = np.abs(plot_2dval)
                            
                        
                    else:
                        
                        datamap = input_dat
                    
                    print(datamap.max())
                    print(datamap.min())
                    
                    
                    
                    
                    CPB = cm.viridis
                    Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())
                    
                    if ii < 2:
                        
                        self.paper_contour(plot_2dval = datamap, R_coord = RadLoc, Z_coord = VertLoc, 
                                quantity = dataname, itername = aa, 
                    log_bar = True, color_dic = color_dic, A_dic = A_dic, axs = axs[ii, 0], 
                    cmap = CPB, norm = Lnorm, levels = 20)
                    
                    else:
                        
                        ik = ii % 2
                        self.paper_contour(plot_2dval = datamap, R_coord = RadLoc, Z_coord = VertLoc, 
                                quantity = dataname, itername = aa, 
                    log_bar = True, color_dic = color_dic, A_dic = A_dic, axs = axs[ik, 1], 
                    cmap = CPB, norm = Lnorm, levels = 20)
                        
                    
                    
                    
                    
                    if aa == 'org':
                            
                        axs[ii, 0].set_xlim(0, 2)
                        
                    elif aa == 'dot3':
                        
                        axs[ii, 0].set_xlim(0.3, 2.3)
                    
                    elif aa == 'dot5':
                        
                        ik = ii % 2
                        axs[ik, 1].set_xlim(0.5, 2.5)
                    
                    elif aa == 'dot7':
                        
                        ik = ii % 2
                        axs[ik, 1].set_xlim(0.7, 2.7)
                    
                    if ii < 2:
                        
                        axs[ii, 0].add_artist(text_list[ii])
                        axs[ii, 0].set_xlabel('R [m]')              
                        axs[ii, 0].set_ylim(-2, -0.7)
                    
                    else:
                        
                        ik = ii % 2
                        
                        axs[ik, 1].add_artist(text_list[ii])
                        axs[ik, 1].set_xlabel('R [m]')              
                        axs[ik, 1].set_ylim(-2, -0.7)
                    
                    
                    
                fig.supylabel('Z [m]')
                fig.subplots_adjust(right=0.8)
                cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
                # fig.supxlabel('R [m]')
                if dataname == 'neutral density':
                    fig.suptitle('Atomic{} contour plot'.format(dataname))
                
                else:
                    fig.suptitle('{} {} part contour plot'.format(dataname, sideswitch))
                    
                    
                smap = cm.ScalarMappable(Lnorm, CPB)    
                fig.colorbar(smap, cax= cbar_ax)
                # plt.tight_layout()
        
        
        
        



    
    
    