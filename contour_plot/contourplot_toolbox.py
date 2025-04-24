# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 19:23:56 2025

@author: ychuang
"""


import matplotlib.pyplot as plt
from matplotlib import colors, cm
import matplotlib.tri as tri
from matplotlib.colors import LogNorm
from fit_data.fitting_method import fit_method_collection
import numpy as np
from numpy import ma
from matplotlib.offsetbox import AnchoredText




class contour_plot_method_collect:
        
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
       
    
    def contour_plot(self, plot_2dval, R_coord, Z_coord, quantity):
        CMAP = cm.viridis
        NORM= plt.Normalize(plot_2dval.min(), plot_2dval.max())
        
        plt.figure()
        plt.contourf(R_coord, Z_coord, plot_2dval, levels= 20, cmap=CMAP,norm=NORM)
        plt.title('{} contour plot'.format(quantity))
        
        
        SM= cm.ScalarMappable(NORM,CMAP)    
        plt.colorbar(SM)
        plt.show()
    
    
    
    def paper_contour(self, plot_2dval, R_coord, Z_coord, quantity, itername, 
                     log_bar, color_dic, A_dic, axs, cmap, norm, levels):
        
        
        vessel = self.data['vessel'][itername]
        
        
        if log_bar:
            if np.all(plot_2dval == 0):
                print('data_file is an zero matrix')
                
            elif np.any(plot_2dval == 0):
                plot_2dval = ma.masked_where(plot_2dval <= 0, plot_2dval)
                
                datamap = np.abs(plot_2dval)
                
                
                # CPB = cm.viridis
                # Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())
                axs.contourf(R_coord, Z_coord, datamap, levels= levels, 
                             cmap = cmap, norm = norm)
                
                
            else:
                
                datamap = np.abs(plot_2dval)
                
                # CPB = cm.viridis
                # Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())
                
                axs.contourf(R_coord, Z_coord, datamap,levels= levels, 
                             cmap = cmap, norm = norm)
                
        else:
            

            # NORM = axs.Normalize(plot_2dval.min(), plot_2dval.max())
            
            # CMAP = cm.viridis
            
            axs.contourf(R_coord, Z_coord, plot_2dval, levels= levels, cmap= cmap, norm = norm)
            
            
        axs.plot(vessel[:,0]/1000, vessel[:,1]/1000, color = 'black')
        # axs.plot(vessel[:,0]/1000, vessel[:,1]/1000, color = color_dic[itername])
    
    
    
    
    def eirene_contourplot_method(self, simu_direc, data, itername, axs, datname, norm_type):
        
        
        Nodes=np.fromfile('{}/fort.33'.format(simu_direc),sep=' ') #Alternatively use fort.33
        NN=int(Nodes[0])
        XNodes=Nodes[1:NN+1]
        YNodes=Nodes[NN+1:]
        
        
        CMAP = cm.viridis
        
        
        
        if norm_type == 'masklog':
            
            data_mask = ma.masked_where(data <= 0, data) 
            datamap = np.abs(data_mask)       
            NORM_data = colors.LogNorm(np.nanmin(datamap), np.nanmax(datamap))
        
        elif norm_type == 'allpositivelog':
                        
            data_mask = ma.masked_where(data == 0, data) 
            datamap = np.abs(data_mask)      
            NORM_data = colors.LogNorm(np.nanmin(datamap), np.nanmax(datamap))
        
        elif norm_type == 'std_normalize':
                        
            normalized_data = (data - np.mean(data)) / np.std(data)
            NORM_data = colors.Normalize(normalized_data.min(), normalized_data.max())

        
        elif norm_type == 'std_lognorm':
                        
            normalized_data = (data - np.mean(data)) / np.std(data)
            data_mask = ma.masked_where(normalized_data <= 0, normalized_data) 
            datamap = np.abs(data_mask)
            NORM_data = colors.LogNorm(np.nanmin(datamap), np.nanmax(datamap))
        
        
        elif norm_type == 'max_normalize':
            
            data_mask = ma.masked_where(data == 0, data)
            normalized_data = data_mask / np.max(data_mask)
            NORM_data = colors.Normalize(normalized_data.min(), normalized_data.max())
        
        
        elif norm_type == 'natural':
            
            NORM_data = colors.Normalize(data.min(), data.max())
            
                      
        else:
            
            NORM_data = colors.Normalize(data.min(), data.max())
            
        
        # Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())
    
        Triangles = np.loadtxt('{}/fort.34'.format(simu_direc), 
            skiprows=1, usecols=(1,2,3)) #Alternatively use fort.34
        # print(Triangles -1)
    
        TP = tri.Triangulation(XNodes, YNodes, triangles= (Triangles -1))
    
        # CMAP = cm.viridis
        
        
        """
        
        xlim = (0, 100)
        ylim = (0, 200)
        
        
        # Find which triangles have any point outside the limits
        mask = np.any(
            (XNodes[TP.triangles] < xlim[0]) | (XNodes[TP.triangles] > xlim[1]) |
            (YNodes[TP.triangles] < ylim[0]) | (YNodes[TP.triangles] > ylim[1]), axis=1)
        
        # Apply the mask
        TP.set_mask(mask)
        
        """
        
        
        
        
        
        ad = float(itername[0])
        ap = float(itername[1])
               
        axs.tripcolor(TP, data, shading = 'flat', cmap = CMAP, norm = NORM_data)
        axs.set_title('{0} $\Gamma_r ={1:.3f}$*$10^{{20}}$ 1/s, $q_r = {2:.3f}*10^5$ W'.format(datname, ad, ap))

        # plt.title('Neutral density contour plot outerleg')
    
        SM_data= cm.ScalarMappable(NORM_data, CMAP)
        
        plt.colorbar(SM_data)
        
        # plt.xlim(0, 100)              
        # plt.ylim(0, 200)
    
    
    
    def subplot_eirene_contour_method(self, simu_direc, data, itername, axs, cmap, norm):
        
        
        Nodes=np.fromfile('{}/fort.33'.format(simu_direc),sep=' ') #Alternatively use fort.33
        NN=int(Nodes[0])
        XNodes=Nodes[1:NN+1]
        YNodes=Nodes[NN+1:]
        
        
        data_mask = ma.masked_where(data <= 0, data)
        
        datamap = np.abs(data_mask)
        # CMAP = cm.viridis
        # NORM_data = colors.LogNorm(np.nanmin(datamap), np.nanmax(datamap))
        # Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())
    
        Triangles = np.loadtxt('{}/fort.34'.format(simu_direc), 
            skiprows=1, usecols=(1,2,3)) #Alternatively use fort.34
        # print(Triangles -1)
    
        TP = tri.Triangulation(XNodes, YNodes, triangles= (Triangles -1))
    
        # CMAP = cm.viridis
        
               
        IM = axs.tripcolor(TP, data, shading = 'flat', cmap = cmap, norm = norm)
        axs.set_title('Nd {}_{}'.format(itername[0], itername[1]))
        # plt.title('Neutral density contour plot outerleg')
        
        axs.set_xlim(0, 100)              
        axs.set_ylim(0, 200)
    
    
    
    def twcontourp(self, plot_2dval, R_coord, Z_coord, quantity, axs, norm_type, lv):
        CMAP = cm.viridis

        
        if norm_type == 'masklog':
            
            data_mask = ma.masked_where(plot_2dval <= 0, plot_2dval) 
            datamap = np.abs(data_mask)       
            NORM = colors.LogNorm(np.nanmin(datamap), np.nanmax(datamap))
        
        elif norm_type == 'allpositivelog':
                        
            data_mask = ma.masked_where(plot_2dval == 0, plot_2dval) 
            datamap = np.abs(data_mask)      
            NORM = colors.LogNorm(np.nanmin(datamap), np.nanmax(datamap))
        
        elif norm_type == 'std_normalize':
                        
            normalized_data = (plot_2dval - np.mean(plot_2dval)) / np.std(plot_2dval)
            NORM = colors.Normalize(normalized_data.min(), normalized_data.max())
        
        elif norm_type == 'std_lognorm':
                        
            normalized_data = (plot_2dval - np.mean(plot_2dval)) / np.std(plot_2dval)
            data_mask = ma.masked_where(normalized_data <= 0, normalized_data) 
            datamap = np.abs(data_mask)
            NORM = colors.LogNorm(np.nanmin(datamap), np.nanmax(datamap))
            
        elif norm_type == 'max_normalize':
            
            data_mask = ma.masked_where(plot_2dval == 0, plot_2dval)
            normalized_data = data_mask / np.max(data_mask)
            NORM = colors.Normalize(normalized_data.min(), normalized_data.max())
        
        elif norm_type == 'max_lognorm':
            
            data_mask = ma.masked_where(plot_2dval == 0, plot_2dval)
            normalized_data = data_mask / np.max(data_mask)
            datamap = np.abs(normalized_data)
            NORM = colors.LogNorm(np.nanmin(datamap), np.nanmax(datamap))
        
        
        elif norm_type == 'natural':
            
            NORM = colors.Normalize(plot_2dval.min(), plot_2dval.max())


        axs.contourf(R_coord, Z_coord, plot_2dval, levels= lv, cmap=CMAP,norm=NORM)
        axs.set_title('{} contour plot'.format(quantity))
        
        
        SM= cm.ScalarMappable(NORM,CMAP)    
        plt.colorbar(SM)
        plt.show()
    
    
    
    
    
    