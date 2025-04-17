# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 18:54:47 2025

@author: ychuang
"""

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from numpy import ma
from matplotlib.offsetbox import AnchoredText
from matplotlib import colors, cm
from twscan_module.twinscan_prepare import twscan_assist




class Eirene_contour_save:


    
    def __init__(self, DF, data, twa: twscan_assist):
        
        self.DF = DF
        self.data = data
        self.twa = twa


    def eirene_contourplot_tool(self, simudir, data, plot33):
        
        
        data_mask = ma.masked_where(data <= 0, data)
        
        datamap = np.abs(data_mask)
        CMAP = cm.viridis
        NORM_data = colors.LogNorm(np.nanmin(datamap), np.nanmax(datamap))
        # Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())
        
        
        if plot33:
            
            Nodes=np.fromfile('{}/fort.33'.format(simudir),sep=' ') #Alternatively use fort.33
            NN=int(Nodes[0])
            XNodes=Nodes[1:NN+1]
            YNodes=Nodes[NN+1:]
        
        
            numberlist = np.zeros(NN)
            for i in range(NN):
                numberlist[i] = i
        
            plt.figure()
            plt.scatter(XNodes[:500], YNodes[:500])
            
         
        Triangles = np.loadtxt('{}/fort.34'.format(simudir), 
            skiprows=1, usecols=(1,2,3)) #Alternatively use fort.34
        # print(Triangles -1)
    
        TP = tri.Triangulation(XNodes, YNodes, triangles= (Triangles -1))
    
        # CMAP = cm.viridis
    
        plt.figure()
        plt.tripcolor(TP, data, shading = 'flat', cmap = CMAP, norm = NORM_data)
        plt.title('Neutral density contour plot')
        # plt.title('Neutral density contour plot outerleg')
    
        SM_data= cm.ScalarMappable(NORM_data, CMAP)    
        plt.colorbar(SM_data)
    
    
    
    def eirene_contour_plot(self):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == False:
            
            simu_dir = self.data['dirdata']['simudir']
            dat = self.data['ft46']['pdena'][:, 0]
            abs_dat = np.absolute(dat)
            
            self.eirene_contourplot_method(simudir = simu_dir, data = abs_dat)
            
        
        
        elif withshift == False and withseries == True:
            
            print('we are working on it!')
            
            
            
        else:
            print('eirene_contour_plot function needs to be checked!')
    
    
    
    
    