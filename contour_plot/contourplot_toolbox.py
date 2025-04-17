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