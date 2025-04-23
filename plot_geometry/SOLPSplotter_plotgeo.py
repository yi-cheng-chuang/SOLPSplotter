# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 18:31:57 2025

@author: ychuang
"""

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import colors, cm, ticker


class plot_geo:
    
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
        
    
        
    
    def contour_plot(self, plot_2dval, R_coord, Z_coord, quantity):
        CMAP = cm.viridis
        NORM= plt.Normalize(plot_2dval.min(), plot_2dval.max())
        
        fig, ax = plt.subplots()
        ax.contourf(R_coord, Z_coord, plot_2dval, levels= 20, cmap=CMAP,norm=NORM)
        ax.set_title('{} contour plot'.format(quantity))
        
        
        SM= cm.ScalarMappable(NORM,CMAP)    
        fig.colorbar(SM, ax=ax)
        plt.show()
    
    def contour_log_plot(self, plot_2dval, R_coord, Z_coord, quantity):
        CMAP = cm.viridis
        NORM= colors.LogNorm(plot_2dval.min(), plot_2dval.max())
        
        fig, ax = plt.subplots()
        ax.contourf(R_coord, Z_coord, plot_2dval, levels= 20, cmap=CMAP,norm=NORM)
        ax.set_title('{} contour plot'.format(quantity))
        
        
        SM= cm.ScalarMappable(NORM,CMAP)    
        fig.colorbar(SM, ax=ax)
        plt.show()
    
    def plot_g(self):
        
        psi = self.data['gfile']['gcomp']['psiN']
        rarray = self.data['gfile']['gcomp']['gR']
        zarray = self.data['gfile']['gcomp']['gZ']
        psirz = self.data['gfile']['g']['psirz']
        psi_sym = self.data['gfile']['gcomp']['sym_psi']
        psiN_sym = self.data['gfile']['gcomp']['psiN_sym']
        
        self.contour_log_plot(plot_2dval = psiN_sym, R_coord = rarray, 
                         Z_coord = zarray, quantity = 'PsiN_sym')
        
        self.contour_plot(plot_2dval = psi_sym, R_coord = rarray, 
                         Z_coord = zarray, quantity = 'Psi_sym')
    
    
    def plot_sec(self):
               
        
        RadLoc = self.data['grid']['RadLoc']
        VertLoc = self.data['grid']['VertLoc']

        
        pol_list = [22, 23, 24, 25]
        
        # Generate 10 evenly spaced integers between 1 and 100
        numbers = np.linspace(61, 83, num=23, dtype=int)
        
        # Convert to a list
        integer_list = numbers.tolist()
        
        # print(integer_list)
        

        plt.figure()
        for in_pol in integer_list:
            crloc = RadLoc[:, int(in_pol)]
            czloc = VertLoc[:, int(in_pol)]
            plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
        plt.xlabel('R: [m]')
        plt.ylabel('Z: [m]')
    
    
    
    
    
    
    
    
    
    
    
    