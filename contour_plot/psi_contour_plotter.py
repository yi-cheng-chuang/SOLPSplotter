# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 17:58:04 2025

@author: ychuang
"""



import matplotlib.pyplot as plt
from matplotlib import colors, cm, ticker
import matplotlib.tri as tri
from matplotlib.colors import LogNorm
import numpy as np
from numpy import ma
from matplotlib.offsetbox import AnchoredText




class psiN_plotContour:
    
       
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    def psiN_populator(self):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        
        if withshift == False and withseries == False:
            
            
            psi2dmap = self.data['gfile']['gcomp']['interp_dic']['RBS']
            VVFILE = self.data['vessel']
            
            r = np.linspace(0.0, 2.0, 200)
            z = np.linspace(-2.0, 2.0, 200)
            R, Z = np.meshgrid(r, z, indexing='ij')
            V = psi2dmap(r, z, grid=True)
            
            fig, ax = plt.subplots()
            
            # Base contours
            cs = ax.contour(R, Z, V, levels=30, linewidths=2.0)
            
            # ψ_N = 1 contour
            cs1 = ax.contour(R, Z, V, levels=[1.0], colors='red', linewidths=2.8)
            
            # Put the label exactly where you want it (pick a clear spot in your plot)
            label_texts = ax.clabel(
                cs1,
                levels=[1.0],
                fmt={1.0: r'$\psi_N = 1$'},
                manual=[(1.55, -0.80)],   # <-- change (R, Z) to your preferred location
                inline=True,
                inline_spacing=5,
                fontsize=10,
            )
            
            # Optional: add a white halo so the label is readable over lines
            import matplotlib.patheffects as pe
            for txt in label_texts:
                txt.set_path_effects([pe.withStroke(linewidth=3, foreground="white")])
            
            ax.plot(VVFILE[:,0]/1000, VVFILE[:,1]/1000, color='black')
            ax.set_xlabel('R'); ax.set_ylabel('Z')
            ax.set_aspect('equal', adjustable='box')
            fig.colorbar(cs, ax=ax, label='$\psi_N$')
            plt.show()

    