# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 16:38:30 2024

@author: ychuang
"""

from SOLPSplotter_plot import Opacity_study
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import matplotlib.tri as tri

import fitting_method as fm
import numpy as np





def contour_plot(plot_2dval, R_coord, Z_coord, quantity):
    CMAP = cm.viridis
    NORM= plt.Normalize(plot_2dval.min(), plot_2dval.max())
    
    plt.figure(figsize=(6,12))
    plt.contourf(R_coord, Z_coord, plot_2dval, levels= 20, cmap=CMAP,norm=NORM)
    plt.title('{} contour plot'.format(quantity))
    
    
    SM= cm.ScalarMappable(NORM,CMAP)    
    plt.colorbar(SM)
    plt.show()

def load_vessel_method(fdir):
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