# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 09:45:57 2023

@author: user
"""

from B2plotter_class import B2plotter
import xarray as xr
import matplotlib.pyplot as plt
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate
from scipy.optimize import curve_fit
from matplotlib.path import Path
from scipy.stats import binned_statistic
import numpy as np


class psi_RRsep_calc(B2plotter):
    def __init__(self, DEV, withshift, withseries, DefaultSettings):
        B2plotter.__init__(self, DEV, withshift, withseries, DefaultSettings)
    
    
    def plot_psi(self, pol_list):
        self.data['poloidal_index'] = pol_list
        
        for j in pol_list:
            self.calcpsi_1D(pol_loc= j)
        
             
        ln = len(pol_list)
        pol_loc = np.zeros(ln)
        i = 0
        
        for ii in self.data['poloidal_index']:
            pol_loc[i] = int(ii)
            i = i + 1
        
        rad_range = self.data['DefaultSettings']['YDIM']
        
        
        if self.withshift == False and self.withseries == False:
            
            plt.figure(figsize=(7,7))
            for xc in pol_list:
                psi_plot = self.data['psi']['psi_{}_val'.format(xc)][:, 3]
                for p in range(rad_range):
                    plt.scatter(int(xc), psi_plot[p], color = 'r', label= 'psiN')
            plt.xlabel('Poloidal index')
            plt.ylabel('psiN')
            plt.title('poloidal index verses psiN')
        
        else:
            print('plot_psi function needs improvement')
    
    
    def plot_psi_surface(self):
        
        import numpy as np
        from mpl_toolkits.mplot3d import Axes3D  
        # Axes3D import has side effects, it enables using projection='3d' in add_subplot
        import matplotlib.pyplot as plt
        
        
        fig = plt.figure(figsize=(7,7))
        ax = fig.add_subplot(111, projection='3d')
        x = self.data['gfile']['gcomp']['gR']
        y = self.data['gfile']['gcomp']['gZ']
        X, Y = np.meshgrid(x, y)
        Z = np.array(self.data['gfile']['gcomp']['psiN'])
        z_one = np.array(np.ones((65, 65)))
        # Z = zs.reshape(X.shape)
        
        ax.plot_surface(X, Y, Z)
        
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        
        plt.show()

    
    
    def plot_RZ(self, pol_list):
        Attempt = self.data['dircomp']['Attempt']
        DRT = self.data['dirdata']['outputdir']['Output']
        DRT2 = self.data['dirdata']['outputdir']['Output2']
        XDIM = self.data['b2fgeo']['nx'] + 2
        YDIM = self.data['b2fgeo']['ny'] + 2
        
        
        
        RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                    usecols = (3)).reshape((YDIM, XDIM))
        VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                      usecols = (3)).reshape((YDIM,XDIM))
        
        
        if self.withshift == False and self.withseries == False:
            plt.figure(figsize=(7,7))
            for aa in pol_list:
                crloc = RadLoc[:, int(aa)]
                czloc = VertLoc[:, int(aa)]
                plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
            plt.xlabel('R: [m]')
            plt.ylabel('Z: [m]')
            
            
            plt.figure(figsize=(7,7))
            for ab in range(YDIM):
                for ac in pol_list:
                    crp = RadLoc[ab, int(ac)]
                    czp = VertLoc[ab, int(ac)]
                    plt.plot(crp, czp, 'o-', color = 'b', label = 'R&Zlocation')
            plt.xlabel('R: [m]')
            plt.ylabel('Z: [m]')
        else:
            print('plot_RZ function needs improvement')
        
    
    def plot_seperatrix(self):
        psi_1d = self.data['psi']['psival'][0, :]
        # self.data['psi']['psi1d'] = psi_1d
        
        pol_range = int(self.data['b2fgeo']['nx'] + 2)
        rad_range = int(self.data['b2fgeo']['ny'] + 2)
        
        
        index_low = []
        index_high = []
        index = np.zeros(2)
        for y in range(rad_range):
            if psi_1d[y] <= 1:
                index_low.append(y)
            if psi_1d[y] >= 1:
                index_high.append(y)
    
        
        index[0] = index_low[-1]
        index[1] = index_high[0]
        
        index_dic = {'index_low': index_low, 'index_high': index_high, 
                     'index': index}
        self.data['index'] = index_dic
        
                
            
        
        
            
        