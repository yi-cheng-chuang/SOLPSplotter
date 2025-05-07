# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 12:15:25 2025

@author: ychuang
"""


import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np


class B2_boundary_contour:
    
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    
    def plot_B2boundary(self):
        
        
        X = self.data['grid']['RadLoc']
        Y = self.data['grid']['VertLoc']
            
        vessel = self.data['vessel']
        LFS_Rsep = self.data['calc_minor_radius']['LFS']['LFS_R_sep']
        HFS_Rsep = self.data['calc_minor_radius']['HFS']['HFS_R_sep']
        majloc = self.data['calc_minor_radius']['majR_loc']

        # Plot the curvilinear grid
        sep_R_up = X[18, :]
        sep_Z_up = Y[18, :]
        
        sep_R_down = X[17, :]
        sep_Z_down = Y[17, :]
        
        sep_R = 0.5*(sep_R_up + sep_R_down)
        sep_Z = 0.5*(sep_Z_up + sep_Z_down)
        
        

        cut_filter = True
        plot_vessel = True

        plt.figure()

        if cut_filter:
            plt.plot(X[0, :13], Y[0, :13], 'b', linewidth=1)
            plt.plot(X[0, 13:36], Y[0, 13:36], 'b', linewidth=1)
            plt.plot(X[0, 36:47], Y[0, 36:47], 'b', linewidth=1)
            plt.plot(X[0, 49:62], Y[0, 49:62], 'b', linewidth=1)
            plt.plot(X[0, 62:85], Y[0, 62:85], 'b', linewidth=1)
            plt.plot(X[0, 85:], Y[0, 85:], 'b', linewidth=1)
            plt.plot(X[-1, :47], Y[-1, :47], 'b', linewidth=1, label = 'B2.5 grid boundary')
            plt.plot(X[-1, 49:], Y[-1, 49:], 'b', linewidth=1)
            # plt.plot(sep_R, sep_Z,'-', color = 'green', linewidth=1, label = 'separatrix')
            plt.plot(sep_R, sep_Z,'-', color = 'green', linewidth=1, label= 'separatrix')
            
                
        else:

            plt.plot(X[0, :], Y[0, :], 'b', linewidth=1)
            plt.plot(X[-1, :], Y[-1, :], 'b', linewidth=1)


        if plot_vessel:
            
            plt.plot(vessel[:,0]/1000, vessel[:,1]/1000, color = 'black', linewidth = 2)

        else:
            pass

        
        # Plot a vertical line at x=2
        plt.axvline(x= round(LFS_Rsep, 2), color='red', linestyle='--')
        plt.axvline(x= round(HFS_Rsep, 2), color='red', linestyle='--')
        plt.axvline(x= round(majloc, 2), color='red', linestyle='--')
        
        # Add legend and show
        plt.legend()



        
        plt.gca().set_aspect('equal')
        plt.xlabel("R: [m]")
        plt.ylabel("Z: [m]")
        plt.show()
        
        
    
    