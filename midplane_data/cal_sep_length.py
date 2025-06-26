# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 20:39:10 2025

@author: ychuang
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


class seplength_calculator:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    def calc_sepx_length(self, cr, cz, plot):
        # Define some points:
        p = len(cr)
        points = np.zeros((p, 2))
        points[:, 0] = cr
        points[:, 1] = cz
        # points = np.array([[0, 1, 8, 2, 2],
        #                    [1, 0, 6, 7, 2]]).T  # a (nbre_points x nbre_dim) array
        
        "Linear length along the line:"
        
        # crzdiff = np.diff(points, axis=0)
        # crzsum = np.sum( np.diff(points, axis=0)**2, axis=1 )
        distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )))
        distance = np.insert(distance, 0, 0)/distance[-1]
        
        arclength = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )))
        arclength = np.insert(arclength, 0, 0)
        
        
        "Interpolation for different methods:"
        
        interpolations_methods = ['slinear', 'quadratic', 'cubic']
        alpha = np.linspace(0, 1, 38)
        
        interpolated_points = {}
        for method in interpolations_methods:
            interpolator =  interp1d(distance, points, kind=method, axis=0)
            interpolated_points[method] = interpolator(alpha)
        
        "Graph:"
        
        if plot:
            plt.figure()
            for method_name, curve in interpolated_points.items():
                plt.plot(*curve.T, '-', label=method_name);
            
            plt.plot(*points.T, 'ok', label='original points');
            plt.axis('equal'); plt.legend(); plt.xlabel('r'); plt.ylabel('z');
            
        

        return arclength, interpolated_points




    def sep_length(self, direct):
        
        
        if self.DF.withshift:
            
            sepl_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                
                
                nx = self.data['b2fgeo'][aa]['nx']
                ny = self.data['b2fgeo'][aa]['ny']
                
                rad_grid = self.data['grid']['RadLoc'][aa]
                vert_grid = self.data['grid']['VertLoc'][aa]
                
                
                # R_loc = rad_grid[18, :][::-1]
                
                # print(R_loc)
                # print(R_loc)
                
                
                # Z_loc = vert_grid[18, :][::-1]
                
                if direct == 'in_to_out':
                    
                    R_loc = rad_grid[18, 1:nx+1]
                    Z_loc = vert_grid[18, 1:nx+1]
                
                elif direct == 'out_to_in':
                    
                    R_loc = rad_grid[18, 1:nx+1][::-1]               
                    Z_loc = vert_grid[18, 1:nx+1][::-1]
                    
            
                          
                arclength, interpolated_points = self.calc_sepx_length(cr = R_loc, 
                                                        cz = Z_loc, plot = False)
            
                sepl_dic[aa] = {'sepl': arclength, 'interp': interpolated_points}
            
            self.data['sepl'] = sepl_dic
        
        else:
            pass
            
            
        
        
        



    

