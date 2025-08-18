# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 23:34:10 2025

@author: ychuang
"""



from load_coordinate.load_coord_method import load_coordgeo_method
from fit_data.fitting_method import fit_method_collection
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from scipy.interpolate import RectBivariateSpline


class core_boundary_calc:
    
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    
    def calc_core_boundary_length_method(self, cr, cz, plot):
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
    
    
    
    
    def calc_core_boundary_length(self, plot):
        
        if self.DF.DEV == 'cross_machine':
            
            if self.DF.Dnames == 'mastu_mast':
                
                dev_list = ['mast', 'mastu']
                
                for aa in dev_list:
                    
                    if aa == 'mast':
                        
                        rad = self.data['{}_grid'.format(aa)]['RadLoc'][0, 25:73]
                        vert = self.data['{}_grid'.format(aa)]['VertLoc'][0, 25:73]
                        
                        arc_len, inter_p = self.calc_core_boundary_length_method(cr = rad, cz = vert, plot = plot)
                        
                        hx_arr = self.data['{}_b2fgeo'.format(aa)]['hx'][25:73, 0]
                        sum_hx = sum(hx_arr)
                        
                        print('The {} arc length is: {:.3f}'.format(aa, arc_len[-1]))
                        print('The {} summation hx is:'.format(aa))
                        print(sum_hx)
                    
                    
                    elif aa == 'mastu':
                        
                        rad_a = self.data['{}_grid'.format(aa)]['RadLoc'][0, 13:36]
                        vert_a = self.data['{}_grid'.format(aa)]['VertLoc'][0, 13:36]
                        rad_b = self.data['{}_grid'.format(aa)]['RadLoc'][0, 62:85]
                        vert_b = self.data['{}_grid'.format(aa)]['VertLoc'][0, 62:85]
                        
                        arc_len_a, inter_p_a = self.calc_core_boundary_length_method(cr = rad_a, cz = vert_a, plot = plot)
                        arc_len_b, inter_p_b = self.calc_core_boundary_length_method(cr = rad_b, cz = vert_b, plot = plot)
                        
                        print(type(arc_len_a))
                        arc_len = arc_len_a[-1] + arc_len_b[-1]
                        
                        hx_arr_a = self.data['{}_b2fgeo'.format(aa)]['hx'][13:36, 0]
                        sum_hx_a = sum(hx_arr_a)
                        
                        hx_arr_b = self.data['{}_b2fgeo'.format(aa)]['hx'][62:85, 0]
                        sum_hx_b = sum(hx_arr_b)
                        
                        sum_hx = sum_hx_a + sum_hx_b
                        
                        print('The {} arc length is: {:.3f}'.format(aa, arc_len))
                        print('The {} summation hx is:'.format(aa))
                        print(sum_hx)
                    
                    
            
    def calc_core_boundary_flux(self):
        
        
        if self.DF.DEV == 'cross_machine':
            
            if self.DF.Dnames == 'mastu_mast':
                
                dev_list = ['mast', 'mastu']
                
                for aa in dev_list:
                    
                    if aa == 'mast':
                        
                        fnay = self.data['{}_b2wdat'.format(aa)]['b2npc_fnays'][0][25:73, 1]
                        fnay_tot = sum(fnay)
                        
                        fhey = self.data['{}_b2wdat'.format(aa)]['b2nph9_fhey'][25:73, 1]
                        fhey_tot = sum(fhey)
                        
                        hx = self.data['{}_b2wdat'.format(aa)]['hx'][25:73, 1]
                        hz = self.data['{}_b2wdat'.format(aa)]['hz'][25:73, 1]
                        area_sum = sum(np.multiply(hx, hz))
                        
                        
                        
                        print('the total radial particle flux at core boundary is {:.3e}'.format(fnay_tot))
                        print('the total radial heat flux at core boundary is {:.3e}'.format(fhey_tot))
                        print('the total area at core boundary is {:.3e}'.format(area_sum))
                        
                        
                        mast_dat_dic = {'mast_fnay': fnay_tot, 'mast_fhey': fhey_tot, 'mast_area': area_sum}
                    
                    
                    elif aa == 'mastu':
                        
                        fnay_a = self.data['{}_b2wdat'.format(aa)]['b2npc_fnays'][0][13:36, 1]
                        fnay_tot_a = sum(fnay_a)
                        
                        fnay_b = self.data['{}_b2wdat'.format(aa)]['b2npc_fnays'][0][62:85, 1]
                        fnay_tot_b = sum(fnay_b)
                        
                        fnay_tot = fnay_tot_a + fnay_tot_b
                        
                        
                        
                        fhey_a = self.data['{}_b2wdat'.format(aa)]['b2nph9_fhey'][13:36, 1]
                        fhey_tot_a = sum(fhey_a)
                        
                        fhey_b = self.data['{}_b2wdat'.format(aa)]['b2nph9_fhey'][62:85, 1]
                        fhey_tot_b = sum(fhey_b)
                        
                        fhey_tot = fhey_tot_a + fhey_tot_b
                        
                        
                        
                        hx_a = self.data['{}_b2wdat'.format(aa)]['hx'][13:36, 1]
                        hz_a = self.data['{}_b2wdat'.format(aa)]['hz'][13:36, 1]
                        area_sum_a = sum(np.multiply(hx_a, hz_a))
                        
                        hx_b = self.data['{}_b2wdat'.format(aa)]['hx'][62:85, 1]
                        hz_b = self.data['{}_b2wdat'.format(aa)]['hz'][62:85, 1]
                        area_sum_b = sum(np.multiply(hx_b, hz_b))
                        
                        area_sum = area_sum_a + area_sum_b
                        
                        print('the total radial particle flux at core boundary is {:.3e}'.format(fnay_tot))
                        print('the total radial heat flux at core boundary is {:.3e}'.format(fhey_tot))
                        print('the total area at core boundary is {:.3e}'.format(area_sum))
            
            
                        mastu_dat_dic = {'mastu_fnay_inner': fnay_tot_a, 'mastu_fnay_outer': fnay_tot_b, 
                                         'mastu_fhey_inner': fhey_tot_a, 'mastu_fhey_outer': fhey_tot_b,
                                         'mastu_area_inner': area_sum_a, 'mastu_area_outer': area_sum_b, }
            
            
            
            
                mast_area = mast_dat_dic['mast_area']
                mastu_area_inner = mastu_dat_dic['mastu_area_inner']
                mastu_area_outer = mastu_dat_dic['mastu_area_outer']
                
                ratio_inner = mastu_area_inner / mast_area
                ratio_outer = mastu_area_outer / mast_area
                
                mastu_pflux_inner = mast_dat_dic['mast_fnay']* ratio_inner
                mastu_pflux_outer = mast_dat_dic['mast_fnay']* ratio_outer
                mastu_hflux_inner = mast_dat_dic['mast_fhey']* ratio_inner
                mastu_hflux_outer = mast_dat_dic['mast_fhey']* ratio_outer
                
                
                print('the set radial inner particle flux at core boundary for mastu is {:.3e}'.format(mastu_pflux_inner))
                print('the set radial outer particle flux at core boundary for mastu is {:.3e}'.format(mastu_pflux_outer))
                print('the set radial inner heat flux at core boundary for mastu is {:.3e}'.format(mastu_hflux_inner))
                print('the set radial outer heat flux at core boundary for mastu is {:.3e}'.format(mastu_hflux_outer))
                

            
            
            
        





        
        