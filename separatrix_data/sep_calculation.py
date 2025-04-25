# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 00:18:11 2025

@author: ychuang
"""

from load_coordinate.load_coord_method import load_coordgeo_method
from fit_data.fitting_method import fit_method_collection
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np


class poloidal_separatrix_calc:
    
    
    def __init__(self, DF, data, lcm: load_coordgeo_method, fmc: fit_method_collection):
        
        self.DF = DF
        self.data = data
        self.fmc = fmc
        self.lcm = lcm



#-----------------------------------------------------------------------------

    #calculate the poloidal angle function    
        
    def calc_pol_angle_method(self, RadLoc, VertLoc, gfile, shift, pol_list, plot_angle):
        
        
        mag_axis_z = gfile['zmaxis']
        # print(mag_axis_z)
        
        mag_axis_r = gfile['rmaxis'] + shift
        
                         
        xprloc = RadLoc[:, 72][19]
        xpzloc = VertLoc[:, 72][19]
        
        xpoint_R = xprloc - mag_axis_r
        xpoint_Z = xpzloc - mag_axis_z
        xpoint_angle = np.arctan2(xpoint_Z, xpoint_R)*180 / np.pi
        

        angle_list = []
        
        for pol_index in pol_list:
            rloc = RadLoc[:, int(pol_index)][19]
            zloc = VertLoc[:, int(pol_index)][19]
            
            x = rloc - mag_axis_r
            y = zloc - mag_axis_z
            
            
            angle = np.arctan2(y, x)*180 / np.pi
            if angle <= xpoint_angle:
                angle = angle + 360
            angle_list.append(angle)
            
        pol_loc = np.asarray(list(map(int, pol_list)))
        
        if plot_angle: 
            plt.figure(figsize=(7,7))
            plt.plot(pol_loc, angle_list, 'o', color = 'r', label= 'poloidal angle')
            plt.xlabel('poloidal index')
            plt.title('poloidal angle verses poloidal index')
            plt.legend()
                
        return angle_list, xpoint_angle
    
    
    
    def calc_pol_angle(self, pol_list, plot_angle):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == False:
            
            rad_grid = self.data['grid']['RadLoc']
            vert_grid = self.data['grid']['VertLoc']     
            g_data = self.data['gfile']['g']
            shift_val = self.data['dircomp']['shift_value']
            
                        
            angle_list, xpoint_angle = self.calc_pol_angle_method(RadLoc = rad_grid, VertLoc = vert_grid, 
                gfile = g_data, shift = shift_val, plot_angle = plot_angle, pol_list = pol_list)
            
            print('xpoint angle is {:.2f}'.format(xpoint_angle))
            
            angle_dic = {'angle_list': angle_list, 'xpoint_angle': xpoint_angle}
            
            self.data['angle'] = angle_dic
            
        elif withshift == True and withseries == False:        
            angle_dic = {}
            angle_list_dic = {}
            xpoint_angle_dic = {}
            index_angle_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                
                rad_grid = self.data['grid']['RadLoc'][aa]
                vert_grid = self.data['grid']['VertLoc'][aa]    
                g_data = self.data['gfile']['g']
                shift_val = self.data['dircomp']['shift_dic'][aa]
                
                angle_list, xpoint_angle = self.calc_pol_angle_method(RadLoc = rad_grid, VertLoc = vert_grid, 
                    gfile = g_data, shift = shift_val, plot_angle = plot_angle, pol_list = pol_list)
                
                angle_list_dic[aa] = angle_list
                xpoint_angle_dic[aa] = xpoint_angle
                
                
                print('xpoint angle is {:.2f} for {} case'.format(xpoint_angle, aa))
                
                   
            angle_dic = {'angle_list': angle_list_dic, 'xpoint_angle': xpoint_angle_dic}
            
            self.data['angle'] = angle_dic
            
           
        elif withshift == False and withseries == True:
            
            den_list = list(self.data['dircomp']['Attempt'].keys())
            aa = den_list[0]
            
            rad_grid = self.data['grid']['RadLoc']
            vert_grid = self.data['grid']['VertLoc']  
            g_data = self.data['gfile']['g']
            shift_val = self.data['dircomp']['shift_value']
            
            
            angle_list, xpoint_angle = self.calc_pol_angle_method(RadLoc = rad_grid, VertLoc = vert_grid, 
                gfile = g_data, shift = shift_val, plot_angle = plot_angle, pol_list = pol_list)
            
            angle_dic = {'angle_list': angle_list, 'xpoint_angle': xpoint_angle}
            
            self.data['angle'] = angle_dic
        
        elif withshift == True and withseries == True:
            print('calc_pol_angle function is not there yet!')
        
        else:
            print('calc_pol_angle function has a bug')
    
            


    # Poloidally create weight for separatrix


    def calc_polsep_method(self, pol_list, psiN):
         
        # pol_range = int(geo['nx'] + 2)
        # rad_range = int(geo['ny'] + 2)
        
        RadLoc = self.data['grid']['RadLoc']
 
        # psiN = self.data['psi']['psival'][aa][1:37, 1:97]
        
        "Calculate the weight for psiN = 1"
        
        "Calculate weight"
        
        st = int(pol_list[0])
        ed = int(pol_list[-1]) + 1
        
        high_psi = psiN[20, st:ed]
        low_psi = psiN[18, st:ed]
        
        list_len = len(high_psi)
    
        weight_mid = np.zeros(list_len)
        for x in range(list_len):
            weight_mid[x] = (1 - low_psi[x])/ (high_psi[x] - low_psi[x])
        

        return weight_mid
    
    
    def calc_polsep(self, pol_list):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == False:
            print('calc_polsep function is not there yet!')
        
        
        elif withshift == True and withseries == False:        
            

            weight_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                
                psiN = self.data['psi']['psival'][aa][1:37, 1:97]
                
                               
                weight_mid = self.calc_polsep_method(pol_list = pol_list, psiN = psiN)
                
                weight_dic[aa] = weight_mid
            
            
            self.data['polpsi_weight'] = weight_dic
            
           
        elif withshift == False and withseries == True:
            
            print('calc_polsep function is not there yet!')
            
            
        elif withshift == True and withseries == True:
            print('calc_RRsep function is not there yet!')
        
        else:
            print('calc_RRsep function has a bug')








