# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 19:31:29 2025

@author: ychuang
"""


from SOLPS_input.header import *


class midplane_AM_NT:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    def AM_NT_midprof(self, itername, AM_flag):
        
        if AM_flag == 'atom':
            
            den = 'dab2'
            temp = 'tab2'
            
        elif AM_flag == 'mol':
            
            den = 'dmb2'
            temp = 'tmb2'
        
        ev = 1.6021766339999999 * pow(10, -19)
        
        
        if self.withshift == True and self.withseries == False:
            
            neu_data = self.data['ft44'][itername][den]
            neu_pro = np.transpose(neu_data[:, :, 0])
            atom_temp_data = self.data['ft44'][itername][temp]
            
            atom_temp = np.transpose(atom_temp_data[:, :, 0])
            atom_temp_pro = atom_temp / ev
        
        elif self.withshift == False and self.withseries == True:
            
            if self.series_flag == 'twin_scan':
                
                nf = itername[0]
                tf = itername[1]
                
                neu_data = self.data['ft44'][nf][tf][den]
                neu_pro = np.transpose(neu_data[:, :, 0])
                atom_temp_data = self.data['ft44'][nf][tf][temp]
                atom_temp = np.transpose(atom_temp_data[:, :, 0])
                atom_temp_pro = atom_temp / ev
                              
            else:
                
                neu_data = self.data['ft44'][itername][den]
                neu_pro = np.transpose(neu_data[:, :, 0])
                atom_temp_data = self.data['ft44'][itername][temp]
                atom_temp = np.transpose(atom_temp_data[:, :, 0])
                atom_temp_pro = atom_temp / ev
        
        if self.withshift == True and self.withseries == False:
        
            leftcut = self.data['b2fgeo'][itername]['leftcut'][0]
            rightcut = self.data['b2fgeo'][itername]['rightcut'][0]
            weight = self.data['midplane_calc'][itername]['weight']
            psi_coord = self.data['midplane_calc'][itername]['psi_solps_mid']
        
        elif self.withshift == False and self.withseries == True:
        
            leftcut = self.data['b2fgeo']['leftcut'][0]
            rightcut = self.data['b2fgeo']['rightcut'][0]
            weight = self.data['midplane_calc']['weight'][1:37]
            psi_coord = self.data['midplane_calc']['psi_solps_mid'][1:37]
        
        else:
            print('NeuDen_plotmethod, please check withshift and withseries flag')
        
        weight_B = np.ones(len(weight))- weight
        
        
        if dat_size == 'small':
            
            ap = self.data['midplane_calc']['average_pair']
            pair = (ap[0]-1, ap[1]-1)
        
        else:
            
            print('I do not use full data size for my study')
        
        mid_neu_pro = np.multiply(neu_pro[:, 58], weight) + np.multiply(neu_pro[:, 60], weight_B)
        mid_atom_temp_pro = np.multiply(atom_temp_pro[:, 58], weight) + np.multiply(atom_temp_pro[:, 60], weight_B)
        
        return psi_coord, mid_neu_pro, mid_atom_temp_pro








    