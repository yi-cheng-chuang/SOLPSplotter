# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 16:55:28 2025

@author: ychuang
"""

import numpy as np



class target_dataload:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data



    def tarNTdata(self, itername, data_struc):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        dat_size = self.DF.data_size
        
        
        if withshift == True and withseries == False:
            
            if dat_size == 'small':
                nx = data_struc['nx']
                ny = data_struc['ny']
                b2fstate = self.data['b2fstate']
                ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
                Te_J = b2fstate['te'][1:nx+1, 1:ny+1]
                
                ev = 1.6021766339999999 * pow(10, -19)
                te_dat = Te_J / ev
                
                source = self.data['b2wdat']['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                vol = self.data['b2wdat']['vol'][1:nx+1, 1:ny+1]
                sx = np.divide(source, vol)
                
                neuden_dat = self.data['ft44']['dab2'][:, :, 0]
                # neuden_dat = np.transpose(data[:, :, 0])
                              
                psi_coord = self.data['psi']['psival'][1:ny+1, 1:nx+1]
            
            else:
                
                print('I do not use full data size for my study')
        
        
        elif withshift == True and withseries == False:
            
                
            if dat_size == 'small':
                nx = data_struc['nx']
                ny = data_struc['ny']
                b2fstate = self.data['b2fstate'][itername]
                ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
                Te_J = b2fstate['te'][1:nx+1, 1:ny+1]
                
                ev = 1.6021766339999999 * pow(10, -19)
                te_dat = Te_J / ev
                
                source = self.data['b2wdat'][itername]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                vol = self.data['b2wdat'][itername]['vol'][1:nx+1, 1:ny+1]
                sx = np.divide(source, vol)
                
                neuden_dat = self.data['ft44'][itername]['dab2'][:, :, 0]
                # neuden_dat = np.transpose(data[:, :, 0])
                              
                psi_coord = self.data['psi']['psival'][1:ny+1, 1:nx+1]
            
            else:
                
                print('I do not use full data size for my study')
            
        
        elif withshift == False and withseries == True:
            
            if self.DF.series_flag == 'twin_scan':
                
                
                nf = itername[0]
                tf = itername[1]
                    
                                  
                if dat_size == 'small':
                    nx = data_struc['nx']
                    ny = data_struc['ny']
                    b2fstate = self.data['b2fstate'][nf][tf]
                    ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
                    Te_J = b2fstate['te'][1:nx+1, 1:ny+1]
                    
                    ev = 1.6021766339999999 * pow(10, -19)
                    te_dat = Te_J / ev
                                    
                    psi_coord = self.data['psi']['psival'][1:ny+1, 1:nx+1]
                    
                    source = self.data['b2wdat'][nf][tf]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                    vol = self.data['b2wdat'][nf][tf]['vol'][1:nx+1, 1:ny+1]
                    sx = np.divide(source, vol)
                    
                    neuden_dat = self.data['ft44'][nf][tf]['dab2'][:, :, 0]
                    # neuden_dat = np.transpose(data[:, :, 0])
                
                else:
                    
                    print('I do not use full data size for my study')
                
            else:
                    
                if dat_size == 'small':
                    nx = data_struc['nx']
                    ny = data_struc['ny']
                    b2fstate = self.data['b2fstate'][itername]
                    ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
                    Te_J = b2fstate['te'][1:nx+1, 1:ny+1]
                
                
                    
                    ev = 1.6021766339999999 * pow(10, -19)
                    te_dat = Te_J / ev
                    
                    source = self.data['b2wdat'][itername]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                    vol = self.data['b2wdat'][itername]['vol'][1:nx+1, 1:ny+1]
                    sx = np.divide(source, vol)
                    
                    neuden_dat = self.data['ft44'][itername]['dab2'][:, :, 0]
                    # neuden_dat = np.transpose(data[:, :, 0])
                                  
                    psi_coord = self.data['psi']['psival'][1:ny+1, 1:nx+1]
            
        
        
            
            
        psi_dic = {'inner target': psi_coord[:, 0], 'outer target': psi_coord[:, nx-1]}
        ne_dic = {'inner target': ne_dat[0, :], 'outer target': ne_dat[nx-1, :]}
        te_dic = {'inner target': te_dat[0, :], 'outer target': te_dat[nx-1, :]}
        sx_dic = {'inner target': sx[0, :], 'outer target': sx[nx-1, :]}
        neuden_dic = {'inner target': neuden_dat[0, :], 'outer target': neuden_dat[nx-1, :]}
        
        
        target_dic = {'psiN': psi_dic, 'ne': ne_dic, 'te': te_dic, 'source': sx_dic, 
                      'neuden': neuden_dic}
        
        
        return target_dic
    
    
    

    
    
