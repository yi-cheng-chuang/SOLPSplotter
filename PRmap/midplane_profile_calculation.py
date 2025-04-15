# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 16:12:28 2025

@author: ychuang
"""

from load_simulation_data.read_B2simulation_data import load_B2simu_data
from load_directory.load_dirdata_method import load_dir_method
import numpy as np




class midplane_radial:
    
    def __init__(self, DF, data, lbd: load_B2simu_data, ldm: load_dir_method):
        
        self.DF = DF
        self.data = data
        self.lbd = lbd
        self.ldm = ldm
    
    
    def calc_midplane_profile_method(self, itername, data_struc):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        nx = data_struc['nx']
        ny = data_struc['ny']
        
        if withshift == False and withseries == False:
            
            b2fstate = self.data['b2fstate']
                
            if self.DF.data_size == 'small':
                
                data = self.data['ft44']['dab2']
                neu_pro = np.transpose(data[:, :, 0])
                
                source = self.data['b2wdat']['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                vol = self.data['b2wdat']['vol'][1:nx+1, 1:ny+1]
                sx = np.divide(source, vol)
            
            else:
                
                print('I do not use full data size for my study')
            
            

        elif withshift == True and withseries == False:
            
            
            b2fstate = self.data['b2fstate'][itername]
                
            if self.DF.data_size == 'small':
                
                data = self.data['ft44'][itername]['dab2']
                neu_pro = np.transpose(data[:, :, 0])
                
                source = self.data['b2wdat'][itername]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                vol = self.data['b2wdat'][itername]['vol'][1:nx+1, 1:ny+1]
                sx = np.divide(source, vol)
            
            else:
                
                print('I do not use full data size for my study')
        
        elif withshift == False and withseries == True:
            
            series_flag = self.DF.series_flag
            
            
            if series_flag == 'twin_scan':
                
                nf = itername[0]
                tf = itername[1]
                
                b2fstate = self.data['b2fstate'][nf][tf]
                    
                if self.DF.data_size == 'small':
                    
                    data = self.data['ft44'][nf][tf]['dab2']
                    neu_pro = np.transpose(data[:, :, 0])
                    
                    source = self.data['b2wdat'][nf][tf]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                    vol = self.data['b2wdat'][nf][tf]['vol'][1:nx+1, 1:ny+1]
                    sx = np.divide(source, vol)
                
                else:
                    
                    print('I do not use full data size for my study')
                
            else:
                
                b2fstate = self.data['b2fstate'][itername]
                    
                if self.DF.data_size == 'small':
                    
                    data = self.data['ft44'][itername]['dab2']
                    neu_pro = np.transpose(data[:, :, 0])
                    
                    source = self.data['b2wdat'][itername]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                    vol = self.data['b2wdat'][itername]['vol'][1:nx+1, 1:ny+1]
                    sx = np.divide(source, vol)
                
                else:
                    
                    print('I do not use full data size for my study')
      
            
        if self.DF.data_size == 'small':
            ne_pro = b2fstate['ne'][1:nx+1, 1:ny+1].transpose()
            Te_J = b2fstate['te'][1:nx+1, 1:ny+1].transpose()
        
        else:
            
            print('I do not use full data size for my study')
            
        
        
        ev = 1.6021766339999999 * pow(10, -19)
        te_pro = Te_J / ev
        
        # neu_pro = np.transpose(data[:, :, 0])
        
        
        
        if withshift == False and withseries == False:
            
                           
            if self.DF.data_size == 'small':
                
                psi_coord = self.data['midplane_calc']['psi_solps_mid'][1:ny+1]
                weight = self.data['midplane_calc']['weight'][1:ny+1]
            
            else:
                
                print('I do not use full data size for my study')
                
                
        
        elif withshift == True and withseries == False:
        
                
            if self.DF.data_size == 'small':
                
                psi_coord = self.data['midplane_calc'][itername]['psi_solps_mid'][1:ny+1]
                weight = self.data['midplane_calc'][itername]['weight'][1:ny+1]
            
            else:
                
                print('I do not use full data size for my study')
            
            
        elif withshift == False and withseries == True:
            
                
            if self.DF.data_size == 'small':
                
                psi_coord = self.data['midplane_calc']['psi_solps_mid'][1:ny+1]
                weight = self.data['midplane_calc']['weight'][1:ny+1]
            
            else:
                
                print('I do not use full data size for my study')
        
        else:
            print('nete_TSplotmethod geo cut has a bug!')
        
        weight_B = np.ones(len(weight))- weight
            
        
        if self.DF.data_size == 'small':
            
            ap = self.data['midplane_calc']['average_pair']
            pair = (ap[0]-1, ap[1]-1)
        
        else:
            
            print('I do not use full data size for my study')
        
        
        mid_ne_pro = np.multiply(ne_pro[:, pair[0]], weight) + np.multiply(ne_pro[:, pair[1]], weight_B)
        mid_te_pro = np.multiply(te_pro[:, pair[0]], weight) + np.multiply(te_pro[:, pair[1]], weight_B)
        mid_neu_pro = np.multiply(neu_pro[:, pair[0]], weight) + np.multiply(neu_pro[:, pair[1]], weight_B)
        mid_S_pro = np.multiply(sx[pair[0], :], weight) + np.multiply(sx[pair[1], :], weight_B)
        
        
        midplane_profile_dic = {'psiN': psi_coord, 'mid_ne': mid_ne_pro, 'mid_te': mid_te_pro,
                                'mid_nd': mid_neu_pro, 'mid_S': mid_S_pro}
        
        
        return midplane_profile_dic

    
    def oneDscan_midplane(self, iterlist, dat_dic, dat_struc):
        
        
        series_compare = self.DF.series_compare
        
        if series_compare == True:
            
            print('we will improve this in the future!')
        
        else:
            
            
            for aa in iterlist:

                midprofiles = self.calc_midplane_profile_method(itername = aa, data_struc = dat_struc)
                
                dat_dic[aa] = midprofiles
                
                
        return dat_dic
    
        
    def twoDscan_midplane(self, iterlist, iterlist_a, iterlist_b, dat_dic, dat_struc):
        

            
        for tp in iterlist:
            aa = tp[0]
            ab = tp[1]
            
            midprofiles = self.calc_midplane_profile_method(itername = tp, data_struc = dat_struc)
            dat_dic[aa][ab] = midprofiles
  
        return dat_dic
    
    
       
    def calc_midplane_profile(self):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        
        if withshift == False and withseries == False:
            
            nx = self.data['b2fgeo']['nx']
            ny = self.data['b2fgeo']['ny']
            dat_struc = {'nx': nx, 'ny': ny}
                
            midprofile_dic = self.calc_midplane_profile_method(itername = None, data_struc = dat_struc)
            
            self.data['midplane_profile'] = midprofile_dic
        
        elif withshift == True and withseries == False:
            
            
            midprofile_dic = {}
            
            
            for aa in self.data['dircomp']['multi_shift']:
                
                
                nx = self.data['b2fgeo'][aa]['nx']
                ny = self.data['b2fgeo'][aa]['ny']
                dat_struc = {'nx': nx, 'ny': ny}
                
                
                midprofiles = self.calc_midplane_profile_method(itername = aa, data_struc = dat_struc)
                
                midprofile_dic[aa] = midprofiles
            
            self.data['midplane_profile'] = midprofile_dic
        
        elif withshift == False and withseries == True:
            
            midprofile_dic = {}
            
            nx = self.data['b2fgeo']['nx']
            ny = self.data['b2fgeo']['ny']
            dat_struc = {'nx': nx, 'ny': ny}
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            if self.DF.series_flag == 'twin_scan':
                
                ds_key, ts_key = self.lbd.twokeylists(printvalue= False)
                
                mid_dic = self.ldm.two_layer_dic(key_a = ds_key, key_b = ts_key)
                
                midprofile_dic = self.twoDscan_midplane(iterlist = scan, iterlist_a = ds_key, iterlist_b = ts_key, 
                                                        dat_dic = mid_dic, dat_struc = dat_struc)
                
            
            else:
                
                mid_dic = {}
                midprofile_dic = self.oneDscan_midplane(iterlist = scan, dat_dic = mid_dic, dat_struc = dat_struc)
            
            self.data['midplane_profile'] = midprofile_dic
        
        elif withshift == True and withseries == True:
            print('calc_midplane_profile is not there yet!')
        
        
        else:
            print('calc_midplane_profile has a bug')
    



"""
backup

def nete_midprof_cp(self, itername, cptag, data_struc):
    
    withshift = self.DF.withshift
    withseries = self.DF.withseries
    

    if withshift == True and withseries == False:
        
        
        
        b2fstate = self.data['b2fstate'][itername][cptag]
        
        if data_struc['size'] == 'full':
            neu_pro = self.data['outputdata']['NeuDen'][itername]
        elif data_struc['size'] == 'small':
            data = self.data['ft44'][itername][cptag]['dab2']
            neu_pro = np.transpose(data[:, :, 0])
    
    
    
    nx = data_struc['nx']
    ny = data_struc['ny']
    
    if data_struc['size'] == 'full':
        ne_pro = b2fstate['ne'].transpose()
        Te_J = b2fstate['te'].transpose()
        
    elif data_struc['size'] == 'small':
        ne_pro = b2fstate['ne'][1:nx+1, 1:ny+1].transpose()
        Te_J = b2fstate['te'][1:nx+1, 1:ny+1].transpose()
        
    
    
    ev = 1.6021766339999999 * pow(10, -19)
    te_pro = Te_J / ev
    
    # neu_pro = np.transpose(data[:, :, 0])
    
    if withshift == True and withseries == False:
        
        if data_struc['size'] == 'full':
            psi_coord = self.data['midplane_calc'][itername]['psi_solps_mid']
            weight = self.data['midplane_calc'][itername]['weight']
            
        elif data_struc['size'] == 'small':
            psi_coord = self.data['midplane_calc'][itername]['psi_solps_mid'][1:ny+1]
            weight = self.data['midplane_calc'][itername]['weight'][1:ny+1]
    else:
        print('nete_midprof_cp is not there yet!')
    
    weight_B = np.ones(len(weight))- weight
         
    mid_ne_pro = np.multiply(ne_pro[:, 58], weight) + np.multiply(ne_pro[:, 60], weight_B)
    mid_te_pro = np.multiply(te_pro[:, 58], weight) + np.multiply(te_pro[:, 60], weight_B)
    mid_neu_pro = np.multiply(neu_pro[:, 58], weight) + np.multiply(neu_pro[:, 60], weight_B)
    
    return psi_coord, mid_ne_pro, mid_te_pro, mid_neu_pro

"""
