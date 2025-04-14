# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 20:05:34 2025

@author: ychuang
"""

"""
AMS is atomic and molecular neutrals and Source!



"""
from SOLPS_input.header import *




class midplane_ndsource_withcut:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
       

    def ndSmid_withcut_method(self, itername, data_struc, cut_range):

        dat_size = self.DF.data_size
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        if withshift == False and withseries == False:
            
            if dat_size == 'small':
                nx = data_struc['nx']
                ny = data_struc['ny']
                data = self.data['ft44']['dab2']
                neu_pro = np.transpose(data[:, :, 0])
                weight = self.data['midplane_calc']['weight'][1:ny+1]
                psi_coord = self.data['midplane_calc']['psi_solps_mid'][1:ny+1]
                
                source = self.data['b2wdat']['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                vol = self.data['b2wdat']['vol'][1:nx+1, 1:ny+1]
                sx = np.divide(source, vol)
            
            else:
                
                print('I do not use full data size for my study')

        elif withshift == True and withseries == False:
            
            if dat_size == 'small':
                nx = data_struc['nx']
                ny = data_struc['ny']
                data = self.data['ft44'][itername]['dab2']
                neu_pro = np.transpose(data[:, :, 0])
                weight = self.data['midplane_calc'][itername]['weight'][1:ny+1]
                psi_coord = self.data['midplane_calc'][itername]['psi_solps_mid'][1:ny+1]
                
                
                source = self.data['b2wdat'][itername]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                vol = self.data['b2wdat'][itername]['vol'][1:nx+1, 1:ny+1]
                sx = np.divide(source, vol)
            
            else:
                
                print('I do not use full data size for my study')
            
        
        elif withshift == False and withseries == True:
            
            if self.DF.series_flag == 'twin_scan':
                
                nf = itername[0]
                tf = itername[1]
                    
                                  
                if dat_size == 'small':
                    nx = data_struc['nx']
                    ny = data_struc['ny']
                    data = self.data['ft44'][nf][tf]['dab2']
                    neu_pro = np.transpose(data[:, :, 0])
                    weight = self.data['midplane_calc']['weight'][1:ny+1]
                    psi_coord = self.data['midplane_calc']['psi_solps_mid'][1:ny+1]
                    
                    source = self.data['b2wdat'][nf][tf]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                    vol = self.data['b2wdat'][nf][tf]['vol'][1:nx+1, 1:ny+1]
                    sx = np.divide(source, vol)

                
                else:
                    
                    print('I do not use full data size for my study')
                
            else:
                
                if dat_size == 'small':
                    nx = data_struc['nx']
                    ny = data_struc['ny']
                    data = self.data['ft44'][itername]['dab2']
                    neu_pro = np.transpose(data[:, :, 0])
                    weight = self.data['midplane_calc']['weight'][1:ny+1]
                    psi_coord = self.data['midplane_calc']['psi_solps_mid'][1:ny+1]
                    
                    
                    source = self.data['b2wdat'][itername]['b2npc_sna'][0][1:nx+1, 1:ny+1]                
                    vol = self.data['b2wdat'][itername]['vol'][1:nx+1, 1:ny+1]
                    sx = np.divide(source, vol)
                
                else:
                    
                    print('I do not use full data size for my study')
            

        weight_B = np.ones(len(weight))- weight
            
        
        if dat_size == 'small':
            
            ap = self.data['midplane_calc']['average_pair']
            pair = (ap[0]-1, ap[1]-1)
        
        else:
            
            print('I do not use full data size for my study')
        
        mid_neu_pro = np.multiply(neu_pro[:, pair[0]], weight) + np.multiply(neu_pro[:, pair[1]], weight_B)
        mid_S_pro = np.multiply(sx[:, pair[0]], weight) + np.multiply(sx[:, pair[1]], weight_B)
        
        psi_list = []
        nd_list = []
        S_list = []
        st = cut_range[0]
        ed = cut_range[1]
        
        for ind, coord in enumerate(psi_coord):
            
            if coord >= st and coord <= ed:
                psi_list.append(coord)
                nd_list.append(mid_neu_pro[ind])
                S_list.append(mid_S_pro[ind])
        
        
        ndSmid_cut_dic = {'psi_cut': psi_list, 'nd_cut': nd_list, 'source_cut': S_list}
        
        
        
        return ndSmid_cut_dic
        
    
    
    def oneDscan_ndScut(self, iterlist, dat_dic, dat_struc, cut_range):
        
        
        series_compare = self.DF.series_compare
        
        if series_compare == True:
            
            print('we will improve this in the future!')
        
        else:
            
            
            for aa in iterlist:

                ndcut = self.ndSmid_withcut_method(itername = aa, data_struc = dat_struc, cut_range = cut_range)
                
                dat_dic[aa] = ndcut
                
                
        return dat_dic
    
    
    
    def twoDscan_ndScut(self, iterlist, iterlist_a, iterlist_b, dat_dic, dat_struc, cut_range):

        
        for tp in iterlist:
            aa = tp[0]
            ab = tp[1]
            
            ndcut = self.ndSmid_withcut_method(itername = tp, data_struc = dat_struc, cut_range = cut_range)
            dat_dic[aa][ab] = ndcut
      
        
        return dat_dic
    
        
    def calc_ndSmid_cut(self):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        
        if withshift == False and withseries == False:
            
            nx = self.data['b2fgeo']['nx']
            ny = self.data['b2fgeo']['ny']
            dat_struc = {'nx': nx, 'ny': ny}
            cut_list = [0.95, 1.1]
                
            ndSmidcut_dic = self.ndSmid_withcut_method(itername = None, data_struc = dat_struc, cut_range = cut_list)
            
            self.data['ndSmid_cutprofile'] = ndSmidcut_dic
        
        elif withshift == True and withseries == False:
            
            
            ndSmidcut_dic = {}
            
            
            for aa in self.data['dircomp']['multi_shift']:
                
                
                nx = self.data['b2fgeo'][aa]['nx']
                ny = self.data['b2fgeo'][aa]['ny']
                dat_struc = {'nx': nx, 'ny': ny}
                cut_list = [0.95, 1.1]
                
                
                ndSmidcut = self.ndSmid_withcut_method(itername = aa, data_struc = dat_struc, cut_range = cut_list)
                
                ndSmidcut_dic[aa] = ndSmidcut
            
            self.data['ndSmid_cutprofile'] = ndSmidcut_dic
        
        elif withshift == False and withseries == True:
            
            midprofile_dic = {}
            
            nx = self.data['b2fgeo']['nx']
            ny = self.data['b2fgeo']['ny']
            dat_struc = {'nx': nx, 'ny': ny}
            cut_list = [0.95, 1.1]
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            if self.DF.series_flag == 'twin_scan':
                
                ds_key, ts_key = self.lbd.twokeylists(printvalue= False)
                
                ndcut_dic = self.ldm.two_layer_dic(key_a = ds_key, key_b = ts_key)
                
                ndSmidcut_dic = self.twoDscan_ndScut(iterlist = scan, iterlist_a = ds_key, 
                                    iterlist_b = ts_key, dat_dic = ndcut_dic, dat_struc = dat_struc, 
                                    cut_range = cut_list)
                
            else:
                
                ndc_dic = {}
                ndSmidcut_dic = self.oneDscan_ndScut(iterlist = scan, dat_dic = ndc_dic, 
                                                   dat_struc = dat_struc, cut_range = cut_list)
            
            self.data['ndSmid_cutprofile'] = ndSmidcut_dic
        
        elif withshift == True and withseries == True:
            print('calc_midplane_profile is not there yet!')
        
        
        else:
            print('calc_midplane_profile has a bug')
    


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
        
        mid_neu_pro = np.multiply(neu_pro[:, 58], weight) + np.multiply(neu_pro[:, 60], weight_B)
        mid_atom_temp_pro = np.multiply(atom_temp_pro[:, 58], weight) + np.multiply(atom_temp_pro[:, 60], weight_B)
        
        return psi_coord, mid_neu_pro, mid_atom_temp_pro




"""

    def twinscan_Smid(self, iter_index, data_struc):
        
        result_dic = self.data['radial_fit_data'][iter_index]
        
        
        if self.series_flag == 'twin_scan':
            
            nf = iter_index[0]
            tf = iter_index[1]
            
            nx = data_struc['nx']
            ny = data_struc['ny']
            source = self.data['iout_data']['source'][nf][tf][:, :]
            weight = self.data['midplane_calc']['weight'][1:ny+1]
            psi_coord = self.data['midplane_calc']['psi_solps_mid'][1:ny+1]

            
        else:
            
            print('twinscan_Smid is not there yet!')
            

        weight_B = np.ones(len(weight))- weight
               
        mid_S_pro = np.multiply(source[:, 58], weight) + np.multiply(source[:, 60], weight_B)
        
        psi_list = []
        S_list = []
        
        for ind, coord in enumerate(psi_coord):
            
            if coord >= 0.95 and coord <= 1.1:
                psi_list.append(coord)
                S_list.append(mid_S_pro[ind])
        
        
        
        return psi_list, S_list


"""






