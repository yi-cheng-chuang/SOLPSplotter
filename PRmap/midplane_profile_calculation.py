# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 16:12:28 2025

@author: ychuang
"""

 
import numpy as np




class midplane_radial:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    def calc_midplane_profile_method(self, itername, data_struc):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        if withshift == False and withseries == False:
            
            b2fstate = self.data['b2fstate']
                
            if self.DF.data_size == 'small':
                
                data = self.data['ft44']['dab2']
                neu_pro = np.transpose(data[:, :, 0])
            
            else:
                
                print('I do not use full data size for my study')
            
            

        elif withshift == True and withseries == False:
            
            
            b2fstate = self.data['b2fstate'][itername]
                
            if self.DF.data_size == 'small':
                
                data = self.data['ft44'][itername]['dab2']
                neu_pro = np.transpose(data[:, :, 0])
            
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
                
                else:
                    
                    print('I do not use full data size for my study')
                
            else:
                
                b2fstate = self.data['b2fstate'][itername]
                    
                if self.DF.data_size == 'small':
                    
                    data = self.data['ft44'][itername]['dab2']
                    neu_pro = np.transpose(data[:, :, 0])
                
                else:
                    
                    print('I do not use full data size for my study')
                
                
        
        nx = data_struc['nx']
        ny = data_struc['ny']
        
            
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
        
        
        midplane_profile_dic = {'psiN': psi_coord, 'mid_ne': mid_ne_pro, 'mid_te': mid_te_pro,
                                'mid_nd': mid_neu_pro}
        
        
        return midplane_profile_dic
    
    
    def oneDscan_midplan_profile(self):
        
        series_flag = self.DF.series_flag
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        dat_struc = {'nx': nx, 'ny': ny}
        
        midprofiles = self.calc_midplane_profile_method(itername = aa, data_struc = dat_struc)
        
        midprofile_dic[aa] = midprofiles
    
    
    def twoDscan_b2wdat(self, iterlist, iterlist_a, iterlist_b, two_layer_dic):
        


        b2wdat_dic = self.ldm.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
        # print(state_dic)
        # dim_dic = self.ldm.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
        
        for tp in iterlist:
            aa = tp[0]
            ab = tp[1]
            
            file_loc = '{}/'.format(self.data['dirdata']['simudir'][aa][ab])
            na_dat = self.data['b2fstate'][aa][ab]['na']
                       
            b2wdat = read_b2wdat(b2wdatLoc = file_loc, 
                                      nSpec = np.shape(na_dat)[2])
            b2wdat_dic[aa][ab] = vars(b2wdat)
            # dim_dic[aa][ab] = {'nx': dim[0], 'ny': dim[1], 'ns': dim[2]}
        
        return b2wdat_dic
    
    
    
    
    
    def calc_midplane_profile(self):
        
        # self.load_ft44()
        
        dat_size = self.DF.data_size
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
            
            
            for aa in list(self.data['dircomp']['Attempt'].keys()):
                
                
                
                
                
                
                scan = list(self.data['dircomp']['Attempt'].keys())
                
                if self.DF.series_flag == 'twin_scan':
                    
                    mcds = self.data['dircomp']
                    
                    ds_key = []
                    ts_key = []
                    
                    print('this is mcds')
                    print(type(mcds['denscan_list'][3]))
                    
                    for x in mcds['denscan_list']:
                        ds_key.append('{:.3f}'.format(x))
                        
                    print(ds_key)
                    for x in mcds['tempscan_list']:
                        ts_key.append('{:.3f}'.format(x))
                        

                    print(ts_key)
                    
                    
                    b2wdat_dic = self.ldm.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
                    # print(state_dic)
                    # dim_dic = self.ldm.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
                    
                    for tp in iterlist:
                        aa = tp[0]
                        ab = tp[1]
                        
                        file_loc = '{}/'.format(self.data['dirdata']['simudir'][aa][ab])
                        na_dat = self.data['b2fstate'][aa][ab]['na']
                                   
                        b2wdat = read_b2wdat(b2wdatLoc = file_loc, 
                                                  nSpec = np.shape(na_dat)[2])
                        b2wdat_dic[aa][ab] = vars(b2wdat)
                        # dim_dic[aa][ab] = {'nx': dim[0], 'ny': dim[1], 'ns': dim[2]}

                b2wdat_dic = self.twoDscan_b2wdat(iterlist = scan, 
                                    iterlist_a = ds_key, iterlist_b = ts_key)
                
                
                
                else:
                    b2wdat_dic = self.oneDscan_b2wdat(iterlist = scan)
                
            
            self.data['midplane_profile'] = midprofile_dic
        
        elif withshift == True and withseries == True:
            print('calc_midplane_profile is not there yet!')
        
        
        else:
            print('calc_midplane_profile has a bug')
    
    
    
    def ndmid_withcut_method(self, itername, data_struc):

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
            
            else:
                
                print('I do not use full data size for my study')

        elif withshift == True and withseries == False:
            
            if data_struc['size'] == 'small':
                nx = data_struc['nx']
                ny = data_struc['ny']
                data = self.data['ft44'][itername]['dab2']
                neu_pro = np.transpose(data[:, :, 0])
                weight = self.data['midplane_calc'][itername]['weight'][1:ny+1]
                psi_coord = self.data['midplane_calc'][itername]['psi_solps_mid'][1:ny+1]
            
            else:
                
                print('I do not use full data size for my study')
            
        
        elif withshift == False and withseries == True:
            
            if self.series_flag == 'twin_scan':
                
                nf = itername[0]
                tf = itername[1]
                    
                                  
                if data_struc['size'] == 'small':
                    nx = data_struc['nx']
                    ny = data_struc['ny']
                    data = self.data['ft44'][nf][tf]['dab2']
                    neu_pro = np.transpose(data[:, :, 0])
                    weight = self.data['midplane_calc']['weight'][1:ny+1]
                    psi_coord = self.data['midplane_calc']['psi_solps_mid'][1:ny+1]
                
                else:
                    
                    print('I do not use full data size for my study')
                
            else:
                
                if data_struc['size'] == 'small':
                    nx = data_struc['nx']
                    ny = data_struc['ny']
                    data = self.data['ft44'][itername]['dab2']
                    neu_pro = np.transpose(data[:, :, 0])
                    weight = self.data['midplane_calc']['weight'][1:ny+1]
                    psi_coord = self.data['midplane_calc']['psi_solps_mid'][1:ny+1]
                
                else:
                    
                    print('I do not use full data size for my study')
            

        weight_B = np.ones(len(weight))- weight
            
        
        if data_struc['size'] == 'small':
            
            ap = self.data['midplane_calc']['average_pair']
            pair = (ap[0]-1, ap[1]-1)
        
        else:
            
            print('I do not use full data size for my study')
        
        mid_neu_pro = np.multiply(neu_pro[:, pair[0]], weight) + np.multiply(neu_pro[:, pair[1]], weight_B)
        
        psi_list = []
        nd_list = []
        
        for ind, coord in enumerate(psi_coord):
            
            if coord >= 0.95 and coord <= 1.1:
                psi_list.append(coord)
                nd_list.append(mid_neu_pro[ind])
        
        
        
        return psi_list, nd_list
        
    
    
    
    
    
    
    
    
    
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
