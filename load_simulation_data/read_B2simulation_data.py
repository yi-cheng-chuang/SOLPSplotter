# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:11:08 2024

@author: ychuang
"""

 
from load_directory.load_dirdata_method import load_dir_method
from load_simulation_data.load_B2_data_method import read_b2fstate, read_b2fplasmf, read_b2wdat, read_iout_method
import numpy as np



class load_B2simu_data:
    
    def __init__(self, DF, data, ldm: load_dir_method):
        
        self.DF = DF
        self.data = data
        self.ldm = ldm
  
    
    def twokeylists(self, printvalue):
        
        mcds = self.data['dircomp']
        
        ds_key = []
        ts_key = []
        
        if printvalue:
            print('this is mcds')
            print(type(mcds['denscan_list'][3]))
        else:
            pass
            
        
        for x in mcds['denscan_list']:
            ds_key.append('{:.3f}'.format(x))
        
        if printvalue:
            print(ds_key)
        else:
            pass
            
        
        for x in mcds['tempscan_list']:
            ts_key.append('{:.3f}'.format(x))
            
        if printvalue:
            print(ts_key)
        else:
            pass
        
        return ds_key, ts_key
        
    
    
    def oneDscan_dat(self, iterlist, dat_type):
        
        
        series_compare = self.DF.series_compare
        
        if series_compare == True:
            
            if dat_type == 'b2fstate':
                
                state_dic = {}
                dim_dic = {}
                
                for aa in iterlist:
                    
                    state_cpdic = {'fixed': {}, 'flux': {}}
                    for kk in ['fixed', 'flux']:
                        file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa][kk], 'b2fstate')
                        state, dim = read_b2fstate(b2fstateLoc = file_loc)
                        state_cpdic[kk] = vars(state)
                    
                    
                    state_dic[aa] = state_cpdic
                    dim_dic[aa] = {'nx': dim[0], 'ny': dim[1], 'ns': dim[2]}
                
                return state_dic, dim_dic
            
            
            elif dat_type == 'b2wdat':
                
                print('we will improve this in the future!')
                
        
        else:
            
            if dat_type == 'b2fstate':
                
                dat_dic = {}
                dim_dic = {}
                
                for aa in iterlist:
                    file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], 'b2fstate')
                    state, dim = read_b2fstate(b2fstateLoc = file_loc)
                    dat_dic[aa] = vars(state)
                    dim_dic[aa] = {'nx': dim[0], 'ny': dim[1], 'ns': dim[2]}
                    
                
                return dat_dic, dim_dic
            
            
            elif dat_type == 'b2wdat':
                
                b2wdat_dic = {}
                
                for aa in iterlist:
                                   
                    file_loc = '{}/'.format(self.data['dirdata']['simudir'][aa])
                    na_dat = self.data['b2fstate'][aa]['na']
                    b2wdat = read_b2wdat(b2wdatLoc = file_loc, 
                                              nSpec = np.shape(na_dat)[2])
                    b2wdat_dic[aa] = vars(b2wdat)


                return b2wdat_dic
                
  
    
    def twoDscan_dat(self, iterlist, iterlist_a, iterlist_b, dat_type, dat_dic, dim_dic):
        
        
        if dat_type == 'b2fstate':
            
            for tp in iterlist:
                aa = tp[0]
                ab = tp[1]
                
                file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa][ab], 'b2fstate')
                state, dim = read_b2fstate(b2fstateLoc = file_loc)
                dat_dic[aa][ab] = vars(state)
                dim_dic[aa][ab] = {'nx': dim[0], 'ny': dim[1], 'ns': dim[2]}
            
            
            return dat_dic, dim_dic
        
        elif dat_type == 'b2wdat':
            
            
            for tp in iterlist:
            
                aa = tp[0]
                ab = tp[1]
                
                file_loc = '{}/'.format(self.data['dirdata']['simudir'][aa][ab])
                na_dat = self.data['b2fstate'][aa][ab]['na']
                           
                b2wdat = read_b2wdat(b2wdatLoc = file_loc, 
                                          nSpec = np.shape(na_dat)[2])
                dat_dic[aa][ab] = vars(b2wdat)
            
            return dat_dic


    
    
    def load_b2fstate(self):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        
        if withshift == False and withseries == False:
            
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'], 'b2fstate')
            state, dim = read_b2fstate(b2fstateLoc = file_loc)
            state_dic = vars(state)
            dim_dic = {'nx': dim[0], 'ny': dim[1], 'ns': dim[2]}
            self.data['b2fstate'] = state_dic
            self.data['DefaultSettings']['dims'] = dim_dic
            # self.b2fstate = state
        
        elif withshift == True and withseries == False:
            
            scan = self.data['dircomp']['multi_shift']
            
            state_dic, dim_dic = self.oneDscan_dat(iterlist = scan, dat_type = 'b2fstate')

            self.data['b2fstate'] = state_dic
            self.data['DefaultSettings']['dims'] = dim_dic
        
        elif withshift == False and withseries == True:
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            series_flag = self.DF.series_flag
            
            
            if series_flag == 'twin_scan':
                
                
                ds_key, ts_key = self.twokeylists(printvalue= False)
                
                state_dic = self.ldm.two_layer_dic(key_a = ds_key, key_b = ts_key)
                # print(state_dic)
                dim_dic = self.ldm.two_layer_dic(key_a = ds_key, key_b = ts_key)
                
                state_dic, dim_dic = self.twoDscan_dat(iterlist = scan, 
                                    iterlist_a = ds_key, iterlist_b = ts_key, dat_type = 'b2fstate',
                                    dat_dic = state_dic, dim_dic = dim_dic)
            else:
                state_dic, dim_dic = self.oneDscan_dat(iterlist = scan, dat_type = 'b2fstate')
                

            self.data['b2fstate'] = state_dic
            self.data['DefaultSettings']['dims'] = dim_dic
            
        
        else:
            print('load_b2fstate function is not there yet!')
    

       
    
    def load_b2wdat(self):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        if withshift == False and withseries == False:
            
            file_loc = '{}/'.format(self.data['dirdata']['simudir'])
            na_dat = self.data['b2fstate']['na']
            
            
            b2wdat = read_b2wdat(b2wdatLoc = file_loc, 
                                      nSpec = np.shape(na_dat)[2])
            b2wdat_dic = vars(b2wdat)
                
                   

            self.data['b2wdat'] = b2wdat_dic
        
        
        
        elif withshift == True and withseries == False:
            
            b2wdat_dic = {}

            for aa in self.data['dircomp']['multi_shift']:
                
                file_loc = '{}/'.format(self.data['dirdata']['simudir'][aa])
                na_dat = self.data['b2fstate']['org']['na']
                
                
                b2wdat = read_b2wdat(b2wdatLoc = file_loc, 
                                          nSpec = np.shape(na_dat)[2])
                b2wdat_dic[aa] = vars(b2wdat)   

            self.data['b2wdat'] = b2wdat_dic
        
        
        
        elif withshift == False and withseries == True:
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            if self.DF.series_flag == 'twin_scan':
                
                ds_key, ts_key = self.twokeylists(printvalue= False)
                
                wdat_dic = self.ldm.two_layer_dic(key_a = ds_key, key_b = ts_key)
                dim_dic = self.ldm.two_layer_dic(key_a = ds_key, key_b = ts_key)
                
                b2wdat_dic = self.twoDscan_dat(iterlist = scan, 
                                    iterlist_a = ds_key, iterlist_b = ts_key, dat_type = 'b2wdat',
                                    dat_dic = wdat_dic, dim_dic = dim_dic)
            else:
                b2wdat_dic = self.oneDscan_dat(iterlist = scan, dat_type = 'b2wdat')
                

            self.data['b2wdat'] = b2wdat_dic
            # self.data['DefaultSettings']['dims'] = dim_dic
            
        
        else:
            print('load_b2fstate function is not there yet!')
        
    
    
    
    
    def load_b2fplasmf(self):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        
        if withshift == False and withseries == False:
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'], 'b2fplasmf')
            dim_dic = self.data['DefaultSettings']['dims']
            n_pol = dim_dic['nx']
            n_rad = dim_dic['ny']
            n_sp = dim_dic['ns']
            fplasma = read_b2fplasmf(fileName = file_loc, nx = n_pol, 
                                          ny = n_rad, ns = n_sp)
            fplasma_dic = vars(fplasma)
            
            self.data['b2fplasmf'] = fplasma_dic
            # print('the next line is b2fplasmf')
            # print(type(k))
        
        elif withshift == True and withseries == False:
            fplasma_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], 'b2fplasmf')
                dim_dic = self.data['DefaultSettings']['dims']
                n_pol = dim_dic[aa]['nx']
                n_rad = dim_dic[aa]['ny']
                n_sp = dim_dic[aa]['ns']
                fplasma = read_b2fplasmf(fileName = file_loc, nx = n_pol, 
                                              ny = n_rad, ns = n_sp)
                fplasma_dic[aa] = vars(fplasma)
                
                
            self.data['b2fplasmf'] = fplasma_dic
        
        elif withshift == False and withseries == True:
            fplasma_dic = {}
            
            for aa in list(self.data['dircomp']['Attempt'].keys()):
                file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], 'b2fplasmf')
                dim_dic = self.data['DefaultSettings']['dims']
                n_pol = dim_dic[aa]['nx']
                n_rad = dim_dic[aa]['ny']
                n_sp = dim_dic[aa]['ns']
                fplasma = read_b2fplasmf(fileName = file_loc, nx = n_pol, 
                                              ny = n_rad, ns = n_sp)
                fplasma_dic[aa] = vars(fplasma)
                
                
            self.data['b2fplasmf'] = fplasma_dic
            
                
        else:
            print('load_b2fplasmf function is not there yet!')
        
        
    
    
    def twinscan_iout(self, iterlist, iterlist_a, iterlist_b, filename):
        

        iout_dic = self.ldm.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
        # print(iout_dic)
        ny = self.data['b2fgeo']['ny']
        nx = self.data['b2fgeo']['nx']
        
        for tp in iterlist:
            aa = tp[0]
            ab = tp[1]
            
            file_loc = '{}/{}/{}'.format(self.data['dirdata']['simudir'][aa][ab], 'output', filename)
            iout = read_iout_method(fdir = file_loc, fname = filename, 
                                  ny = ny, nx = nx)
            
            iout_dic[aa][ab] = iout
        return iout_dic
    
    
    
    def load_iout(self, filename, simple_quant):
        filename_list = filename.split('.')
        quant = simple_quant
        
        # if simple_quant:
        #     pass
        # else:
        #     quant = quant.split('_')[1]
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        series_flag = self.DF.series_flag
        
        
        if withshift == False and withseries == False:
            file_loc = '{}/{}/{}'.format(self.data['dirdata']['simudir'], 'output', filename)           
            ny = self.data['b2fgeo']['ny']
            nx = self.data['b2fgeo']['nx']
            iout = read_iout_method(fdir = file_loc, fname = filename, 
                                  ny = ny, nx = nx)
            
            self.data['iout_data'][quant] = iout

        
        elif withshift == True and withseries == False:
            iout_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                ny = self.data['b2fgeo'][aa]['ny']
                nx = self.data['b2fgeo'][aa]['nx']
                file_loc = '{}/{}/{}'.format(self.data['dirdata']['simudir'][aa], 'output', filename)           
                iout = read_iout_method(fdir = file_loc, fname = filename, 
                                      ny = ny, nx = nx)
                
                iout_dic[aa] = iout
            
            self.data['iout_data'][quant] = iout_dic
                
                
        
        elif withshift == False and withseries == True:
            
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            
            if series_flag == 'twin_scan':
                
                ds_key, ts_key = self.twokeylists(printvalue= False)
                
                
                iout_dic = self.twinscan_iout(iterlist = scan, iterlist_a = ds_key, 
                              iterlist_b = ts_key, filename = filename)
            else:
                
                iout_dic = {}
                for aa in list(self.data['dircomp']['Attempt'].keys()):
                    ny = self.data['b2fgeo']['ny']
                    nx = self.data['b2fgeo']['nx']
                    
                    file_loc = '{}/{}/{}'.format(self.data['dirdata']['simudir'][aa], 'output', filename)           
                    iout = read_iout_method(fdir = file_loc, fname = filename, 
                                          ny = ny, nx = nx)
                    
                    iout_dic[aa] = iout
                
            
            
            self.data['iout_data'][quant] = iout_dic
            
            
        
        else:
            print('load_b2fstate function is not there yet!')
        
        
        return quant
    
    def load_iout_ratio(self, file_tuple, itername):
        
        
        series_flag = self.DF.series_flag
        
        if len(file_tuple) > 2:
            print('input more than two files!')
        
        else:
            for ftu in file_tuple:
                self.load_iout(filename = ftu[0], simple_quant = ftu[1])
            
            if itername == None:
                
                name1 = file_tuple[0][1]
                data1 = self.data['iout_data'][name1]
                
                name2 = file_tuple[1][1]
                data2 = self.data['iout_data'][name2]
                
                ratio_data = np.divide(data1, data2)
                quantname = '{}_divide_{}'.format(name1, name2)
                
                self.data['iout_data'][quantname] = ratio_data
                # self.data['iout_data']['{}_abs'.format(quantname)] = np.abs(ratio_data)
            
            elif series_flag == 'twin_scan':
                
                
                name1 = file_tuple[0][1]
                nf = itername[0]
                tf = itername[1]
                
                print(itername)
                
                data1 = self.data['iout_data'][name1][nf][tf]
                
                name2 = file_tuple[1][1]
                data2 = self.data['iout_data'][name2][nf][tf]
                
                ratio_data = np.divide(data1, data2)
                quantname = '{}_divide_{}'.format(name1, name2)
            
            
            else:
                
                name1 = file_tuple[0][1]
                data1 = self.data['iout_data'][name1][itername]
                
                name2 = file_tuple[1][1]
                data2 = self.data['iout_data'][name2][itername]
                
                ratio_data = np.divide(data1, data2)
                quantname = '{}_divide_{}'.format(name1, name2)
                
                # self.data['iout_data'][quantname][itername] = ratio_data
                # self.data['iout_data']['{}_abs'.format(quantname)] = np.abs(ratio_data)
                

            return ratio_data, quantname
        
            
            
    
        
    def load_iout_multi(self, name1, name2, itername, input_name):
        
        if input_name == None and itername == None:
            
            data1 = self.data['iout_data'][name1]
            
            data2 = self.data['iout_data'][name2]
            
            multi_data = np.multiply(data1, data2)
            quantname = '{}_multi_{}'.format(name1, name2)
            
            self.data['iout_data'][quantname] = multi_data
            # self.data['iout_data']['{}_abs'.format(quantname)] = np.abs(ratio_data)
        
        elif input_name != None and itername == None:
            
            data1 = self.data['iout_data'][name1]
            
            data2 = self.data['iout_data'][name2]
            
            multi_data = np.multiply(data1, data2)
            quantname = '{}'.format(input_name)
            
        
        elif input_name != None and itername != None:
            data1 = self.data['iout_data'][name1][itername]
            
            data2 = self.data['iout_data'][name2][itername]
            
            multi_data = np.multiply(data1, data2)
            quantname = '{}'.format(input_name)
        
        return multi_data, quantname
    
    
    def load_iout_divide(self, name1, name2, itername, input_name):
        
        if input_name == None and itername == None:
            
            data1 = self.data['iout_data'][name1]
            
            data2 = self.data['iout_data'][name2]
            
            multi_data = np.divide(data1, data2)
            quantname = '{}_divide_{}'.format(name1, name2)
            
            self.data['iout_data'][quantname] = multi_data
            # self.data['iout_data']['{}_abs'.format(quantname)] = np.abs(ratio_data)
        
        elif input_name != None and itername == None:
            data1 = self.data['iout_data'][name1]
            
            data2 = self.data['iout_data'][name2]
            
            multi_data = np.divide(data1, data2)
            quantname = '{}'.format(input_name)
        
            
        elif input_name != None and itername != None:
            data1 = self.data['iout_data'][name1][itername]
            
            data2 = self.data['iout_data'][name2][itername]
            
            multi_data = np.divide(data1, data2)
            quantname = '{}'.format(input_name)
        
        return multi_data, quantname
    
    
    
    
    def load_iout_name_ratio(self, setname, name1, name2, stdname, itername):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        if withshift == False and withseries == False:
        
            data1 = self.data['iout_data'][name1]
            
            data2 = self.data['iout_data'][name2]
            
            ratio_data = np.divide(data1, data2)
            quantname = '{}_divide_{}'.format(name1, name2)
            
            self.data['iout_data'][quantname] = ratio_data
            # self.data['iout_data']['{}_abs'.format(quantname)] = np.abs(ratio_data)
            
            return quantname
        
        elif withshift == True and withseries == False:
            
            if setname == 'iout_data':
                
                data1 = self.data[setname][name1]
                
                data2 = self.data[setname][name2][stdname]
                
                ratio_data = np.divide(data1, data2)*100
                quantname = '{}_change_percent'.format(name2)
                
                # self.data['iout_data'][quantname] = ratio_data
                # self.data['iout_data']['{}_abs'.format(quantname)] = np.abs(ratio_data)
                
                return ratio_data, quantname
            
            elif setname == 'ft44':
                
                data1 = self.data[setname][name1]
                
                data2 = self.data[setname][stdname][name2]
                
                ratio_data = np.divide(data1, data2)*100
                quantname = '{}_change_percent'.format(name2)
                
                # self.data['iout_data'][quantname] = ratio_data
                # self.data['iout_data']['{}_abs'.format(quantname)] = np.abs(ratio_data)
                
                return ratio_data, quantname
                
            
           
        
        
        else:
            print('load_iout_name_ratio function is not there yet!')
        
    
        
        
    
    
    def load_differ(self, setname, name, itername1, itername2):
        
        if setname == 'iout_data':
            data1 = self.data[setname][name][itername1]
            
            data2 = self.data[setname][name][itername2]
            
            diff_data = data1 - data2
            quantname = '{}_{}_minus_{}'.format(name, itername1, itername2)
            
            self.data[setname][quantname] = diff_data
            # self.data['iout_data']['{}_abs'.format(quantname)] = np.abs(ratio_data)
        
        elif setname == 'ft44':
            
            data1 = self.data[setname][itername1][name]
            
            data2 = self.data[setname][itername2][name]
            
            diff_data = data1 - data2
            quantname = '{}_{}_minus_{}'.format(name, itername1, itername2)
            
            self.data[setname][quantname] = diff_data
            # self.data['iout_data']['{}_abs'.format(quantname)] = np.abs(ratio_data)
            
            
        
        
        return quantname
            
            
            
"""

state_dic = {}
dim_dic = {}

for aa in self.data['dircomp']['multi_shift']:
    file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], 'b2fstate')
    state, dim = lbdm.read_b2fstate(b2fstateLoc = file_loc)
    state_dic[aa] = vars(state)
    dim_dic[aa] = {'nx': dim[0], 'ny': dim[1], 'ns': dim[2]}



state_dic = {}
dim_dic = {}

for aa in list(self.data['dircomp']['Attempt'].keys()):
    file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], 'b2fstate')
    state, dim = lbdm.read_b2fstate(b2fstateLoc = file_loc)
    state_dic[aa] = vars(state)
    dim_dic[aa] = {'nx': dim[0], 'ny': dim[1], 'ns': dim[2]}

self.data['b2fstate'] = state_dic
self.data['DefaultSettings']['dims'] = dim_dic
# self.b2fstate = state


for aa in self.data['dircomp']['multi_shift']:
    
    file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], '{}'.format(ftname))
    ft44 = lEdm.read_ft44(fileName = file_loc)
    ft44_dic[aa] = vars(ft44)



"""     
        
            
        
            
            
            

            
            
    
    
    
                
                
            

                    
                
                
            
        
    
    
    