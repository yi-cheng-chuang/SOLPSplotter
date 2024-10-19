# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:11:08 2024

@author: ychuang
"""


from SOLPSplotter_load_expdata import load_expdata
import load_mast_expdata_method as lmem
import load_B2_data_method as lbdm
import load_Eirene_data_method as lEdm
import numpy as np



class load_simu_data(load_expdata):
    
    def __init__(self, DefaultSettings, loadDS):
        load_expdata.__init__(self, DefaultSettings, loadDS)
        # Employee.__init__(self, first, last, pay)
        
        "Parameters"
        if isinstance(DefaultSettings['Parameters'], dict):
            self.Parameters = DefaultSettings['Parameters']
        else:
            print('parameter has to be a dictionary')
            
        if DefaultSettings['Parameters'] is None:
            print('There is no parameters input')
        else:
            self.Parameters = DefaultSettings['Parameters']
            
        Plist = []
        for pkey, pvalue in self.Parameters.items():
            Plist.append(pkey)
        
        self.data['paramkey'] = Plist
        self.data['Parameter'] = DefaultSettings['Parameters']

        
    "Add and remove elements from parameter or defaultsettings"
    def add_dic(self, new_set, assign='default'):
        if assign == 'param':
            self.Parameters = new_set | self.Parameters
            Plist = []
            for pkey, pvalue in self.Parameters.items():
                Plist.append(pkey)
        
        # elif assign == 'default':
        #     self.DefaultSettings = new_set | self.DefaultSettings
        #     keylist = []
        #     for key, value in self.DefaultSettings.items():
        #         keylist.append(key)
        else:
            print('assign parameter is incorrect')
            
        
    def remove_dic(self, new_set, assign='param'):
        if assign == 'param':
            if new_set.keys() in self.data['defaultkey']:
                del self.Parameters[new_set.keys()]
            Plist = []
            for pkey, pvalue in self.Parameters.items():
                Plist.append(pkey)
                
        # elif assign == 'default':
        #     if new_set.keys() in self.data['defaultkey']:
        #         del self.DefaultSettings[new_set.keys()]
        #     keylist = []
        #     for key, value in self.DefaultSettings.items():
        #         keylist.append(key)
        else:
            print('assign parameter incorrect')

    
    def load_output_data_method(self, param, itername):
        if self.withshift == False and self.withseries == False:
            BASEDRT = self.data['dirdata']['outputdir']['Output']
            Attempt = self.data['dircomp']['Attempt']
            XGrid = int(self.data['b2fgeo']['nx'])
            # X = self.data['gridsettings']['X']
            # YSurf = int(self.data['b2fgeo']['ny'])
            # Y = self.data['gridsettings']['Y']
            XMin= 1
            XMax= XGrid
            # print(XGrid)
            XDIM = int(self.data['DefaultSettings']['XDIM'])
            YDIM = int(self.data['DefaultSettings']['YDIM'])
        elif self.withshift == True and self.withseries == False:
            BASEDRT = self.data['dirdata']['outputdir'][itername]['Output']
            Attempt = self.data['dircomp']['Attempt'][itername]
            XGrid = int(self.data['b2fgeo'][itername]['nx'])
            # print(XGrid)
            XDIM = int(self.data['DefaultSettings']['XDIM'][itername])
            YDIM = int(self.data['DefaultSettings']['YDIM'][itername])
            
        elif self.withshift == False and self.withseries == True:
            
            if self.series_flag == 'twin_scan':
                
                nf = itername[0]
                tf = itername[1]
                
                simu_dir = self.data['dirdata']['simudir'][nf][tf]
                
                BASEDRT = '{}/{}'.format(simu_dir, 'Output')
            else:
                BASEDRT = self.data['dirdata']['outputdir'][itername]['Output']

            Attempt = self.data['dircomp']['Attempt'][itername]
            XGrid = int(self.data['b2fgeo']['nx'])
            # print(XGrid)
            XDIM = int(self.data['DefaultSettings']['XDIM'])
            YDIM = int(self.data['DefaultSettings']['YDIM'])
        
        # DRT = '{}/Attempt{}'.format(BASEDRT, str(Attempt))   #Generate path
        # self.data['outputdata'][param] = xr.DataArray(np.zeros((YSurf,XGrid,N)), 
        #                                  coords=[Y,X,Attempts], 
        # dims=['Radial_Location','Poloidal_Location','Attempt'], name = param)
        
        n = 0
        output_data = np.zeros([YDIM, XDIM], dtype= np.float32)
        test = param in self.Parameters.keys()       
        if test:
            # print('yes, {} is in parameter'.format(param))
            RawData = np.loadtxt('{}/{}{}'.format(BASEDRT, param, str(Attempt)),usecols = (3))
            # temp_dic = {param: RawData}
            # self.data['temp'] = temp_dic
        elif test == False:
            print('no, {} is not in parameter'.format(param))
        else:
            print('there might be a bug')
            
            
        if len(RawData) > 0:        
            if RawData.size == XDIM*YDIM:
                # self.data['outputdata'][param].values[:,:,n] = RawData.reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
                output_data = RawData.reshape((YDIM,XDIM))
            elif RawData.size == XDIM*YDIM*2:
                raw_split = np.array_split(RawData, 2)
                param_dic = {'D_0': raw_split[0].reshape((YDIM,XDIM)), 
                             'D_1': raw_split[1].reshape((YDIM,XDIM))}
                output_data = param_dic
                # print('let work on it')
                
                
            elif RawData.size != XDIM*YDIM:
                print('rawdata size is not equal to {}, it is {}'.format(str(XDIM*YDIM), str(RawData.size)))
            # elif RawData.size == XDIM*YDIM*2:
            #     self.data['outputdata'][param].values[:,:,n] = RawData.reshape((2*YDIM,XDIM))[1+YDIM:2*YDIM-1,XMin:XMax+1]
        else:
            print('we have a problem loading rawdata')
        
        return output_data
    
    def one_dim_scan_output(self, iterlist, param):
        
        param_data_dic = {}
        for aa in iterlist:
            param_data_dic[aa] = self.load_output_data_method(param = param, itername = aa)
        
        return param_data_dic
    
    
    def two_dim_scan_output(self, iterlist, iterlist_a, iterlist_b, param):
        


        param_data_dic = lmem.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
        
        for tp in iterlist:
            aa = tp[0]
            ab = tp[1]
            
            param_data_dic[aa][ab] = self.load_output_data_method(param = param, itername = tp)
        
        return param_data_dic
       
    def load_output_data(self, param):
        if self.withshift == False and self.withseries == False:
            output = self.load_output_data_method(param = param, itername = None)
            self.data['outputdata'][param] = output
            
            
        elif self.withshift == True and self.withseries == False:
            
            
            scan = self.data['dircomp']['multi_shift']
            
            param_data_dic = self.one_dim_scan_output(iterlist = scan, param = param)
                      
            self.data['outputdata'][param] = param_data_dic
            
            
        elif self.withshift == False and self.withseries == True:
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            if self.series_flag == 'twin_scan':
                
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
                
                
                param_data_dic = self.two_dim_scan_output(iterlist = scan, iterlist_a = ds_key, 
                                    iterlist_b = ts_key, param = param)
            else:
                param_data_dic = self.one_dim_scan_output(iterlist = scan, param = param)
        
            self.data['outputdata'][param] = param_data_dic
        
        elif self.withshift == True and self.withseries == True:
            print('load_output_data is not there yet, to be continue...')
        
        else:
            print('There is a bug')
    
    
    
    def one_dim_scan_b2fstate(self, iterlist):
        
        if self.series_compare == True:
            
            state_dic = {}
            dim_dic = {}
            
            for aa in iterlist:
                
                state_cpdic = {'fixed': {}, 'flux': {}}
                for kk in ['fixed', 'flux']:
                    file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa][kk], 'b2fstate')
                    state, dim = lbdm.read_b2fstate(b2fstateLoc = file_loc)
                    state_cpdic[kk] = vars(state)
                
                
                state_dic[aa] = state_cpdic
                dim_dic[aa] = {'nx': dim[0], 'ny': dim[1], 'ns': dim[2]}
        
        else:
            
            state_dic = {}
            dim_dic = {}
            
            for aa in iterlist:
                file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], 'b2fstate')
                state, dim = lbdm.read_b2fstate(b2fstateLoc = file_loc)
                state_dic[aa] = vars(state)
                dim_dic[aa] = {'nx': dim[0], 'ny': dim[1], 'ns': dim[2]}
        
        
        
        
        return state_dic, dim_dic
    
    

    
    
    
    
    def two_dim_scan_b2fstate(self, iterlist, iterlist_a, iterlist_b):
        


        state_dic = lmem.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
        print(state_dic)
        dim_dic = lmem.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
        
        for tp in iterlist:
            aa = tp[0]
            ab = tp[1]
            
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa][ab], 'b2fstate')
            state, dim = lbdm.read_b2fstate(b2fstateLoc = file_loc)
            state_dic[aa][ab] = vars(state)
            dim_dic[aa][ab] = {'nx': dim[0], 'ny': dim[1], 'ns': dim[2]}
        
        return state_dic, dim_dic
    
    
    
    def load_b2fstate(self):
        if self.withshift == False and self.withseries == False:
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'], 'b2fstate')
            state, dim = lbdm.read_b2fstate(b2fstateLoc = file_loc)
            state_dic = vars(state)
            dim_dic = {'nx': dim[0], 'ny': dim[1], 'ns': dim[2]}
            self.data['b2fstate'] = state_dic
            self.data['DefaultSettings']['dims'] = dim_dic
            # self.b2fstate = state
        
        elif self.withshift == True and self.withseries == False:
            
            scan = self.data['dircomp']['multi_shift']
            
            state_dic, dim_dic = self.one_dim_scan_b2fstate(iterlist = scan)

            self.data['b2fstate'] = state_dic
            self.data['DefaultSettings']['dims'] = dim_dic
        
        elif self.withshift == False and self.withseries == True:
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            if self.series_flag == 'twin_scan':
                
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
                
                state_dic, dim_dic = self.two_dim_scan_b2fstate(iterlist = scan, 
                                    iterlist_a = ds_key, iterlist_b = ts_key)
            else:
                state_dic, dim_dic = self.one_dim_scan_b2fstate(iterlist = scan)
                

            self.data['b2fstate'] = state_dic
            self.data['DefaultSettings']['dims'] = dim_dic
            
        
        else:
            print('load_b2fstate function is not there yet!')
    
    
    
    
    def load_b2fplasmf(self):
        if self.withshift == False and self.withseries == False:
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'], 'b2fplasmf')
            dim_dic = self.data['DefaultSettings']['dims']
            n_pol = dim_dic['nx']
            n_rad = dim_dic['ny']
            n_sp = dim_dic['ns']
            fplasma = lbdm.read_b2fplasmf(fileName = file_loc, nx = n_pol, 
                                          ny = n_rad, ns = n_sp)
            fplasma_dic = vars(fplasma)
            
            self.data['b2fplasmf'] = fplasma_dic
            # print('the next line is b2fplasmf')
            # print(type(k))
        
        elif self.withshift == True and self.withseries == False:
            fplasma_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], 'b2fplasmf')
                dim_dic = self.data['DefaultSettings']['dims']
                n_pol = dim_dic[aa]['nx']
                n_rad = dim_dic[aa]['ny']
                n_sp = dim_dic[aa]['ns']
                fplasma = lbdm.read_b2fplasmf(fileName = file_loc, nx = n_pol, 
                                              ny = n_rad, ns = n_sp)
                fplasma_dic[aa] = vars(fplasma)
                
                
            self.data['b2fplasmf'] = fplasma_dic
        
        elif self.withshift == False and self.withseries == True:
            fplasma_dic = {}
            
            for aa in list(self.data['dircomp']['Attempt'].keys()):
                file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], 'b2fplasmf')
                dim_dic = self.data['DefaultSettings']['dims']
                n_pol = dim_dic[aa]['nx']
                n_rad = dim_dic[aa]['ny']
                n_sp = dim_dic[aa]['ns']
                fplasma = lbdm.read_b2fplasmf(fileName = file_loc, nx = n_pol, 
                                              ny = n_rad, ns = n_sp)
                fplasma_dic[aa] = vars(fplasma)
                
                
            self.data['b2fplasmf'] = fplasma_dic
            
                
        else:
            print('load_b2fplasmf function is not there yet!')
        
        
    """
    there are two option of fort file to load, one is fort.44.i, the other
    one is fort.46.i
          
    """
    
    def one_dim_scan_ft46(self, iterlist):
        
        ft46_dic = {}
        
        for aa in iterlist:
            
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], 'fort.46.i')
            ft46 = lEdm.read_ft46(fileName = file_loc)
            ft46_dic[aa] = vars(ft46)
            
            
        return ft46_dic
        
    
    
    def two_dim_scan_ft46(self, iterlist, iterlist_a, iterlist_b):
        

        ft46_dic = lmem.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
        
        for tp in iterlist:
            aa = tp[0]
            ab = tp[1]
            
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa][ab], 'fort.46.i')
            ft46 = lEdm.read_ft46(fileName = file_loc)
            ft46_dic[aa][ab] = vars(ft46)
        
        return ft46_dic
    
    
    
    def load_ft46(self):
        
        ftname = 'fort.46.i'
        
        if self.withshift == False and self.withseries == False:
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'], '{}'.format(ftname))
            ft46 = lEdm.read_ft46(fileName = file_loc)
            ft46_dic = vars(ft46)
            
            self.data['ft46'] = ft46_dic
            # print('the next line is b2fplasmf')
            # print(type(k))
        
        elif self.withshift == True and self.withseries == False:
            ft46_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                
                file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], '{}'.format(ftname))
                ft46 = lEdm.read_ft46(fileName = file_loc)
                ft46_dic[aa] = vars(ft46)
                
                
            self.data['ft46'] = ft46_dic
        
        elif self.withshift == False and self.withseries == True:
            # ft46_dic = {}
            
            # for aa in list(self.data['dircomp']['Attempt'].keys()):
                
            #     file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], '{}'.format(ftname))
            #     ft46 = lEdm.read_ft46(fileName = file_loc)
            #     ft46_dic[aa] = vars(ft46)
                
                
            # self.data['ft46'] = ft46_dic
            
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            if self.series_flag == 'twin_scan':
                
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
                
                ft46_dic = self.two_dim_scan_ft46(iterlist = scan, iterlist_a = ds_key, 
                                       iterlist_b = ts_key)
                
                # state_dic, dim_dic = self.two_dim_scan_b2fstate(iterlist = scan, 
                #                     iterlist_a = ds_key, iterlist_b = ts_key)
                
            else:
                
                ft46_dic = self.one_dim_scan_ft46(iterlist = scan)
                # state_dic, dim_dic = self.one_dim_scan_b2fstate(iterlist = scan)
            
            
            self.data['ft46'] = ft46_dic
            
            
        else:
            print('load_b2fplasmf function is not there yet!')
    
    
    
    def one_dim_scan_ft44(self, iterlist):
        
        
        if self.series_compare == True:
            
            ft44_dic = {}
            
            for aa in iterlist:
                
                ft44_cp = {'fixed': {}, 'flux': {}}
                for kk in ['fixed', 'flux']:
                    file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa][kk], 'fort.44.i')
                    ft44 = lEdm.read_ft44(fileName = file_loc)
                    ft44_cp[kk] = vars(ft44)
                
                ft44_dic[aa] = ft44_cp
                
        
        else:
            
            ft44_dic = {}
            
            for aa in iterlist:
                
                file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa], 'fort.44.i')
                ft44 = lEdm.read_ft44(fileName = file_loc)
                ft44_dic[aa] = vars(ft44)
        
        
        
        return ft44_dic
    
    
    def two_dim_scan_ft44(self, iterlist, iterlist_a, iterlist_b):
        


        ft44_dic = lmem.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
        
        for tp in iterlist:
            aa = tp[0]
            ab = tp[1]
            
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'][aa][ab], 'fort.44.i')
            ft44 = lEdm.read_ft44(fileName = file_loc)
            ft44_dic[aa][ab] = vars(ft44)
        
        return ft44_dic
     
    
    
    def load_ft44(self):
        
        ftname = 'fort.44.i'
        
        if self.withshift == False and self.withseries == False:
            file_loc = '{}/{}'.format(self.data['dirdata']['simudir'], '{}'.format(ftname))
            ft44 = lEdm.read_ft44(fileName = file_loc)
            ft44_dic = vars(ft44)
            
            self.data['ft44'] = ft44_dic
            # print('the next line is b2fplasmf')
            # print(type(k))
        
        elif self.withshift == True and self.withseries == False:
            
            
            scan = self.data['dircomp']['multi_shift']
            
            ft44_dic = self.one_dim_scan_ft44(iterlist = scan)
            
            self.data['ft44'] = ft44_dic
        
        elif self.withshift == False and self.withseries == True:
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            if self.series_flag == 'twin_scan':
                
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
                
                ft44_dic = self.two_dim_scan_ft44(iterlist = scan, 
                                    iterlist_a = ds_key, iterlist_b = ts_key)
            
            else:
                ft44_dic = self.one_dim_scan_ft44(iterlist = scan)
                
                            
            self.data['ft44'] = ft44_dic
            
                
        else:
            print('load_b2fplasmf function is not there yet!')
    
    def twinscan_iout(self, iterlist, iterlist_a, iterlist_b, filename):
        

        iout_dic = lmem.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
        # print(iout_dic)
        ny = self.data['b2fgeo']['ny']
        nx = self.data['b2fgeo']['nx']
        
        for tp in iterlist:
            aa = tp[0]
            ab = tp[1]
            
            file_loc = '{}/{}/{}'.format(self.data['dirdata']['simudir'][aa][ab], 'output', filename)
            iout = lbdm.read_iout_method(fdir = file_loc, fname = filename, 
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
        
        
        if self.withshift == False and self.withseries == False:
            file_loc = '{}/{}/{}'.format(self.data['dirdata']['simudir'], 'output', filename)           
            ny = self.data['b2fgeo']['ny']
            nx = self.data['b2fgeo']['nx']
            iout = lbdm.read_iout_method(fdir = file_loc, fname = filename, 
                                  ny = ny, nx = nx)
            
            self.data['iout_data'][quant] = iout

        
        elif self.withshift == True and self.withseries == False:
            iout_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                ny = self.data['b2fgeo'][aa]['ny']
                nx = self.data['b2fgeo'][aa]['nx']
                file_loc = '{}/{}/{}'.format(self.data['dirdata']['simudir'][aa], 'output', filename)           
                iout = lbdm.read_iout_method(fdir = file_loc, fname = filename, 
                                      ny = ny, nx = nx)
                
                iout_dic[aa] = iout
            
            self.data['iout_data'][quant] = iout_dic
                
                
        
        elif self.withshift == False and self.withseries == True:
            
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            
            if self.series_flag == 'twin_scan':
                
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
                iout_dic = self.twinscan_iout(iterlist = scan, iterlist_a = ds_key, 
                              iterlist_b = ts_key, filename = filename)
            else:
                
                iout_dic = {}
                for aa in list(self.data['dircomp']['Attempt'].keys()):
                    ny = self.data['b2fgeo']['ny']
                    nx = self.data['b2fgeo']['nx']
                    
                    file_loc = '{}/{}/{}'.format(self.data['dirdata']['simudir'][aa], 'output', filename)           
                    iout = lbdm.read_iout_method(fdir = file_loc, fname = filename, 
                                          ny = ny, nx = nx)
                    
                    iout_dic[aa] = iout
                
            
            
            self.data['iout_data'][quant] = iout_dic
            
            
        
        else:
            print('load_b2fstate function is not there yet!')
        
        
        return quant
    
    def load_iout_ratio(self, file_tuple, itername):
        
        
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
            
            elif self.series_flag == 'twin_scan':
                
                
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
        
        
        if self.withshift == False and self.withseries == False:
        
            data1 = self.data['iout_data'][name1]
            
            data2 = self.data['iout_data'][name2]
            
            ratio_data = np.divide(data1, data2)
            quantname = '{}_divide_{}'.format(name1, name2)
            
            self.data['iout_data'][quantname] = ratio_data
            # self.data['iout_data']['{}_abs'.format(quantname)] = np.abs(ratio_data)
            
            return quantname
        
        elif self.withshift == True and self.withseries == False:
            
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
        
            
        
            
            
            

            
            
    
    
    
                
                
            

                    
                
                
            
        
    
    
    