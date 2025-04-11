# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 14:16:15 2025

@author: ychuang
"""

import numpy as np
from load_directory.load_dirdata_method import load_dir_method


class load_output_data:
    
    def __init__(self, DF, data, ldm: load_dir_method):
        
        self.DF = DF
        self.data = data
        self.ldm = ldm
        
    
    
    def load_output_data_method(self, param, itername):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == False:
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
        elif withshift == True and withseries == False:
            BASEDRT = self.data['dirdata']['outputdir'][itername]['Output']
            Attempt = self.data['dircomp']['Attempt'][itername]
            XGrid = int(self.data['b2fgeo'][itername]['nx'])
            # print(XGrid)
            XDIM = int(self.data['DefaultSettings']['XDIM'][itername])
            YDIM = int(self.data['DefaultSettings']['YDIM'][itername])
            
        elif withshift == False and withseries == True:
            
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
        


        param_data_dic = self.ldm.two_layer_dic(key_a = iterlist_a, key_b = iterlist_b)
        
        for tp in iterlist:
            aa = tp[0]
            ab = tp[1]
            
            param_data_dic[aa][ab] = self.load_output_data_method(param = param, itername = tp)
        
        return param_data_dic
       
    def load_output_data(self, param):
        
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        series_flag = self.DF.series_flag
        
        
        
        if withshift == False and withseries == False:
            output = self.load_output_data_method(param = param, itername = None)
            self.data['outputdata'][param] = output
            
            
        elif withshift == True and withseries == False:
            
            
            scan = self.data['dircomp']['multi_shift']
            
            param_data_dic = self.one_dim_scan_output(iterlist = scan, param = param)
                      
            self.data['outputdata'][param] = param_data_dic
            
            
        elif withshift == False and withseries == True:
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            if series_flag == 'twin_scan':
                
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
        
        elif withshift == True and withseries == True:
            print('load_output_data is not there yet, to be continue...')
        
        else:
            print('There is a bug')