# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 22:41:18 2024

@author: user
"""

from SOLPSplotter_PRmap import RP_mapping
import opacity_plot_method as opm
import matplotlib.pyplot as plt
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate
from scipy.optimize import curve_fit
import numpy as np


class sep_data_process(RP_mapping):
    
    def __init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters):
        RP_mapping.__init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters)
    
    
    
    """
    
    This is the first version of filter, we can later isolate a 
    general function from it so that can be used for b2fstate
    
    """
        
    
    def b2f_filter_method(self, datafile, dim_setting):
        nxny_list = []
        zero_nxny_list = []
               
        nxnyns_list = []
        zero_nxnyns_list = []
        
        nxny_corner_list = []
        zero_nxnycorner_list = []
        
        fluxdim_ns_list = []
        zero_fluxdimns_list = []
        
        nxny_corner_ns_list = []
        zero_nxnycornerns_list = []
        
        rest_list = []
        
        nx = dim_setting['nx']
        ny = dim_setting['ny']
        ns = dim_setting['ns']
        nc = 4
               
        # b2plasmf = datafile
        for plasmf_k in list(datafile.keys()):
            if np.shape(datafile[plasmf_k]) == (nx+2, ny+2):
                if np.all(datafile[plasmf_k] == 0):
                    zero_nxny_list.append(plasmf_k)
                
                else:
                    nxny_list.append(plasmf_k)
            
            elif np.shape(datafile[plasmf_k]) == (nx+2, ny+2, ns):
                for ns_a in range(ns):
                    if np.all(datafile[plasmf_k][:, :, ns_a] == 0):
                        item_dic = {'item': plasmf_k, 'ns': ns_a}
                        zero_nxnyns_list.append(item_dic)
                    
                    else:
                        item_dic = {'item': plasmf_k, 'ns': ns_a}
                        nxnyns_list.append(item_dic)
                    
                
            
            elif np.shape(datafile[plasmf_k]) == (nx+2, ny+2, nc):
                
                for nc_a in range(4):
                    if np.all(datafile[plasmf_k][:, :, nc_a] == 0):
                        
                        item_dic = {'item': plasmf_k, 'nc': nc_a}
                        zero_nxnycorner_list.append(item_dic)
                    
                    else:
                        
                        item_dic = {'item': plasmf_k, 'nc': nc_a}
                        nxny_corner_list.append('{}%{}'.format(plasmf_k, nc_a))
            
            
            
            elif np.shape(datafile[plasmf_k]) == (nx+2, ny+2, 2, 2):
                
                for ns_b in range(ns):
                    for nf_a in range(2):
                        if np.all(datafile[plasmf_k][:, :, nf_a, ns_b] == 0):
                            
                            item_dic = {'item': plasmf_k, 'nf': nf_a, 'ns': ns_b}
                            zero_fluxdimns_list.append(item_dic)
                        
                        else:
                            
                            item_dic = {'item': plasmf_k, 'nf': nf_a, 'ns': ns_b}
                            fluxdim_ns_list.append(item_dic)
                          
            elif np.shape(datafile[plasmf_k]) == (nx+2, ny+2, nc, ns):
                for ns_c in range(ns):
                    for nc_b in range(nc):
                        if np.all(datafile[plasmf_k][:, :, nc_b, ns_c] == 0):
                            
                            item_dic = {'item': plasmf_k, 'nc': nc_b, 'ns': ns_c}
                            zero_nxnycornerns_list.append(item_dic)
                        
                        else:
                            
                            item_dic = {'item': plasmf_k, 'nc': nc_b, 'ns': ns_c}
                            nxny_corner_ns_list.append(item_dic)
                           
            
            else:
                rest_list.append(plasmf_k)
                
                    
        
        key_order_dic = {'nxny': nxny_list, 'zero_nxny': zero_nxny_list, 
                         'nxnyns': nxnyns_list, 'zero_nxnyns': zero_nxnyns_list,
                         'nxnycorner': nxny_corner_list, 'zero_nxnycorner': zero_nxnycorner_list,
                         'fluxdim_ns': fluxdim_ns_list, 'zero_fluxdimns': zero_fluxdimns_list, 
                         'nxny_corner_ns': nxny_corner_ns_list, 
                         'zero_nxnycornerns': zero_nxnycornerns_list, 'the_rest': rest_list}
        
        
        return key_order_dic
    
                    
    
    
                    
    def b2f_file_filter(self, b2f_file, b2f_name):
        if self.withshift == False and self.withseries == False:
                       
            # b2plasmf = self.data['b2fplasmf']
            dim = self.data['DefaultSettings']['dims']
            key_order = self.b2f_filter_method(datafile = b2f_file, 
                                                             dim_setting = dim)
            
            self.data['b2fplasmf_key'] = key_order
            
        elif self.withshift == True and self.withseries == False:
            
            key_order_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                                
                b2plasmf = b2f_file[aa]
                dim = self.data['DefaultSettings']['dims'][aa]
                key_order_dic[aa] = self.b2f_filter_method(datafile = b2plasmf, 
                                                                 dim_setting = dim)
                
                
            self.data['b2fplasmf_key'] = key_order_dic
        
        elif self.withshift == False and self.withseries == True:
            
            key_order_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                                
                b2plasmf = b2f_file[aa]
                dim = self.data['DefaultSettings']['dims'][aa]
                key_order_dic[aa] = self.b2f_filter_method(datafile = b2plasmf, 
                                                                 dim_setting = dim)
                
                
            self.data['{}_key'.format(b2f_name)] = key_order_dic
    
    
    def nxny_sep_data_process(self, specshape_list, shape_spec):
        
        if self.withshift == False and self.withseries == False:
            
            
            # nxny_list = self.data['b2fplasmf_key']['nxny']
            nxny_list = specshape_list[shape_spec]
            
            sep_index = self.data['DefaultSettings']['sep_index_dsa']
            
            sep_data_dic = {}
            
            for item in nxny_list:
                
                if shape_spec == 'nxny':
                    
                    data_U = self.data['b2fplasmf'][item][:, sep_index]
                    data_L = self.data['b2fplasmf'][item][:, sep_index - 1]
                    
                    sep_data = np.mean([data_U, data_L], axis= 0)
                    
                    sep_data_dic[item] = sep_data
                
                elif shape_spec == 'nxnyns':
                    
                    ns = item['ns']
                    plasma_k = item['item']
                                       
                    data_U = self.data['b2fplasmf'][plasma_k][:, sep_index, ns]
                    data_L = self.data['b2fplasmf'][plasma_k][:, sep_index - 1, ns]
                    
                    sep_data = np.mean([data_U, data_L], axis= 0)
                    
                    sep_data_dic['{}%ns{}'.format(plasma_k, ns)] = sep_data
                
                elif shape_spec == 'fluxdim_ns':
                    
                    nf = item['nf']
                    ns = item['ns']
                    plasma_k = item['item']
                                       
                    data_U = self.data['b2fplasmf'][plasma_k][:, sep_index, nf, ns]
                    data_L = self.data['b2fplasmf'][plasma_k][:, sep_index - 1, nf, ns]
                    
                    sep_data = np.mean([data_U, data_L], axis= 0)
                    
                    sep_data_dic['{}%nf{}%ns{}'.format(plasma_k, nf, ns)] = sep_data
                
                
                
            
            self.data['{}_sep_data'.format(shape_spec)] = sep_data_dic
        
        elif self.withshift == True and self.withseries == False:
            
            sep_data_all = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                
            
                # nxny_list = self.data['b2fplasmf_key'][aa]['nxny']
                
                nxny_list = specshape_list[aa][shape_spec]
                               
                sep_index = self.data['DefaultSettings']['sep_index_dsa'][aa]
                
                sep_data_dic = {}
                
                for item in nxny_list:
                    
                    if shape_spec == 'nxny':
                        
                        data_U = self.data['b2fplasmf'][aa][item][:, sep_index]
                        data_L = self.data['b2fplasmf'][aa][item][:, sep_index - 1]
                        
                        sep_data = np.mean([data_U, data_L], axis= 0)
                        
                        sep_data_dic[item] = sep_data
                    
                
                    elif shape_spec == 'nxnyns':
                        
                        ns = item['ns']
                        plasma_k = item['item']
                                           
                        data_U = self.data['b2fplasmf'][aa][plasma_k][:, sep_index, ns]
                        data_L = self.data['b2fplasmf'][aa][plasma_k][:, sep_index - 1, ns]
                        
                        sep_data = np.mean([data_U, data_L], axis= 0)
                        
                        sep_data_dic['{}%ns{}'.format(plasma_k, ns)] = sep_data
                    
                    elif shape_spec == 'fluxdim_ns':
                        
                        nf = item['nf']
                        ns = item['ns']
                        plasma_k = item['item']
                                           
                        data_U = self.data['b2fplasmf'][aa][plasma_k][:, sep_index, nf, ns]
                        data_L = self.data['b2fplasmf'][aa][plasma_k][:, sep_index - 1, nf, ns]
                        
                        sep_data = np.mean([data_U, data_L], axis= 0)
                        
                        sep_data_dic['{}%nf{}%ns{}'.format(plasma_k, nf, ns)] = sep_data
                    
                    
                
                sep_data_all[aa] = sep_data_dic
                
            
            self.data['{}_sep_data'.format(shape_spec)] = sep_data_all
            
            
                   
        else:
            
            print('nxny_sep_process function is not there yet!')
            
            
                
            
            
            
        