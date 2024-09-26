# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 00:22:11 2024

@author: ychuang
"""


from SOLPSplotter_ndscan import neuden_scan
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt



class flux_plot(neuden_scan):
    
    def __init__(self, DefaultSettings, loadDS):
        neuden_scan.__init__(self, DefaultSettings, loadDS)




    def load_b2fstateflux(self, itername, data_struc):
        
       
        if self.withshift == False and self.withseries == True:
            
            if self.series_flag == 'twin_scan':
                
                nf = itername[0]
                tf = itername[1]
                
                b2fstate = self.data['b2fstate'][nf][tf]
                
                nx = data_struc['nx']
                ny = data_struc['ny']
                
                if data_struc['size'] == 'full':
                    fna = b2fstate['fna'][:, :, 1, 1]
                    fhex = b2fstate['fhe'][:, :, 0]
                    fhey = b2fstate['fhe'][:, :, 1]
                    fhix = b2fstate['fhi'][:, :, 0]
                    fhiy = b2fstate['fhi'][:, :, 1]
                    
                elif data_struc['size'] == 'small':
                    ne_pro = b2fstate['ne'][1:nx+1, 1:ny+1].transpose()
                    Te_J = b2fstate['te'][1:nx+1, 1:ny+1].transpose()
                    fna = b2fstate['fna'][1:nx+1, 1:ny+1, 1, 1]
                    fhex = b2fstate['fhe'][1:nx+1, 1:ny+1, 0]
                    fhey = b2fstate['fhe'][1:nx+1, 1:ny+1, 1]
                    fhix = b2fstate['fhi'][1:nx+1, 1:ny+1, 0]
                    fhiy = b2fstate['fhi'][1:nx+1, 1:ny+1, 1]
                
                
                flux_dic = {'fna': fna, 'fhex': fhex, 'fhey': fhey, 'fhix': fhix,
                            'fhiy': fhiy}
                
                return flux_dic
        
        
        
        
    def load_b2time(self, plotvars, scan_style, dat_size):
        
        
        if self.withshift == False and self.withseries == True:
            
            if self.series_flag == 'twin_scan':
                
           
                dircomp = self.data['dircomp']
                
                if scan_style == 'tempscan':
                    
                    key_a = 'denscan_list'
                    key_b = 'tempscan_list'
                
                elif scan_style == 'denscan':
                    
                    key_a = 'tempscan_list'
                    key_b = 'denscan_list'
                
                else:
                    print('twinscan_plot_method, please check the scan_style!')
                
                keylist_a = []
                
                
                for x in dircomp[key_a]:
                    keylist_a.append('{:.3f}'.format(x))
                
                for ta in keylist_a:
                    
                    keylist_b = []
                    
                    for x in dircomp[key_b]:
                        keylist_b.append('{:.3f}'.format(x))
                    
                        
                    iter_key, color_dic, scan_title, label_dic = self.twinscan_prep(ta = ta, 
                    keylist_b = keylist_b, scan_style = scan_style, dat_size = dat_size)
                    
                    
                    # print('check:')
                    # print(iter_key)
                    # print(color_dic)
                    # print(label_dic)
                    
                    
                    self.plot_b2time_method(iterlist = iter_key, scandetail = scan_title,
                            cl_dic = color_dic, A_dic = label_dic, plotvars = plotvars,
                            scan_style = scan_style, dat_size = dat_size)
                
           
            
          
    
    
    def plot_b2fstateflux_method(self, plotvars, iterlist, cl_dic, A_dic, scan_style, 
                          scandetail, dat_size):
        
        
        
        if self.withshift == False and self.withseries == True:
            
            
            if self.series_flag == 'twin_scan':
            
                fig, axs = plt.subplots()
                
                
                nx = self.data['b2fgeo']['nx']
                ny = self.data['b2fgeo']['ny']
                
                
                if dat_size == 'full':
        
                    dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
                
                elif dat_size == 'small':
                    dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
                
                
                
                for aa in iterlist:
                    
                    
                    if scan_style == 'tempscan':
                        
                        ad = aa[1]
                        axs.set_title('{} Temperature scan with Ne = {}'.format(plotvars, scandetail))
                    
                    elif scan_style == 'denscan':
                        
                        ad = aa[0]
                        axs.set_title('{} Density scan with Te = {} eV'.format(plotvars, scandetail))
                    
                    else:
                        print('neteTSplot_method, please check scan_style')
                    

                    
                    flux_dic = self.load_b2fstateflux(itername = aa, data_struc = dat_struc)
                    

                    


                    axs.plot(timesa, np.squeeze(myvar[plotvars]), 
                label= '{}'.format(A_dic[ad]), ls = '-', color = cl_dic[ad])
                    
            
                axs.set_xlabel("time (s)")
                axs.legend(loc="best")