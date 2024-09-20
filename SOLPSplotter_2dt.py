# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 20:02:33 2024

@author: ychuang
"""




#!/usr/bin/env python
#
# run as, for example: 2dt.py "nesepa nesepm nesepi"
#
# Have not checked yet against >2d arrays
#
# JDL


from SOLPSplotter_ndscan import neuden_scan
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt



class twodt_plot(neuden_scan):
    
    def __init__(self, DefaultSettings, loadDS):
        neuden_scan.__init__(self, DefaultSettings, loadDS)



    def load_b2time_method(self, file_loc, plotvars):
        
        file_in = '{}/{}'.format(file_loc, 'b2time.nc')
    
        plotvars = plotvars.split()
    
        try:
            ncIn = Dataset(file_in)
        except:
            print("Error: Could not open "+file_in)
            exit(0)
            
        timesa = ncIn.variables['timesa'][:]
    
        myvar = {}
        for item in plotvars:
            try:
                myvar[item] = ncIn.variables[item][:]
            except:
                print("Error: Could not load variable "+item)
                exit(0)
        
        
        print(np.shape(myvar[plotvars[0]]))
    
        
        return timesa, myvar
    
    
    
    
    
    
    
    
    
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
                
           
            
          
    
    
    def plot_b2time_method(self, plotvars, iterlist, cl_dic, A_dic, scan_style, 
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
                    
                    
                    sim_dir = self.data['dirdata']['simudir']
                    f_loc = sim_dir[aa[0]][aa[1]]

                    
                    timesa, myvar = self.load_b2time_method(file_loc = f_loc, plotvars = plotvars)
                    
                    # print('check 182:')
                    # print(type(ad))
                    
                    # axs.legend(loc= 'lower left', fontsize=10)
                    axs.plot(timesa, np.squeeze(myvar[plotvars]), 
                label= '{}'.format(A_dic[ad]), ls = '-', color = cl_dic[ad])
                    
            
                axs.set_xlabel("time (s)")
                axs.legend(loc="best")
    
    
    
    





    
    
    

