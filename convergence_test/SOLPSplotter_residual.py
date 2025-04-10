# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 21:59:34 2024

@author: ychuang
"""

import os.path
import numpy as np
import matplotlib.pyplot as plt
import itertools
from SOLPSplotter_ndscan import neuden_scan


class residual_plot(neuden_scan):
    
    def __init__(self, DefaultSettings, loadDS):
        neuden_scan.__init__(self, DefaultSettings, loadDS)



    def load_b2ftrace(self, file_loc):
        
        
        file_in = '{}/{}/{}'.format(file_loc, 'b2mn.exe.dir', 'b2ftrace') 
        if os.path.isfile(file_in):
            print("Reading "+file_in)
        else:
            file_in = '{}/{}'.format(file_loc, 'b2ftrace')
            if os.path.isfile(file_in):
                print("Reading "+file_in)
            else:
                print("Error: Could not find b2ftrace")
                exit(0)
    
        with open(file_in) as f:
            lines = f.readlines()
    
        mydatalist = []
        counter = 0
        read_data = False
        for line in lines:
            if "data" in line:
                counter  += 1
                read_data = True
                continue
            if read_data:
                line = line.split()
                part_list = [float(i) for i in line]
                mydatalist.extend(part_list)
            
        mydata = np.array(mydatalist)
        mydata = mydata.reshape(counter,int(len(mydatalist)/counter))
        
        
        iend = np.size(mydata,0)
        print(iend)
        
        return iend, mydata
    
    
    def plot_b2trace(self, plotvars, scan_style, dat_size):
        
        
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
                    
                    
                    self.plot_b2trace_method(iterlist = iter_key, scandetail = scan_title,
                            cl_dic = color_dic, A_dic = label_dic, plotvars = plotvars,
                            scan_style = scan_style, dat_size = dat_size)
                
           
            
          
    
    
    def plot_b2trace_method(self, plotvars, iterlist, cl_dic, A_dic, scan_style, 
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
                
                
                args_dic = {'conD0': 2, 'conD1': 3, 'momD0': 4, 'momD1': 5, 'totmom': 6,
                            'ee': 7, 'ei': 8, 'phi': 9}
                
                args = args_dic[plotvars]
                
                for aa in iterlist:
                    
                    
                    if scan_style == 'tempscan':
                        
                        ad = aa[1]
                        axs.set_title('{} residual Temperature scan with Ne = {}'.format(plotvars, scandetail))
                    
                    elif scan_style == 'denscan':
                        
                        ad = aa[0]
                        axs.set_title('{} residual Density scan with Te = {} eV'.format(plotvars, scandetail))
                    
                    else:
                        print('neteTSplot_method, please check scan_style')
                    
                    
                    sim_dir = self.data['dirdata']['simudir']
                    f_loc = sim_dir[aa[0]][aa[1]]

                    
                    iend, mydata = self.load_b2ftrace(file_loc = f_loc)
                    
                    # print('check 182:')
                    # print(type(ad))
                    
                    # axs.legend(loc= 'lower left', fontsize=10)
                    axs.plot(mydata[0:iend, args], ls = '-',
                label= '{}'.format(A_dic[ad]), color = cl_dic[ad])
                    
                axs.set_ylabel("norm of residuals")
                axs.set_xlabel("Iteration")
                # axs.set_yscale('log')
                axs.legend(loc="best")