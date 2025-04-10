# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 22:35:37 2024

@author: ychuang
"""



#! /usr/bin/env python
import netCDF4
import os
import matplotlib.pyplot as plt
import numpy
import subprocess
import re
from SOLPSplotter_ndscan import neuden_scan



class PB_plot(neuden_scan):
    
    def __init__(self, DefaultSettings, loadDS):
        neuden_scan.__init__(self, DefaultSettings, loadDS)


    def load_b2tallies(self, file_loc):

        
        file_in = '{}/{}/{}'.format(file_loc, 'b2mn.exe.dir', 'b2tallies.nc') 
        if os.path.isfile(file_in):
            f = netCDF4.Dataset(file_in,'r')
            print("Reading "+file_in)
        else:
            file_in = '{}/{}'.format(file_loc, 'b2tallies.nc')
            if os.path.isfile(file_in):
                f = netCDF4.Dataset(file_in,'r')
                print("Reading "+file_in)
            else:
                
                print("Error: Could not find 'b2tallies.nc'")
                exit(0)
        
        return f
        
        
        
        # if os.access('b2mn.exe.dir/b2tallies.nc', os.R_OK):
        #     f = netCDF4.Dataset('b2mn.exe.dir/b2tallies.nc','r')
        # else:
        #     f = netCDF4.Dataset('b2tallies.nc','r')
        
        
              

    def plot_PB(self, scan_style, dat_size):
        
        
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
                                      
                    self.plot_PB_method(iterlist = iter_key, scandetail = scan_title,
                            cl_dic = color_dic, A_dic = label_dic,
                            scan_style = scan_style, dat_size = dat_size)
                
           
            
          
    
    
    def plot_PB_method(self, iterlist, cl_dic, A_dic, scan_style, 
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
                        axs.set_title('{} Temperature scan with Ne = {}'.format('Normalised particle error', scandetail))
                    
                    elif scan_style == 'denscan':
                        
                        ad = aa[0]
                        axs.set_title('{} Density scan with Te = {} eV'.format('Normalised particle error', scandetail))
                    
                    else:
                        print('neteTSplot_method, please check scan_style')
                    
                    
                    sim_dir = self.data['dirdata']['simudir']
                    f_loc = sim_dir[aa[0]][aa[1]]

                    
                    f = self.load_b2tallies(file_loc = f_loc)
                    
                    vreg = f.dimensions['vreg'].size
                    xreg = f.dimensions['xreg'].size
                    yreg = f.dimensions['yreg'].size
                    ns = f.dimensions['ns'].size
                    time = f.dimensions['time'].size
                    times = f.variables['times']
                    species_names=f.variables['species']
                    species=[b''.join(species_names[i,:]).strip().decode('utf-8') for i in range(species_names.shape[0])]
                    elements=[re.sub('[-+0-9]','',species[i]) for i in range(len(species))]
                    mask=[re.match('[a-zA-Z]+0$',species[i])!=None for i in range(len(species))]
                    s=0
                    bounds=[]
                    for i in range(mask.count(True)-1):
                        e=mask.index(True,s+1)
                        bounds+=[[s,e]]
                        s=e
                    bounds+=[[s,len(mask)]]
                    
                    fnaxreg = f.variables['fnaxreg']
                    fnayreg = f.variables['fnayreg']
                    b2stbr_sna_reg = f.variables['b2stbr_sna_reg']
                    b2sext_sna_reg = f.variables['b2sext_sna_reg']
                    
                    if vreg == 5:
                        FULL_X = numpy.array([0,1,0,0,-1,0,0])
                        FULL_Y = numpy.array([0,1,1,1,0,-1,-1,-1])
                    elif vreg ==2:
                        FULL_X = numpy.array([0,1,-1])
                        FULL_Y = numpy.array([0,1,-1])
                    else:
                        raise ValueError('Value of vreg=%s not currently coded' % vreg)
                    
                    for i in range(len(bounds)):
                        fnax = (fnaxreg[:,bounds[i][0]:bounds[i][1],:].sum(axis=1)*FULL_X).sum(axis=1)
                        fnay = (fnayreg[:,bounds[i][0]:bounds[i][1],:].sum(axis=1)*FULL_Y).sum(axis=1)
                        b2stbr = (b2stbr_sna_reg[:,bounds[i][0]:bounds[i][1],1:].copy())
                        b2stbr[:,0,:] = 0 # remove the neutral
                        if b2stbr.shape[1] > 2: b2stbr[:,-1, :] = 0 # remove the fully stripped species for He and up
                        b2stbr = (b2stbr.sum(axis=1)).sum(axis=1)
                    
                        b2sext = b2sext_sna_reg[:,bounds[i][0]:bounds[i][1],1:].copy()
                        b2sext = (b2sext.sum(axis=1)).sum(axis=1)
                    
                        fna_norm = numpy.max([numpy.max(numpy.abs(fnaxreg[:,bounds[i][0]:bounds[i][1],:])),numpy.max(numpy.abs(fnayreg[:,bounds[i][0]:bounds[i][1],:]))])
                        print(elements[bounds[i][0]], fna_norm, (fnax+fnay+b2stbr+b2sext).mean(), ((fnax+fnay+b2stbr+b2sext)/fna_norm).mean())
                        axs.plot(times[:],((fnax+fnay+b2stbr+b2sext)/fna_norm),  
                        label= '{}'.format(A_dic[ad]), ls = '-', color = cl_dic[ad])
                    
            
                axs.set_xlabel("time (s)")
                axs.set_ylabel('normalised particle error')
                axs.legend(loc="best")