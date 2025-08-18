# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 08:47:51 2025

@author: ychuang
"""

from load_directory.load_dirdata_method import load_dir_method
from load_simulation_data.load_B2_data_method import read_b2fstate, read_b2fplasmf, read_b2wdat, read_iout_method
import numpy as np
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt



class load_B2simu_data:
    
    def __init__(self, DF, data, ldm: load_dir_method):
        
        self.DF = DF
        self.data = data
        self.ldm = ldm
    
    
    def main(dir_dat, plotvars):
        file_in = "{}/b2tallies.nc"

        plotvars = plotvars.split()

        try:
            ncIn = Dataset(file_in)
        except:
            print("Error: Could not open "+ file_in)
            exit(0)
            
        times = ncIn.variables['times'][:]

        myvar = {}
        for item in plotvars:
            try:
                myvar[item] = ncIn.variables[item][:]
            except:
                print("Error: Could not load variable "+item)
                exit(0)

    #    print(np.shape(myvar[plotvars[0]]))
        plt.figure()
        for i,item in enumerate(plotvars):
            plt.plot(times,np.squeeze(myvar[plotvars[i]]),label=item,marker='x')
            
        plt.xlabel("time (s)")
        plt.legend(loc="best")
        plt.show()
        

    
    
    







