# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:00:38 2024

@author: user
"""
from B2plotter_class import B2plotter
import re
import transport_coefficient_adjust_method as tcam
import matplotlib.pyplot as plt
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import numpy as np
import xarray as xr
import math




class transport_coefficient_adjustment(B2plotter):
    
    def __init__(self, DEV, withshift, withseries, DefaultSettings, loadDS):
        B2plotter.__init__(self, DEV, withshift, withseries, DefaultSettings)
        # Employee.__init__(self, first, last, pay)
        self.loadDS = loadDS
    
       
    def mod_transco(self, withmod, de_SOL = 26, ki_SOL = 31, ke_SOL = 25, ):
        
        simu_dir = self.data['dirdata']['infolderdir']['org']['simudir']
        file_loc = '{}/b2.transport.inputfile_new'.format(simu_dir)

        trans_list = tcam.load_transcoefile_method(file_loc, plot= False)
        cod = trans_list['1'].T
        coki = trans_list['3'].T
        coke = trans_list['4'].T
        x= cod[:,0]  #the coordinate here is R-R_sep
        yd= cod[:,1]
        yki = coki[:,1]
        yke = coke[:,1]

        m = len(yd)
        if withmod:
            mod_y = np.zeros(m)
            for j in range(m):
                if j<= 26:
                    mod_y[j] = cod[j,1]
                else:
                    mod_y[j] = 15.0
            cod[:,1] = mod_y

            mod_yki = np.zeros(m)
            for j in range(m):
                if j<= 31:
                    mod_yki[j] = coki[j,1]  
                else:
                    mod_yki[j] = 10.0
            coki[:,1] = mod_yki

            mod_yke = np.zeros(m)
            for j in range(m):
                if j<= 25:
                    mod_yke[j] = coke[j,1]  
                else:
                    mod_yke[j] = 18.0
            coke[:,1] = mod_yke
        else:
            pass


        b = tcam.Generate_transcoefile_method(cod, CoeffID=1, SpeciesID=2, M=[1])
        
        c = b2tp.WriteInputfile(file = 'b2.transport.inputfile_mod_{}{}'.format(shift, n), points= trans_list ,M_1 = True, M=[1])

        log_flag = False
        specieslist = ['1','3','4']
        d = tl.unit_dic()
        i = 0

        for k in specieslist:
            if log_flag:
                plt.yscale('log')
                plt.figure(i + 1)
                plt.plot(trans_list[k][0,:], trans_list[k][1,:], 'o-', color = 'orange', label ='{} transport coefficient'.format(filename))
                plt.xlabel('Radial coordinate: $R- R_{sep}$')
                plt.ylabel(d[k][1])
                plt.title(d[k][0])
                plt.legend()
            else:
                plt.figure(i + 1)
                plt.plot(trans_list[k][0,:], trans_list[k][1,:], 'o-', color = 'orange', label ='{} transport coefficient'.format(filename))
                plt.xlabel('Radial coordinate: $R- R_{sep}$')
                plt.ylabel(d[k][1])
                plt.title(d[k][0])
                plt.legend()
            i = i + 1

        plt.show()
    
    
    