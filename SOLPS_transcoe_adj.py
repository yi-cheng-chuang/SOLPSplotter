# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:00:38 2024

@author: user
"""
from SOLPSplotter_geo import load_geometry
import SOLPS_set as sps
import transport_coefficient_adjust_method as tcam
import matplotlib.pyplot as plt
import numpy as np




class transport_coefficient_adjustment(load_geometry):
    
    def __init__(self, DEV, withshift, withseries, DefaultSettings):
        load_geometry.__init__(self, DEV, withshift, 
                                        withseries, DefaultSettings)
    
       
    def mod_transco(self, withmod, de_SOL, ki_SOL, ke_SOL, log_flag):
        
        
        
        folder = '71_n100000_m12n8e3_nts5_a'
        simu_dir = '{}/{}'.format(self.data['dirdata']['simutop'], folder)
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
                if j<= de_SOL:
                    mod_y[j] = cod[j,1]
                else:
                    mod_y[j] = 12.0
            cod[:,1] = mod_y

            mod_yki = np.zeros(m)
            for j in range(m):
                if j<= ki_SOL:
                    mod_yki[j] = coki[j,1]  
                else:
                    mod_yki[j] = 10.0
            coki[:,1] = mod_yki

            mod_yke = np.zeros(m)
            for j in range(m):
                if j<= ke_SOL:
                    mod_yke[j] = coke[j,1]  
                else:
                    mod_yke[j] = 18.0
            coke[:,1] = mod_yke
        else:
            pass


        tcam.Generate_transcoefile_method(cod, CoeffID=1, SpeciesID=2, M=[1])
        
        shift = 'org'
        n = str(sps.s_number(file_loc, series_flag= None)[0])
        
        tcam.Write_transcoefile_method(file = '{}/b2.transport.inputfile_mod_{}{}'.format(simu_dir, shift, n), points= trans_list ,M_1 = True, M=[1])

        log_flag = False
        specieslist = ['1','3','4']
        transcoe_unit = sps.transport_coe_unit()

        for k in specieslist:
            if log_flag:
                plt.yscale('log')
            else:
                pass
            plt.figure(figsize=(7,7))
            plt.plot(trans_list[k][0,:], trans_list[k][1,:], 'o-', color = 'orange')
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
            plt.ylabel(transcoe_unit[k][1])
            plt.title(transcoe_unit[k][0])
            plt.legend()

        plt.show()
    
    
    