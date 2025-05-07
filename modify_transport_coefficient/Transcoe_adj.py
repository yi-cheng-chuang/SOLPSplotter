# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:00:38 2024

@author: user
"""


from matplotlib.offsetbox import AnchoredText
from modify_transport_coefficient.transport_coefficient_adjust_method import transcoe_method
from load_directory.grab_attempt_number import grab_aptn_method
from SOLPS_input.input_setting import set_figdir
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np



class transport_coefficient_adjustment:
    
    def __init__(self, DF, data, gam: grab_aptn_method, tm: transcoe_method):
        
        self.DF = DF
        self.data = data
        self.gam = gam
        self.tm = tm
    
      
    
    def mod_transco_method(self,file_loc, withmod, de_SOL, ki_SOL, ke_SOL, log_flag):
        
        
        DEV = self.DF.DEV
        
        trans_list = self.tm.load_transcoefile_method(file_loc, plot= False)
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
                    mod_y[j] = 5.0
            cod[:,1] = mod_y

            mod_yki = np.zeros(m)
            for j in range(m):
                if j<= ki_SOL:
                    mod_yki[j] = coki[j,1]  
                else:
                    mod_yki[j] = 20.0
            coki[:,1] = mod_yki

            mod_yke = np.zeros(m)
            for j in range(m):
                if j<= ke_SOL:
                    mod_yke[j] = coke[j,1]  
                else:
                    mod_yke[j] = 20.0
            coke[:,1] = mod_yke
        else:
            pass


        self.tm.Generate_transcoefile_method(cod, CoeffID=1, SpeciesID=2, M=[1])
        
        shift = 'org'
        
        if DEV == 'mast':
            n = str(self.gam.s_number(file_loc)[0])
        
        elif DEV == 'mastu':
            
            n = str(self.gam.mastu_atp_number(file_loc, usage = 'transcoe')[0])
            

        simu_dir = file_loc.rsplit("/",2)[0]
        # k = str(ss.s_number(file_loc, series_flag= None))
        print(simu_dir)
        
        self.tm.Write_transcoefile_method(file = '{}/b2.transport.inputfile_mod_{}{}'.format(simu_dir, shift, n), points= trans_list ,M_1 = True, M=[1])


        specieslist = ['1','3','4']
        transcoe_unit = self.tm.transport_coe_unit()

        for k in specieslist:
            
            plt.figure()
            
            if log_flag:
                plt.yscale('log')
            else:
                pass
            plt.plot(trans_list[k][0,:], trans_list[k][1,:], 'o-', color = 'orange')
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
            # plt.ylabel(transcoe_unit[k][1])
            plt.title(transcoe_unit[k][0])

        plt.show()
        

        
    def transcoe_detailmod_method(self,file_loc, withmod, de_cut, ki_cut, ke_cut, log_flag, ped_adj):
        
        
        DEV = self.DF.DEV
        
        trans_list = self.tm.load_transcoefile_method(file_loc, plot= False)
        cod = trans_list['1'].T
        coki = trans_list['3'].T
        coke = trans_list['4'].T
        x= cod[:,0]  #the coordinate here is R-R_sep
        yd= cod[:,1]
        yki = coki[:,1]
        yke = coke[:,1]
        
        
        D_ped_adj = ped_adj[0]
        ki_ped_adj = ped_adj[1]
        ke_ped_adj = ped_adj[2]
        

        m = len(yd)
        if withmod:
            
            de_ped = de_cut[0]
            de_SOL = de_cut[1]
            
            mod_y = np.zeros(m)
            for j in range(m):
                
                
                if j>= de_ped and j<= de_SOL:
                    mod_y[j] = cod[j,1]
                    
                elif j< de_ped:
                    
                    if D_ped_adj:
                        mod_y[j] = 5.0
                    
                    else:
                        mod_y[j] = cod[j,1]
                        

                    
                elif j> de_SOL:
                    mod_y[j] = 20.0
                
                    
            cod[:,1] = mod_y
            
            
            ki_ped = ki_cut[0]
            ki_SOL = ki_cut[1]
            
            
            mod_yki = np.zeros(m)
            for j in range(m):
                
                
                if j>= ki_ped and j<= ki_SOL:
                    mod_yki[j] = coki[j,1]
                
                elif j< ki_ped:
                    
                    if ki_ped_adj:
                        mod_yki[j] = 14.0
                    
                    else:
                        mod_yki[j] = coki[j,1]
                
                
                elif j> ki_SOL:
                    
                    mod_yki[j] = 20.0
            
            
            coki[:,1] = mod_yki
            
            
            ke_ped = ke_cut[0]
            ke_SOL = ke_cut[1]
            
            
            
            mod_yke = np.zeros(m)
            for j in range(m):
                
                
                if j>= ke_ped and j<= ke_SOL:
                    
                    mod_yke[j] = coke[j,1]
                
                elif j< ke_ped:
                    
                    if ke_ped_adj:
                        mod_yke[j] = 20.0
                    
                    else:
                        mod_yke[j] = coke[j,1]
                
                
                elif j> ke_SOL:
                    
                    mod_yke[j] = 20.0
            
            
            
            
            
            coke[:,1] = mod_yke
        
        
        
        
        else:
            pass


        self.tm.Generate_transcoefile_method(cod, CoeffID=1, SpeciesID=2, M=[1])
        
        shift = 'org'
        
        if DEV == 'mast':
            n = str(self.gam.s_number(file_loc)[0])
        
        elif DEV == 'mastu':
            
            n = str(self.gam.mastu_atp_number(file_loc, usage = 'transcoe')[0])
            

        simu_dir = file_loc.rsplit("/",2)[0]
        # k = str(ss.s_number(file_loc, series_flag= None))
        print(simu_dir)
        
        self.tm.Write_transcoefile_method(file = '{}/b2.transport.inputfile_mod_{}{}'.format(simu_dir, shift, n), points= trans_list ,M_1 = True, M=[1])

        specieslist = ['1','3','4']
        transcoe_unit = self.tm.transport_coe_unit()

        for k in specieslist:
            
            plt.figure()
            
            if log_flag:
                plt.yscale('log')
            else:
                pass
        
            # plt.yscale('log')
            plt.plot(trans_list[k][0,:], trans_list[k][1,:], 'o-', color = 'orange')
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
            # plt.ylabel(transcoe_unit[k][1])
            plt.title(transcoe_unit[k][0])

        plt.show()




        
        
    def mod_transco(self, withmod, de_SOL, ki_SOL, ke_SOL, log_flag, modnew):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        
        if withshift == False and withseries == False:
            simudir = self.data['dirdata']['simudir']
            
            if modnew:
                fileloc = '{}/b2.transport.inputfile_new'.format(simudir)
            else:
                fileloc = '{}/b2.transport.inputfile'.format(simudir)
                
            self.mod_transco_method(file_loc = fileloc, withmod = withmod, de_SOL = de_SOL, 
                                    ki_SOL = ki_SOL, ke_SOL = ke_SOL, log_flag = log_flag)
        
        elif withshift == True and withseries == False:
            simudir = self.data['dirdata']['simudir']['org']
            
            if modnew:
                fileloc = '{}/b2.transport.inputfile_new'.format(simudir)
            else:
                fileloc = '{}/b2.transport.inputfile'.format(simudir)
                
            self.mod_transco_method(file_loc = fileloc, withmod = withmod, de_SOL = de_SOL, 
                                    ki_SOL = ki_SOL, ke_SOL = ke_SOL, log_flag = log_flag)
    
    
    def transco_mod_detail(self, withmod, de_cut, ki_cut, ke_cut, log_flag, modnew, ped_adj):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        
        if withshift == False and withseries == False:
            simudir = self.data['dirdata']['simudir']
            
            if modnew:
                fileloc = '{}/b2.transport.inputfile_new'.format(simudir)
            else:
                fileloc = '{}/b2.transport.inputfile'.format(simudir)
            
            self.transcoe_detailmod_method(file_loc = fileloc, withmod = withmod, de_cut = de_cut,  
                                    ki_cut = ki_cut, ke_cut = ke_cut, log_flag = log_flag, ped_adj = ped_adj)
        
        elif withshift == True and withseries == False:
            simudir = self.data['dirdata']['simudir']['org']
            
            
            if modnew:
                fileloc = '{}/b2.transport.inputfile_new'.format(simudir)
            else:
                fileloc = '{}/b2.transport.inputfile'.format(simudir)
            
            
            self.transcoe_detailmod_method(file_loc = fileloc, withmod = withmod, de_cut = de_cut,  
                                    ki_cut = ki_cut, ke_cut = ke_cut, log_flag = log_flag, ped_adj = ped_adj)
        
        
        

    








    
    
    