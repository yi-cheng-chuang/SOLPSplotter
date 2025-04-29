# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 17:26:21 2025

@author: ychuang
"""

from midplane_data.SOLPSplotter_PRmap import RP_mapping
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np


class minor_radius_measurement:
    
    
    def __init__(self, DF, data, rp: RP_mapping):
        
        self.DF = DF
        self.data = data
        self.rp = rp
    
    
    
    def calc_Rsep(self, coord_dic, psiN_function, rad_range, maxis_Z):
        
        
        czup = coord_dic[:, 0]
        czlow = coord_dic[:, 1]
        crup = coord_dic[:, 2]
        crlow = coord_dic[:, 3]
        
         
        
        weight_mid = np.zeros(rad_range)
        for x in range(rad_range):
            weight_mid[x] = (maxis_Z - czlow[x])/ (czup[x] - czlow[x])
         
        
        mid_R = np.zeros(rad_range)
        for xa in range(rad_range):
            mid_R[xa] = weight_mid[xa]*crup[xa] + (1 - weight_mid[xa])*crlow[xa]
        
        mid_Z = np.zeros(rad_range)
        for xb in range(rad_range):
            mid_Z[xb] = weight_mid[xb]*czup[xb] + (1 - weight_mid[xb])*czlow[xb]
            
            
        psi_solps_mid = np.zeros(rad_range)
        
        # psi_solps_cp = psiNinterp_RBS(crloc, czloc)
        for i in range(rad_range):
            psi_mid = psiN_function(mid_R[i], mid_Z[i])
            psi_solps_mid[i] = psi_mid
        
        sep_index_dic = self.rp.calc_sep_index(psi = psi_solps_mid, sep_value = 1)
        sep_index_high = int(sep_index_dic['index'][1])
        
        print('index high number is {}'.format(sep_index_high))
        
        weight_psi = (psi_solps_mid[sep_index_high] -1)/(psi_solps_mid[sep_index_high] - psi_solps_mid[sep_index_high -1])
        
        R_sep = (1 - weight_psi)*mid_R[sep_index_high] + weight_psi*mid_R[sep_index_high -1]
        
        print('rsep is {:.3f}'.format(R_sep))
        
        return mid_R, mid_Z, R_sep
    
    
   
    
    def minor_radius_method(self, geo, radgrid, vertgrid, gfile_data, psiNinterp_function, plotRR):
               
        pol_range = int(geo['nx'] + 2)
        rad_range = int(geo['ny'] + 2)
        
        RadLoc = radgrid
        VertLoc = vertgrid
        psiNinterp_RBS = psiNinterp_function
        
        
        # mag_axis_z = gfile_data['zmaxis']
        # mag_axis_z = 0.15
        # print(mag_axis_z)
        
        
        "Calculate midplane R for Z at the magnetic axis"
        
        "Calculate weight"
        
        DEV = self.DF.DEV
        
        if DEV == 'mast':
            
            
            mag_axis_z = gfile_data['zmaxis']
            
            
            
            crup = RadLoc[:, 58]
            crlow = RadLoc[:, 60]
            czup = VertLoc[:, 58]
            czlow = VertLoc[:, 60]
            
            average_pair = (58, 60)
        
        elif DEV == 'mastu':
            
                
            mag_axis_z = gfile_data['zmaxis']
            
            sides = ['LFS', 'HFS']
            
            calc_dic = {}
            
            for ss in sides:
                
                if ss == 'LFS':
                    
                    crup = RadLoc[:, 72]
                    crlow = RadLoc[:, 74]
                    czup = VertLoc[:, 72]
                    czlow = VertLoc[:, 74]
                    
                    average_pair = (72, 74)
                    
                    coord_datadic = np.zeros((rad_range, 4))
                    coord_datadic[:, 0] = czup
                    coord_datadic[:, 1] = czlow
                    coord_datadic[:, 2] = crup
                    coord_datadic[:, 3] = crlow
                    
                    
                    mid_R, mid_Z, R_sep = self.calc_Rsep(coord_dic = coord_datadic, 
                            psiN_function = psiNinterp_RBS, rad_range = rad_range, 
                            maxis_Z = mag_axis_z)
                    
                    calc_dic[ss] = {'LFS_midR': mid_R, 'LFS_midZ': mid_Z, 'LFS_R_sep': R_sep}
                
                
                elif ss == 'HFS':
                    
                    crup = RadLoc[:, 23]
                    crlow = RadLoc[:, 25]
                    czup = VertLoc[:, 23]
                    czlow = VertLoc[:, 25]
                    
                    average_pair = (23, 25)
                    
                    coord_datadic = np.zeros((rad_range, 4))
                    coord_datadic[:, 0] = czup
                    coord_datadic[:, 1] = czlow
                    coord_datadic[:, 2] = crup
                    coord_datadic[:, 3] = crlow
                    
                    
                    mid_R, mid_Z, R_sep = self.calc_Rsep(coord_dic = coord_datadic, 
                            psiN_function = psiNinterp_RBS, rad_range = rad_range, 
                            maxis_Z = mag_axis_z)
                    
                    calc_dic[ss] = {'HFS_midR': mid_R, 'HFS_midZ': mid_Z, 'HFS_R_sep': R_sep}
                   
 
        else:
            print("please check the device input")
        
        
        
        
        
        tot_a = calc_dic['LFS']['LFS_R_sep'] - calc_dic['HFS']['HFS_R_sep']
        majR_loc = calc_dic['HFS']['HFS_R_sep'] + 0.5*tot_a
        
        
        calc_dic['tot_a'] = tot_a
        calc_dic['majR_loc'] = majR_loc
        

        
        if DEV == 'mast':
            pol_list = [57, 58, 59, 60, 61]
            
            
        
        elif DEV == 'mastu':
            
            
            sides = ['LFS', 'HFS']
            
            for ss in sides:
                
                if ss == 'LFS':
                    
                    pol_list = [70, 71, 72, 73, 74]
                    
                    mid_R = calc_dic[ss]['LFS_midR']
                    mid_Z = calc_dic[ss]['LFS_midZ']
                    
                    
                    if plotRR:
                        plt.figure()
                        for in_pol in pol_list:
                            crloc = RadLoc[:, int(in_pol)]
                            czloc = VertLoc[:, int(in_pol)]
                            plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
                        plt.plot(mid_R, mid_Z, color = 'r', label= 'R_Rsep')
                        plt.xlabel('R: [m]')
                        plt.ylabel('Z: [m]')
                 
                    
                elif ss == 'HFS':
                    
                    pol_list = [21, 22, 23, 24, 25]
                    
                    mid_R = calc_dic[ss]['HFS_midR']
                    mid_Z = calc_dic[ss]['HFS_midZ']
                    
                    
                    if plotRR:
                        plt.figure()
                        for in_pol in pol_list:
                            crloc = RadLoc[:, int(in_pol)]
                            czloc = VertLoc[:, int(in_pol)]
                            plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
                        plt.plot(mid_R, mid_Z, color = 'r', label= 'R_Rsep')
                        plt.xlabel('R: [m]')
                        plt.ylabel('Z: [m]')
            
            
            
            
            return calc_dic
            


    def calc_minor_radius(self, plotRR):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        DEV = self.DF.DEV
        
        
        
        if withshift == False and withseries == False:
            b2fgeo = self.data['b2fgeo']
            radloc = self.data['grid']['RadLoc']
            vertloc = self.data['grid']['VertLoc']
            gfile = self.data['gfile']['g']
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            
            
            calc_dic = self.minor_radius_method(geo = b2fgeo, radgrid = radloc, vertgrid = vertloc, 
                gfile_data = gfile, psiNinterp_function = psiNinterp_RBS, 
                        plotRR = plotRR)
            
            
            self.data['calc_minor_radius'] = calc_dic
            
   
        
        elif withshift == True and withseries == False:        
            
            calc_minor_radius = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                
                b2fgeo = self.data['b2fgeo'][aa]
                radloc = self.data['grid']['RadLoc'][aa]
                vertloc = self.data['grid']['VertLoc'][aa]
                gfile = self.data['gfile']['g']
                psiNinterp_RBS = self.data['gfile']['gcomp'][aa]['interp_dic']['RBS']
                
                               
                calc_dic = self.minor_radius_method(geo = b2fgeo, 
                            radgrid = radloc, vertgrid = vertloc, gfile_data = gfile, psiNinterp_function = psiNinterp_RBS, 
                            plotRR = plotRR)
                
                calc_minor_radius[aa] = calc_dic
                
            
            
            self.data['calc_minor_radius'] = calc_minor_radius
            
            
                
           
        elif withshift == False and withseries == True:
            den_rep = list(self.data['dircomp']['Attempt'].keys())[0]
            
            b2fgeo = self.data['b2fgeo']
            radloc = self.data['grid']['RadLoc']
            vertloc = self.data['grid']['VertLoc']
            gfile = self.data['gfile']['g']
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            
            
            self.minor_radius_method(geo = b2fgeo, radgrid = radloc, vertgrid = vertloc, 
                gfile_data = gfile, psiNinterp_function = psiNinterp_RBS, 
                        plotRR = plotRR)
            
            self.data['calc_minor_radius'] = calc_dic
            
          
        
        elif self.withshift == True and self.withseries == True:
            print('calc_RRsep function is not there yet!')
        
        else:
            print('calc_RRsep function has a bug')
            