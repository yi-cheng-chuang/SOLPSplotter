# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 01:55:05 2025

@author: ychuang
"""


import matplotlib.pyplot as plt 
import numpy as np


class bin_plot:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data

    
    def bin_with_angle(self):
        
        test_list = np.linspace(250, -100, 8)
        # print(test_list)
        
        angle_list = self.data['angle']['angle_list']
        
        angle_flist = [round(an, 2) for an in angle_list]
        # print(angle_flist)
        
        
        twlayer_list = []
        twind_list = []
        for em in range(len(test_list)-1):
            twlayer_list.append([])
            twind_list.append([])
        
        # print(twlayer_list)

        
        for in_an, an in enumerate(angle_flist):
            
            
            for i in range(len(test_list)-1):
                
                if an >= test_list[i + 1] and an<= test_list[i]:
                    
                    twlayer_list[i].append(an)
                    twind_list[i].append(in_an)
        
        # print(twind_list)
        
        return twind_list
    
    
    
    def bin_check(self, pol_list):
        
        twindex_list = self.bin_with_angle()
        
        # print(twindex_list)
        
        
        kp = len(twindex_list)
        
        # print(f'the length of k_list is: {kp}')
        
        
        twpol_list = []
        
        for em in range(7):
            twpol_list.append([])

            
        for i, kp in enumerate(twindex_list):
            
            for kc in kp:
                
                twpol_list[i].append(int(pol_list[kc]))
        
        
        print(twpol_list)
        
        
        angsep_list = np.linspace(250, -100, 8)
        marker_label = []
        for ai in range(len(angsep_list) -1):
            
            aup = angsep_list[ai]
            adown = angsep_list[ai + 1]
            marker_label.append(f'{aup}~{adown}')
        
        print(marker_label)
        
        
        return twpol_list, marker_label
    
    
    
    def pol_bin_plot(self, pol_list):
        
        vessel = self.data['vessel']
        RadLoc = self.data['grid']['RadLoc']
        VertLoc = self.data['grid']['VertLoc']
        
        
        "shell"
        
        north_R = RadLoc[-1, :]
        north_Z = VertLoc[-1, :]
        
        south_R_core = RadLoc[0, 25:72]
        south_Z_core = VertLoc[0, 25:72]
        
        south_R_pfr1 = RadLoc[0, 73:]
        south_Z_pfr1 = VertLoc[0, 73:]
        
        south_R_pfr2 = RadLoc[0, :24]
        south_Z_pfr2 = VertLoc[0, :24]
        
        sep_R_up = RadLoc[18, :]
        sep_Z_up = VertLoc[18, :]
        
        sep_R_down = RadLoc[17, :]
        sep_Z_down = VertLoc[17, :]
        
        sep_R = 0.5*(sep_R_up + sep_R_down)
        sep_Z = 0.5*(sep_Z_up + sep_Z_down)
        
        twpol_list, marker_label = self.bin_check(pol_list = pol_list)
        
        st = twpol_list[2][0]
        ed = twpol_list[2][-1]
        
        
        st_R = RadLoc[:, st]
        st_Z = VertLoc[:, st]
        
        ed_R = RadLoc[:, ed]
        ed_Z = VertLoc[:, ed]
        
        
        outl = twpol_list[-1][-1]
        innl = twpol_list[0][0]
        
        outl_R = RadLoc[:, outl]
        outl_Z = VertLoc[:, outl]
        
        
        innl_R = RadLoc[:, innl]
        innl_Z = VertLoc[:, innl]
        
        
        
        
        fig, axs = plt.subplots()
        
        "plot the shell"
        
        axs.plot(north_R, north_Z,'-', color = 'blue', label= 'B2 contour')
        axs.plot(south_R_pfr1, south_Z_pfr1,'-', color = 'blue')
        axs.plot(south_R_pfr2, south_Z_pfr2,'-', color = 'blue')
        axs.plot(south_R_core, south_Z_core,'-', color = 'blue')
        axs.plot(st_R, st_Z,'-', color = 'red', label = marker_label[2])
        axs.plot(ed_R, ed_Z,'-', color = 'red')
        axs.plot(sep_R, sep_Z,'-', color = 'green', label = 'separatrix')
        axs.plot(outl_R, outl_Z,'--', color = 'black', label = 'plot limits')
        axs.plot(innl_R, innl_Z,'--', color = 'black')
        axs.legend(loc= 'upper right')

        axs.plot(vessel[:,0]/1000, vessel[:,1]/1000, color = 'black')
        axs.set_xlabel("R: [m]")
        axs.set_ylabel("Z: [m]")
        
        plt.gca().set_aspect('equal')

