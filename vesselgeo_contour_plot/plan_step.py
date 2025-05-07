# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 01:32:17 2025

@author: ychuang
"""



import matplotlib.pyplot as plt 
import numpy as np


class step_planning:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data


    def mastu_Ra_measure(self):
        
        vessel = self.data['vessel']
        RadLoc = self.data['grid']['RadLoc']
        VertLoc = self.data['grid']['VertLoc']
        gfile = self.data['gfile']['g']
        
        
        "shell"
        
        north_R = RadLoc[-1, :]
        north_Z = VertLoc[-1, :]
        
       
        sep_R_up = RadLoc[18, :]
        sep_Z_up = VertLoc[18, :]
        
        sep_R_down = RadLoc[17, :]
        sep_Z_down = VertLoc[17, :]
        
        sep_R = 0.5*(sep_R_up + sep_R_down)
        sep_Z = 0.5*(sep_Z_up + sep_Z_down)
        
      
        
        fig, axs = plt.subplots()
        
        "plot the shell"
        
        axs.plot(RadLoc[-1, :47], VertLoc[-1, :47], 'b', label= 'B2 contour', linewidth=1)
        axs.plot(RadLoc[-1, 49:], VertLoc[-1, 49:], 'b', linewidth=1)
        axs.plot(gfile['rmaxis'], gfile['zmaxis'], 'r', markersize = 5)
        axs.plot(sep_R, sep_Z,'-', color = 'green', label = 'separatrix', linewidth = 2)
        axs.legend(loc= 'upper right')

        axs.plot(vessel[:,0]/1000, vessel[:,1]/1000, color = 'black',linewidth = 2)
        axs.set_xlabel("R: [m]")
        axs.set_ylabel("Z: [m]")
        
        plt.gca().set_aspect('equal')
    
    
    
    def vessel_change(self, plan_major_radius):
        
        gfile = self.data['gfile']['g']
        rmaxis = self.data['calc_minor_radius']['majR_loc']
        vessel = self.data['vessel']
        
        shift = plan_major_radius - rmaxis
        
        
        maxis_origin_R = vessel[:,0]/1000 - rmaxis
        
        
        
        fig, axs = plt.subplots()
        
        "plot the shell"
        

        axs.legend(loc= 'upper right')
        axs.plot(vessel[:,0]/1000, vessel[:,1]/1000, color = 'black',linewidth = 2, label = 'a = 0.54m')
        axs.plot(vessel[:,0]/1000 + shift, vessel[:,1]/1000, color = 'red',linewidth = 2, label = 'a = 0.54m')
        axs.plot(2*maxis_origin_R + plan_major_radius, 2*vessel[:,1]/1000, color = 'orange',linewidth = 2, label = 'a = 1.08m')
        axs.plot(3*maxis_origin_R + plan_major_radius, 3*vessel[:,1]/1000, color = 'green',linewidth = 2, label = 'a = 1.08m')
        axs.plot(3.7*maxis_origin_R + plan_major_radius, 3.7*vessel[:,1]/1000, color = 'blue',linewidth = 2)
        axs.set_xlabel("R: [m]")
        axs.set_ylabel("Z: [m]")
        axs.set_title('Machine size change plan')
        
            
        plt.gca().set_aspect('equal')
    
    
    
    
    def vessel_change(self, plan_major_radius, plan_minor_radius):
        
        gfile = self.data['gfile']['g']
        rmaxis = self.data['calc_minor_radius']['majR_loc']
        LFS_Rsep = self.data['calc_minor_radius']['LFS']['LFS_R_sep']
        HFS_Rsep = self.data['calc_minor_radius']['HFS']['HFS_R_sep']
        tot_a = self.data['calc_minor_radius']['tot_a']
        vessel = self.data['vessel']
        
        shift = plan_major_radius - rmaxis
        
        
        maxis_origin_R = vessel[:,0]/1000 - rmaxis
        LRsep_prim = LFS_Rsep - rmaxis
        HRsep_prim = HFS_Rsep - rmaxis
        
        
        step_a_ratio = plan_minor_radius/ (0.5*tot_a)
        ratio = round(step_a_ratio, 2)
        
        
        
        LRsep_step = ratio*LRsep_prim + plan_major_radius
        HRsep_step = ratio*HRsep_prim + plan_major_radius
        
        print('LFS Rsep is {:.2f}'.format(LRsep_step))
        print('HFS Rsep is {:.2f}'.format(HRsep_step))
        
        
        fig, axs = plt.subplots()
        
        "plot the shell"
        A1 = rmaxis/ (0.5*tot_a)
        A2 = 1.85/0.57
        A3 = 1.7/0.6
        A4 = plan_major_radius/ (1.5*tot_a)
        A5 = plan_major_radius/ plan_minor_radius
        
        
        ratio1 = round(plan_major_radius/ (A2* 0.5* tot_a), 2)
        a1 = round(plan_major_radius/ A2, 2)
        
        ratio2 = round(plan_major_radius/ (A3* 0.5* tot_a), 2)
        a2 = round(plan_major_radius/ A3, 2)
        

        axs.legend(loc= 'upper right')
        axs.plot(vessel[:,0]/1000, vessel[:,1]/1000, color = 'black',linewidth = 2, label = 'A = {:.2f}, a = 0.54m'.format(A1))
        axs.plot(ratio1*maxis_origin_R + plan_major_radius, ratio1*vessel[:,1]/1000, color = 'red',
                 linewidth = 2, label = 'A = {:.2f}, a = {:.2f}m, SPARC'.format(A2, a1))
        axs.plot(ratio2*maxis_origin_R + plan_major_radius, ratio2*vessel[:,1]/1000, 
                 color = 'orange',linewidth = 2, label = 'A = {:.2f}, a = {:.2f}m, DIII-D'.format(A3, a2))
        axs.plot(3*maxis_origin_R + plan_major_radius, 3*vessel[:,1]/1000, 
                 color = 'green',linewidth = 2, label = 'A = {:.2f}, a = 1.62m'.format(A4))
        axs.plot(ratio*maxis_origin_R + plan_major_radius, ratio*vessel[:,1]/1000, 
                 color = 'blue',linewidth = 2, label = 'A = {:.2f}, a = 2m, STEP'.format(A5))
        axs.set_xlabel("R: [m]")
        axs.set_ylabel("Z: [m]")
        axs.set_title('Machine size change plan')
        axs.legend()
        
            
        plt.gca().set_aspect('equal')
        
        
        
        
        
        