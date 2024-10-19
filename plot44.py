# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 21:19:03 2024

@author: ychuang
"""

from SOLPSplotter_iout_dataprocess import iout_process
from matplotlib.offsetbox import AnchoredText
import load_B2_data_method as lBdm
import numpy as np
import matplotlib.pyplot as plt


class plot_ft44(iout_process):
    def __init__(self, DefaultSettings, loadDS):
        iout_process.__init__(self, DefaultSettings, loadDS)




    def polplotft44(self, pol_list):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                neutemp = self.data['ft44'][aa]['tab2'][1:97, 1:37]
                
                ev = 1.6021766339999999 * pow(10, -19)
                neutemp_pro = neutemp / ev
                
                ang_list = self.data['angle']['angle_list'][aa]
                # print(np.shape(nadiff))
                st = int(pol_list[0])
                ed = int(pol_list[-1]) + 1
                
                
                
                weight = self.data['polpsi_weight'][aa]
                
                
                Tn_list = []
                
                for ii, kk in enumerate(pol_list):
                    
                    Tn = neutemp_pro[int(kk), 20]*weight[ii] + (1 - weight[ii])*neutemp_pro[int(kk), 18]
                    Tn_list.append(Tn)
                
                
                axs.plot(ang_list, Tn_list, color = color_dic[aa], 
                              label = 'A = {}'.format(A_dic[aa]))
                
                
                
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                # axs.plot(ang_list, neutemp_pro[st:ed, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
            
            axs.legend(loc= 'best')
            axs.set_title('atomic temperature [eV]')
            axs.set_xlabel('poloidal angle')
    
    
    
    def polplotft44_leg(self):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                neutemp = self.data['ft44'][aa]['tab2'][1:97, 1:37]
                
                ev = 1.6021766339999999 * pow(10, -19)
                neutemp_pro = neutemp / ev
                
                index_list = np.linspace(1, 23, 23)
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs.plot(index_list, neutemp_pro[1:24, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
            
            axs.legend(loc= 'best')
            axs.set_title('atomic temperature inner leg')
            
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                neutemp = self.data['ft44'][aa]['tab2'][1:97, 1:37]
                
                ev = 1.6021766339999999 * pow(10, -19)
                neutemp_pro = neutemp / ev
                
                index_list = np.linspace(1, 22, 22)
                
                
                # axs.plot(ang_list, dndx[st:ed, 1])
                axs.plot(index_list, neutemp_pro[73:97, 18], color = color_dic[aa], label = 'A = {}'.format(A_dic[aa]))
            
            axs.legend(loc= 'best')
            axs.set_title('atomic temperature outer leg')
    
    
    
    def neuden_tar(self, side):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            # side = 'inner target'
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
                
                
                neutemp = self.data['ft44'][aa]['tab2'][1:97, 1:37]
                
                ev = 1.6021766339999999 * pow(10, -19)
                neutemp_pro = neutemp / ev
                
                if side == 'inner target':
                    psiN = self.data['psi']['psival'][aa][1:37, 0]
                    neuden_t = neuden[25, :]
                    
                    
                elif side == 'outer target':
                    
                    psiN = self.data['psi']['psival'][aa][1:37, -1]
                    neuden_t = neuden[60, :]

                
                    
                axs.plot(psiN, neuden_t, '-', color = color_dic[aa], 
                              label = 'A = {}'.format(A_dic[aa]))
                
            
            axs.legend(loc= 'best')
            axs.set_title('atomic density {}'.format(side))
            axs.set_xlabel('$\psi_N$')
    
    
    def neutemp_tar(self):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            side = 'inner target'
            
            fig, axs = plt.subplots()
                      
            for aa in self.data['dircomp']['multi_shift']:
                
                
                
                neutemp = self.data['ft44'][aa]['tab2'][:, :, 0]
                
                ev = 1.6021766339999999 * pow(10, -19)
                neutemp_pro = neutemp / ev
                
                if side == 'inner target':
                    psiN = self.data['psi']['psival'][aa][1:37, 0]
                    neuden_t = neutemp_pro[0, :]
                    
                    
                elif side == 'outer target':
                    
                    psiN = self.data['psi']['psival'][aa][1:37, -1]
                    neuden_t = neutemp_pro[-1, :]

                
                    
                axs.plot(psiN, neuden_t, '-', color = color_dic[aa], 
                              label = 'A = {}'.format(A_dic[aa]))
                
            
            axs.legend(loc= 'best')
            axs.set_title('atomic temperature inner leg')
            
            
            
    
    
    