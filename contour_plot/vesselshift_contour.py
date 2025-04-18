# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 12:22:12 2025

@author: ychuang
"""




import matplotlib.pyplot as plt
from matplotlib import colors, cm, ticker
import matplotlib.tri as tri
from matplotlib.colors import LogNorm
import numpy as np
from numpy import ma
from matplotlib.offsetbox import AnchoredText




class shiftvessel_contour:
    
    
    
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
        
        

    def paper_vessel_method(self, vessel_data, shift_value, meter, 
                           color_dic, A_dic, itername, axs):
        
        if meter:
            axs.plot(vessel_data[:,0]/1000, vessel_data[:,1]/1000, 
        color = color_dic[itername], label= 'A = {}'.format(A_dic[itername]))
            
            
        else:
            axs.plot(vessel_data[:,0], vessel_data[:,1], 
    color = color_dic[itername], label= 'A = {}'.format(A_dic[itername]))
            
    
            
    def plot_vessel(self, itername, independent, meter):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == False:
            
            vessel = self.data['vessel']
            shift = self.data['dircomp']['shift_value']*1000
            
            self.plot_vessel_method(vessel_data = vessel, shift_value = shift, 
                                    independent = independent, meter= meter)
        
        elif withshift == True and withseries == False:
            
            vessel = self.data['vessel'][itername]
            shift = self.data['dircomp']['shift_dic'][itername]*1000
            
            self.plot_vessel_method(vessel_data = vessel, shift_value = shift, 
                                    independent = independent, meter = meter)
        
        elif withshift == False and withseries == True:
            
            vessel = self.data['vessel']
            shift = self.data['dircomp']['shift_value']
            
            self.plot_vessel_method(vessel_data = vessel, shift_value = shift, 
                                    independent = independent, meter = meter)
        
        else:
            
            print('plot_vessel function is not there yet!')
    
    
    
    def shaded_area(self):
        
            
        # CMAP = 'RdBu'
        CMAP = cm.viridis
        
        R_coord = self.data['grid']['RadLoc']['org']
        Z_coord = self.data['grid']['VertLoc']['org']
    
        shade = np.ones([38, 98])
        
        plt.contourf(R_coord, Z_coord, shade, levels= 1, cmap= 'Blues')
        
        
        # cs = ax.contourf(X, Y, z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
        # cbar = fig.colorbar(cs)

      
    
    
    def shift_vessel_in_one(self):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        
        if withshift == True and withseries == False:
            
            fig, axs = plt.subplots()
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            anchored_text = AnchoredText('{}'.format('vessel cross section'), loc='upper right')
            
            self.shaded_area()
            
            # ylabel_text = AnchoredText('{}'.format('Z [m]'), loc='center left')
            
            for aa in self.data['dircomp']['multi_shift']:
                
                vessel = self.data['vessel'][aa]
                shift = self.data['dircomp']['shift_dic'][aa]*1000
                
                self.paper_vessel_method(vessel_data = vessel, shift_value = shift,
            meter = True, color_dic = color_dic, itername = aa, axs = axs, A_dic = A_dic)
                
                axs.add_artist(anchored_text)
                # axs.add_artist(ylabel_text)
                axs.set_xlabel('R [m]')
                axs.set_ylabel('Z [m]')
                axs.legend(loc= 'center right', fontsize=10)
            
            
            axs.set_aspect('equal')
            
            # fig.savefig('vessel.pdf')