# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 00:07:50 2025

@author: ychuang
"""


import matplotlib.pyplot as plt


class plot_setting():
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    
    def set_plot(self):
        
        if self.DF.plot_setting == 'mod_transcoe':
            
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 5, markersize= 9)
            plt.rcParams.update({'font.size': 16})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
        


