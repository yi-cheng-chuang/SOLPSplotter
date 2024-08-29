# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 14:24:08 2024

@author: ychuang
"""


from SOLPSplotter_fit import profile_fit
import matplotlib.pyplot as plt


class SOLPS_converge(profile_fit):
    
    def __init__(self, DefaultSettings, loadDS):
        profile_fit.__init__(self, DefaultSettings, loadDS)
            
    
    def set_plot(self, publish):
        
        plt.rcParams.update({'font.weight': 'normal'})
        plt.rc('lines', linewidth= 3, markersize= 5)
        plt.rcParams.update({'font.size': 14})
        plt.rcParams.update({'figure.facecolor':'w'})
        plt.rcParams.update({'mathtext.default': 'regular'})
    
    
    def profile_combine(self):
        
        exp_psiN = self.data['ExpDict']['psi_normal']
        exp_ne = self.data['ExpDict']['electron_density(10^20/m^3)']
        exp_te = self.data['ExpDict']['electron_temperature(KeV)']
        
        fit_psiN = self.data['experimental_fit']['psiN']
        fit_ne = self.data['experimental_fit']['ne']
        fit_te = self.data['experimental_fit']['te']
        
        
        exp_sol_psiN = []
        exp_sol_ne = []
        exp_sol_te = []
        for ind, epN in enumerate(exp_psiN):
            if epN >= 1:
                exp_sol_psiN.append(epN)
                exp_sol_ne.append(exp_ne[ind])
                exp_sol_te.append(exp_te[ind])
        
        
        fit_edge_psiN = []
        fit_edge_ne = []
        fit_edge_te = []
        for ind, feN in enumerate(fit_psiN):
            if epN <= 1:
                fit_edge_psiN.append(epN)
                fit_edge_ne.append(exp_ne[ind])
                fit_edge_te.append(exp_te[ind])
        
        print('this is the exp sol psiN')
        print(exp_sol_psiN)
        print('this is the fit edge psiN')
        print(fit_edge_psiN)
        