# -*- coding: utf-8 -*-
"""
Created on Tue May  6 19:16:18 2025

@author: ychuang
"""

from lmfit import Model
import numpy as np
import matplotlib.pyplot as plt
from fit_data.fitting_method import fit_method_collection

class lmfit_TS:
    
    def __init__(self, DF, data, fmc: fit_method_collection):
        
        self.DF = DF
        self.data = data
        self.fmc = fmc
    
    
    def lmfit_tanh(self, TS_data, psi):
        
        ne_sort = TS_data['ne']
        te_sort = TS_data['te']
        
        
        Ne = [x *pow(10, -19) for x in ne_sort]
        Te = [x*pow(10, -3) for x in te_sort]
        

        # Create the model
        model = Model(self.fmc.tanh)

        # Provide initial parameter guesses
        params = model.make_params(r0 = 1,h = 0.6,d = 0.01,b = 0.01,m = 3/14)
        # params['b'].set(min=0.1, max=10)

        # Fit to the data
        result_Te = model.fit(Te, params, r= psi)
        result_Ne = model.fit(Ne, params, r= psi)

        # Show results
        print(result_Te.fit_report())

        #Plot


        plt.figure()
        plt.plot(psi, Te, 'bo', label='data')
        plt.plot(psi, result_Te.best_fit, 'r-', label='fit')
        plt.title('check out Te lmfit')
        plt.legend()


        plt.figure()
        plt.plot(psi, Ne, 'bo', label='data')
        plt.plot(psi, result_Ne.best_fit, 'r-', label='fit')
        plt.title('check out Ne lmfit')
        plt.legend()
        
        
        


