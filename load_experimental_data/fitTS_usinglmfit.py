# -*- coding: utf-8 -*-
"""
Created on Tue May  6 19:16:18 2025

@author: ychuang
"""

import pickle
from lmfit import Model
import numpy as np
import matplotlib.pyplot as plt
from fit_data.fitting_method import fit_method_collection
from load_experimental_data.SOLPSplotter_mastuTS import preprocess_mastuTS



class lmfit_TS:
    
    def __init__(self, DF, data, fmc: fit_method_collection, pm: preprocess_mastuTS):
        
        self.DF = DF
        self.data = data
        self.fmc = fmc
        self.pm = pm
    
    
    def lmfit_tanh(self, TS_data, psi, plot):
        
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
        # print(result_Te.fit_report())

        #Plot
        
        if plot:
            
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
            
            
            
            n_tot = 200
            psi_pro = np.linspace(psi.min(), psi.max(), num= n_tot, dtype= float)
            
            
            
            # Evaluate the fitted model at x_new
            y_new = result_Te.eval(r = psi_pro)
            
            # Plot the new prediction
            plt.figure()
            plt.plot(psi, Te, 'bo', label='data')
            plt.plot(psi_pro, y_new, 'g--', label='fit (extrapolated)')
            plt.legend()
            plt.show()
            
        
        
    
    
    
    def load_mastuTS_test(self, plot):
        
        # Replace 'your_file.pkl' with the path to your pickle file
        gbase = self.data['dirdata']['gbase']
        file_path = '{}/TS_pedestal_49404.pkl'.format(gbase)

        # Open and load the pickle file
        with open(file_path, 'rb') as file:
            data = pickle.load(file)
        
        
        print(type(data))
        print(data[49404].keys())
        print(data[49404]['laser_time'])
        print(data[49404]['Te_arrays']['Pulse_0.799_0.830'].keys())
        print(data[49404]['Te_arrays']['Pulse_0.799_0.830']['lasers'].keys())
        time_list = list(data[49404]['Te_arrays']['Pulse_0.799_0.830']['lasers'].keys())

        print(data[49404]['Te_arrays']['Pulse_0.799_0.830']['lasers']['0.80614'].keys())
        dlcfs = data[49404]['Te_arrays']['Pulse_0.799_0.830']['lasers']['0.80614']['dlcfs']
        # print(dlcfs)

        # radius = target_dic['radius']
        # dlcfs = target_dic['dlcfs']
        # te = target_dic['te']
        # ne = target_dic['ne']
        
        
        alldata_dic = {}

        for time in time_list:
            
            target_dic = data[49404]['Te_arrays']['Pulse_0.799_0.830']['lasers'][time]
            # radius = target_dic['radius']
            radius = target_dic['radius']
            te = target_dic['te']
            ne = target_dic['ne']
            # print(type(ne))
            datac_dic = self.pm.LFS_cutoff(dlcfs = radius, ne = ne, te = te)
            
            alldata_dic[time] = datac_dic
                    

        "try to turn into array"
        
        
        func = self.data['midplane_calc']['midR_psi_func']
        
        RBS_func = self.data['gfile']['gcomp']['interp_dic']['RBS']
        
        
        
        for tt in time_list:
            
            
            radius = alldata_dic[tt]['dlcfs']
            
            
            input_len = len(radius)
            exp_psi = np.zeros(input_len)
            z_flat = 0.015*np.ones(input_len)
            
            print('test np ones function')
            print(input_len)
            # print(z_flat)
            
            for i in range(input_len):
                
                exp_psi[i] = RBS_func(radius[i], z_flat[i])
        
        
            self.lmfit_tanh(TS_data = alldata_dic[tt], psi = exp_psi, plot = plot)
        
        


