# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 23:57:42 2025

@author: ychuang
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
from fit_data.fitting_method import fit_method_collection
from scipy.optimize import curve_fit
from lmfit import Model



class preprocess_mastuTS:
    
    def __init__(self, DF, data, fmc: fit_method_collection):
        
        self.DF = DF
        self.data = data
        self.fmc = fmc
    
    
    
    
    def nan_filter(self, pos, ne, te):
        
        ne_output = []
        te_output = []
        pos_output = []
        
        
        for kk, item in enumerate(ne):
            if np.isnan(item) == False and np.isnan(te[kk]) == False:
                ne_output.append(item)
                te_output.append(te[kk])
                pos_output.append(pos[kk])
        
        data_dic = {'dlcfs': pos_output, 'ne': ne_output, 'te': te_output}
        
        return data_dic


    def LFS_cutoff(self, dlcfs, ne, te):
        
        dlcfs_list = []
        ne_list = []
        te_list = []
        
        for kk, it in enumerate(dlcfs):
            
            if it >= 1.25 and np.isnan(ne[kk]) == False and np.isnan(te[kk]) == False:
                dlcfs_list.append(it)
                ne_list.append(ne[kk])
                te_list.append(te[kk])
        
        dlcfs_pro = np.array(dlcfs_list)
        ne_pro = np.array(ne_list)
        te_pro = np.array(te_list)
        
        datac_dic = {'dlcfs': dlcfs_pro, 'ne': ne_pro, 'te': te_pro}
        
        return datac_dic


    def LFS_bin_cutoff(self, dlcfs, ne, te, ne_error, te_error):
        
        dlcfs_list = []
        ne_list = []
        te_list = []
        ne_error_list = []
        te_error_list = []
        
        for kk, it in enumerate(dlcfs):
            
            if it >= -0.2 and np.isnan(ne[kk]) == False and np.isnan(te[kk]) == False:
                dlcfs_list.append(it)
                ne_list.append(ne[kk])
                te_list.append(te[kk])
                ne_error_list.append(ne_error[kk])
                te_error_list.append(te_error[kk])
        
        dlcfs_pro = np.array(dlcfs_list)
        ne_pro = np.array(ne_list)
        te_pro = np.array(te_list)
        ne_error_pro = np.array(ne_error_list)
        te_error_pro = np.array(te_error_list)
        
        datac_dic = {'dlcfs': dlcfs_pro, 'ne': ne_pro, 'te': te_pro, 
                     'ne_error': ne_error_pro, 'te_error': te_error_pro}
        
        return datac_dic
    
    
    
    
    
    
    def load_mastu_TS(self, plot_OD, plot_P, writefile):
        
        # Replace 'your_file.pkl' with the path to your pickle file
        gbase = self.data['dirdata']['gbase']
        file_path = '{}/TS_pedestal_49404.pkl'.format(gbase)

        # Open and load the pickle file
        with open(file_path, 'rb') as file:
            data = pickle.load(file)

        # Display the loaded data
        # print(type(data))
        # print(data[49404].keys())
        # print(data[49404]['Te_arrays']['Pulse_0.820_0.840'].keys())
        # print(data[49404]['Te_arrays']['Pulse_0.820_0.840']['lasers']['0.82222'].keys())

        # target_dic = data[49404]['Te_arrays']['Pulse_0.820_0.840']['lasers']['0.82222']

        # print(target_dic['time'])
        
        
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
        
        # if plot_OD:
            
        #     # Create the plot
        #     plt.figure()
        #     plt.scatter(radius, ne, label='ne', color='b', linewidth=2)

        #     # Add labels and title
        #     plt.xlabel("radius (m)", fontsize=12)
        #     plt.title("ne vs. radius", fontsize=14)


        #     # Create the plot
        #     plt.figure()
        #     plt.scatter(dlcfs, ne, label='ne', color='b', linewidth=2)

        #     # Add labels and title
        #     plt.xlabel("dlcfs (m)", fontsize=12)
        #     plt.title("ne vs. dlcfs", fontsize=14)


        #     # Add legend
        #     plt.legend()

        #     # Show the plot
        #     plt.show()
        
        
        # dlcfs_list = []
        # ne_list = []
        # te_list = []
        
        # for kk, it in enumerate(dlcfs):
            
        #     if it >= -0.2 and np.isnan(ne[kk]) == False and np.isnan(te[kk]) == False:
        #         dlcfs_list.append(it)
        #         ne_list.append(ne[kk])
        #         te_list.append(te[kk])
        
        # dlcfs_pro = np.array(dlcfs_list)
        # ne_pro = np.array(ne_list)
        # te_pro = np.array(te_list)
        
        # print(dlcfs_pro)
        # print(ne_pro)
        
        
        
        alldata_dic = {}

        for time in time_list:
            
            target_dic = data[49404]['Te_arrays']['Pulse_0.799_0.830']['lasers'][time]
            # radius = target_dic['radius']
            radius = target_dic['radius']
            te = target_dic['te']
            ne = target_dic['ne']
            # print(type(ne))
            datac_dic = self.LFS_cutoff(dlcfs = radius, ne = ne, te = te)
            
            alldata_dic[time] = datac_dic
                    

        "prepare for flatten"

        ln = len(time_list)

        an_list = []

        for ti in time_list:
            
            an_list.append(len(alldata_dic[ti]['dlcfs']))
            

        an = max(an_list)

        print(an)

        dlcfs_flat = []
        ne_flat = []
        te_flat = []

        for tm in time_list:
            
            for kb, d in enumerate(alldata_dic[tm]['dlcfs']):
                
                dlcfs_flat.append(d)
                ne_flat.append(alldata_dic[tm]['ne'][kb])
                te_flat.append(alldata_dic[tm]['te'][kb])
            
            
        "try to turn into array"


        # print(dlcfs_flat)

        plt.figure()

        plt.scatter(dlcfs_flat, ne_flat, label= 'ne_flatter')

        # Add labels and title
        plt.xlabel("dlcfs (m)", fontsize=12)
        plt.title("ne vs. dlcfs", fontsize=14)

        # Add legend
        plt.legend()

        # Show the plot
        plt.show()
        
        
        func = self.data['midplane_calc']['midR_psi_func']
        
        RBS_func = self.data['gfile']['gcomp']['interp_dic']['RBS']
        
        
        input_len = len(dlcfs_flat)
        exp_psi = np.zeros(input_len)
        z_flat = 0.15*np.ones(input_len)
        
        print('test np ones function')
        # print(z_flat)
        
        for i in range(input_len):
            
            exp_psi[i] = RBS_func(dlcfs_flat[i], z_flat[i])
        
        
        psi_mid = RBS_func(1.37, 0.15)
        # print(psi_mid)
        # print(exp_psi)
        # print('separatrix is {:.2f}'.format(psi_mid))
        # exp_psi = func(dlcfs_flat)
        
        n_tot = 200
        psi_pro = np.linspace(exp_psi.min(), exp_psi.max(), num= n_tot, dtype= float)
        
        shift_psi = -0.09
        
        new_psi = exp_psi + shift_psi*np.ones(input_len)
        
        print(new_psi)
        
        

        Ne = [x *pow(10, -19) for x in ne_flat]
        Te = [x*pow(10, -3) for x in te_flat]
        
        p0 = [1, 0.6, 0.01, 0.01, 3/14]
        p1 = [1, 0.6, 0.01, 0.01, 3/14]

        popt_ne, pcov_ne = curve_fit(self.fmc.tanh, new_psi, Ne, p0)      
        popt_te, pcov_te = curve_fit(self.fmc.tanh, new_psi, Te, p1)
        
        shift_psipro = np.linspace(new_psi.min(), new_psi.max(), num= n_tot, dtype= float)
        
        # Create the model
        model = Model(self.fmc.tanh)
        
        # Provide initial parameter guesses
        params = model.make_params(r0 = 1,h = 0.6,d = 0.01,b = 0.01,m = 3/14)
        
        # Fit to the data
        result = model.fit(Te, params, r= exp_psi)
        
        # Show results
        print(result.fit_report())
        
        # Plot
        # plt.figure()
        # plt.plot(exp_psi, Te, 'bo', label='data')
        # plt.plot(exp_psi, result.best_fit, 'r-', label='fit')
        # plt.legend()
        
        
        
        
        # Evaluate the fitted model at x_new
        y_new = result.eval(r = psi_pro)
        
        # Plot the new prediction
        plt.figure()
        plt.plot(psi_pro, y_new, 'g--', label='fit (extrapolated)')
        plt.legend()
        plt.show()
        

          
        ne_fit = self.fmc.tanh(psi_pro, popt_ne[0], popt_ne[1], popt_ne[2], 
                         popt_ne[3], popt_ne[4])*pow(10, 19)
        te_fit = self.fmc.tanh(psi_pro, popt_te[0], popt_te[1], popt_te[2], 
                         popt_te[3], popt_te[4])*pow(10, 3)
        
        
        ne_dat = [x*pow(10, -20) for x in ne_fit]
        te_dat = [x*pow(10, -3) for x in te_fit]
        
        exp_dic = {'exp_psi': exp_psi, 'exp_ne': ne_flat, 'exp_te': te_flat}
        
        
        fitprofile = {'exp_psi': exp_psi, 'fit_ne': ne_fit, 
                      'fit_te': te_fit}
        
        self.data['ExpDict']['exp_dic'] = exp_dic
        self.data['ExpDict']['fitprofile'] = fitprofile
        
        te_sep = self.fmc.tanh(1, popt_te[0], popt_te[1], popt_te[2], 
                         popt_te[3], popt_te[4])*pow(10, 3)
        ne_sep = self.fmc.tanh(1, popt_ne[0], popt_ne[1], popt_ne[2], 
                         popt_ne[3], popt_ne[4])*pow(10, 19)
        psi_core = self.data['midplane_calc']['psi_solps_mid']
        
        neb = self.fmc.tanh(psi_core.min(), popt_ne[0], popt_ne[1], popt_ne[2], 
                         popt_ne[3], popt_ne[4])*pow(10, 19)
        teb = self.fmc.tanh(psi_core.min(), popt_te[0], popt_te[1], popt_te[2], 
                         popt_te[3], popt_te[4])*pow(10, 3)
        
        
        
        print('te sep is {:.2}'.format(te_sep))
        print('ne sep is {:.2}'.format(ne_sep))
        print('ne core boundary is {:.3}'.format(neb))
        print('te core boundary is {:.3}'.format(teb))
        
        
        
        if writefile:
            w_datalist = []
            filename = 'fit_49404_82.dat'
            d = self.data['dircomp']
            
            
            if self.DF.terminal == True:
                fdir = '{}/{}'.format(self.data['dirdata']['topdrt'], filename)
                
            elif self.DF.terminal == False:
                fdir = '{}/{}/{}/{}'.format(self.data['dirdata']['basedrt'], 
                                        self.DF.DEV, d['Shot'], filename)
            else:
                print('exp fit file writing has a bug')
            
            for j in range(n_tot):
                w_list =[]
                w_list.append("{: .6f}".format(psi_pro[j]))
                w_list.append("{: .6f}".format(ne_dat[j]))
                w_list.append("{: .6f}".format(te_dat[j]))
                w_writelist = ' '.join(str(y)+ "\t" for y in w_list)
                w_datalist.append(w_writelist)
            
            # for j in range(len(psi_solps[:, 2])):
            #     w_list =[]
            #     w_list.append("{: .6f}".format(psi_solps[:, 2][j]))
            #     w_list.append("{: .6f}".format(ne_fit_solps[j]))
            #     w_list.append("{: .6f}".format(te_fit_solps[j]))
            #     w_writelist = ' '.join(str(y)+ "\t" for y in w_list)
            #     w_datalist.append(w_writelist)
           
            with open(fdir, 'w') as f:
                for l,w_line in enumerate(w_datalist):   
                    f.writelines(w_line + "\n")
        
        
        
        
        
        if plot_P:
            
            # Create the plot
            plt.figure()
            plt.scatter(dlcfs, ne, label='ne', color='b', linewidth=2)

            # Add labels and title
            plt.xlabel("dlcfs (m)", fontsize=12)
            plt.title("ne [$m^{-3}$] vs. dlcfs", fontsize=14)
            plt.legend()

            # Create the plot
            plt.figure()
            plt.scatter(exp_psi, ne_flat, label='ne', color='b', linewidth=2)
            plt.plot(psi_pro, ne_fit, label='ne_fit', color='g', linewidth=2)

            # Add labels and title
            plt.xlabel("psiN", fontsize=12)
            plt.title("ne [$m^{-3}$] vs. psiN", fontsize=14)
            plt.legend()
            
            
            # Create the plot
            plt.figure()
            plt.scatter(dlcfs, te, label='te', color='b', linewidth=2)

            # Add labels and title
            plt.xlabel("dlcfs (m)", fontsize=12)
            plt.title("te [eV] vs. dlcfs", fontsize=14)
            plt.legend()


            # Create the plot
            plt.figure()
            plt.scatter(exp_psi, te_flat, label='te', color='b', linewidth=2)
            plt.plot(psi_pro, te_fit, label='te_fit', color='g', linewidth=2)

            # Add labels and title
            plt.xlabel("psiN", fontsize=12)
            plt.title("te [eV] vs. psiN", fontsize=14)


            # Add legend
            plt.legend()

            # Show the plot
            plt.show()
            
        
        
            
        

    
    
    
    




