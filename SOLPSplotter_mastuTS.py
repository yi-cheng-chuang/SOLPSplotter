# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 23:57:42 2025

@author: ychuang
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
from SOLPSplotter_PRmap import RP_mapping
import fitting_method as fm
from scipy.optimize import curve_fit



class preprocess_mastuTS(RP_mapping):
    
    def __init__(self, DefaultSettings, loadDS):
        RP_mapping.__init__(self, DefaultSettings, loadDS)
    
    
    def load_mastu_TS(self, plot_OD, plot_P, writefile):
        
        # Replace 'your_file.pkl' with the path to your pickle file
        gbase = self.data['dirdata']['gbase']
        file_path = '{}/TS_pedestal_49404.pkl'.format(gbase)

        # Open and load the pickle file
        with open(file_path, 'rb') as file:
            data = pickle.load(file)

        # Display the loaded data
        print(type(data))
        print(data[49404].keys())
        print(data[49404]['Te_arrays']['Pulse_0.820_0.840'].keys())
        print(data[49404]['Te_arrays']['Pulse_0.820_0.840']['lasers']['0.82222'].keys())

        target_dic = data[49404]['Te_arrays']['Pulse_0.820_0.840']['lasers']['0.82222']

        print(target_dic['time'])

        radius = target_dic['radius']
        dlcfs = target_dic['dlcfs']
        te = target_dic['te']
        ne = target_dic['ne']
        
        if plot_OD:
            
            # Create the plot
            plt.figure()
            plt.scatter(radius, ne, label='ne', color='b', linewidth=2)

            # Add labels and title
            plt.xlabel("radius (m)", fontsize=12)
            plt.title("ne vs. radius", fontsize=14)


            # Create the plot
            plt.figure()
            plt.scatter(dlcfs, ne, label='ne', color='b', linewidth=2)

            # Add labels and title
            plt.xlabel("dlcfs (m)", fontsize=12)
            plt.title("ne vs. dlcfs", fontsize=14)


            # Add legend
            plt.legend()

            # Show the plot
            plt.show()
        
        
        dlcfs_list = []
        ne_list = []
        te_list = []
        
        for kk, it in enumerate(dlcfs):
            
            if it >= -0.2 and np.isnan(ne[kk]) == False and np.isnan(te[kk]) == False:
                dlcfs_list.append(it)
                ne_list.append(ne[kk])
                te_list.append(te[kk])
        
        dlcfs_pro = np.array(dlcfs_list)
        ne_pro = np.array(ne_list)
        te_pro = np.array(te_list)
        
        # print(dlcfs_pro)
        # print(ne_pro)
        
        
        func = self.data['midplane_calc']['dsa_psi_func']
        
        exp_psi = func(dlcfs_list)
        
        n_tot = 200
        psi_pro = np.linspace(exp_psi.min(), exp_psi.max(), num= n_tot, dtype=float)
        
        Ne = [x *pow(10, -19) for x in ne_list]
        Te = [x*pow(10, -3) for x in te_list]
        
        p0 = [0.97, 0.6, 0.01, 0.01, 3/14]
        p1 = [0.95, 0.2, 0.02, 0.01, 6/7]
        
        popt_ne, pcov_ne = curve_fit(fm.tanh, exp_psi, Ne, p0)      
        popt_te, pcov_te = curve_fit(fm.tanh, exp_psi, Te, p1)

          
        ne_fit = fm.tanh(psi_pro, popt_ne[0], popt_ne[1], popt_ne[2], 
                         popt_ne[3], popt_ne[4])*pow(10, 19)
        te_fit = fm.tanh(psi_pro, popt_te[0], popt_te[1], popt_te[2], 
                         popt_te[3], popt_te[4])*pow(10, 3)

        
        ne_dat = [x *pow(10, -20) for x in ne_fit]
        te_dat = [x*pow(10, -3) for x in te_fit]
        
        exp_dic = {'exp_psi':exp_psi, 'exp_ne': ne_list, 'exp_te': te_list}
        
        
        
        
        
        fitprofile = {'exp_psi': exp_psi, 'fit_ne': ne_fit, 
                      'fit_te': te_fit}
        
        self.data['ExpDict']['exp_dic'] = exp_dic
        self.data['ExpDict']['fitprofile'] = fitprofile
        
        te_sep = fm.tanh(1, popt_te[0], popt_te[1], popt_te[2], 
                         popt_te[3], popt_te[4])*pow(10, 3)
        ne_sep = fm.tanh(1, popt_ne[0], popt_ne[1], popt_ne[2], 
                         popt_ne[3], popt_ne[4])*pow(10, 19)
        psi_core = self.data['midplane_calc']['psi_solps_mid']
        
        neb = fm.tanh(psi_core.min(), popt_ne[0], popt_ne[1], popt_ne[2], 
                         popt_ne[3], popt_ne[4])*pow(10, 19)
        teb = fm.tanh(psi_core.min(), popt_te[0], popt_te[1], popt_te[2], 
                         popt_te[3], popt_te[4])*pow(10, 3)
        print('te sep is {:.2}'.format(te_sep))
        print('ne sep is {:.2}'.format(ne_sep))
        print('ne core boundary is {:.3}'.format(neb))
        print('te core boundary is {:.3}'.format(teb))
        
        
        
        
        if writefile:
            w_datalist = []
            filename = 'fit_49404_82.dat'
            d = self.data['dircomp']
            
            
            if self.terminal == True:
                fdir = '{}/{}'.format(self.data['dirdata']['topdrt'], filename)
                
            elif self.terminal == False:
                fdir = '{}/{}/{}/{}'.format(self.data['dirdata']['basedrt'], 
                                        self.DEV, d['Shot'], filename)
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
            plt.scatter(exp_psi, ne_pro, label='ne', color='b', linewidth=2)
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
            plt.scatter(exp_psi, te_pro, label='te', color='b', linewidth=2)
            plt.plot(psi_pro, te_fit, label='te_fit', color='g', linewidth=2)

            # Add labels and title
            plt.xlabel("psiN", fontsize=12)
            plt.title("te [eV] vs. psiN", fontsize=14)


            # Add legend
            plt.legend()

            # Show the plot
            plt.show()
            
        
        
            
        

    
    
    
    




