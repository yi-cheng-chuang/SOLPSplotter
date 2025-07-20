# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 19:47:42 2025

@author: ychuang
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
from fit_data.fitting_method import fit_method_collection
from scipy.optimize import curve_fit
from lmfit import Model
from scipy.interpolate import UnivariateSpline



class Amal_Steven_mastuTS:
    
    def __init__(self, DF, data, fmc: fit_method_collection):
        
        self.DF = DF
        self.data = data
        self.fmc = fmc
    
        
    def Steven_LFS_cutoff(self, rad, ne, te, ne_er, te_er):
        
        rad_list = []
        ne_list = []
        te_list = []
        ne_er_list = []
        te_er_list = []
        
        for kk, it in enumerate(rad):
            
            if it >= 1.3 and np.isnan(ne[kk]) == False and np.isnan(te[kk]) == False:
                rad_list.append(it)
                ne_list.append(ne[kk])
                te_list.append(te[kk])
                ne_er_list.append(ne_er[kk])
                te_er_list.append(te_er[kk])
        
        rad_pro = np.array(rad_list)
        ne_pro = np.array(ne_list)
        te_pro = np.array(te_list)
        ne_er_pro = np.array(ne_er_list)
        te_er_pro = np.array(te_er_list)
        
        datac_dic = {'radius': rad_pro, 'ne': ne_pro, 'te': te_pro, 'ne_error': ne_er_pro,
                     'te_error': te_er_pro}
        
        return datac_dic
    
    
    
    
    def Amal_load_mastu_TS(self, plot_OD, plot_P, writefile):
        
        # Replace 'your_file.pkl' with the path to your pickle file
        gbase = self.data['dirdata']['gbase']
        rad_file_path = '{}/TS_data/51704_433ms_thomson_major_radius.txt'.format(gbase)
        den_file_path = '{}/TS_data/51704_433ms_thomson_density.txt'.format(gbase)
        tem_file_path = '{}/TS_data/51704_433ms_thomson_temperature.txt'.format(gbase)
        ne_er_file_path = '{}/TS_data/51704_433ms_thomson_density_error.txt'.format(gbase)
        te_er_file_path = '{}/TS_data/51704_433ms_thomson_temperature_error.txt'.format(gbase)
        


        radius = np.loadtxt(rad_file_path)
        ne = np.loadtxt(den_file_path)
        te = np.loadtxt(tem_file_path)
        ne_er = np.loadtxt(ne_er_file_path)
        te_er = np.loadtxt(te_er_file_path)


        datac_dic = self.Steven_LFS_cutoff(rad = radius, ne = ne, te = te, 
                                           ne_er = ne_er, te_er = te_er)
        
        self.data['ExpDict']['exp_dic'] = datac_dic
        
        func = self.data['midplane_calc']['midR_psi_func']
        
        RBS_func = self.data['gfile']['gcomp']['interp_dic']['RBS']
        
        
        rad_cut = datac_dic['radius']
        ne_cut = datac_dic['ne']
        te_cut = datac_dic['te']
        ne_er_cut = datac_dic['ne_error']
        te_er_cut = datac_dic['te_error']
        
        input_len = len(rad_cut)
        exp_psi = np.zeros(input_len)
        z_flat = 0.015*np.ones(input_len)
        
        print('test np ones function')
        print(input_len)
        # print(z_flat)
        
        for i in range(input_len):
            
            exp_psi[i] = RBS_func(rad_cut[i], z_flat[i])
        
        
        psi_mid = RBS_func(1.389, 0.015)
        # print(psi_mid)
        # print(exp_psi)
        # print('separatrix is {:.2f}'.format(psi_mid))
        # exp_psi = func(dlcfs_flat)
        
        n_tot = 200
        psi_pro = np.linspace(exp_psi.min(), exp_psi.max(), num= n_tot, dtype= float)
        
        shift_psi = -0.04
        
        new_psi = exp_psi + shift_psi*np.ones(input_len)
        
        # print(new_psi)
        
        

        Ne = [x*pow(10, -19) for x in ne_cut]
        Te = [x*pow(10, -3) for x in te_cut]
        
        p0 = [1, 0.6, 0.01, 0.01, 3/14]
        p1 = [1, 0.6, 0.01, 0.01, 3/14]

        popt_ne, pcov_ne = curve_fit(self.fmc.tanh, new_psi, Ne, p0)      
        popt_te, pcov_te = curve_fit(self.fmc.tanh, new_psi, Te, p1)
        
        shift_psipro = np.linspace(new_psi.min(), new_psi.max(), num= n_tot, dtype= float)
        
          
        ne_fit = self.fmc.tanh(shift_psipro, popt_ne[0], popt_ne[1], popt_ne[2], 
                         popt_ne[3], popt_ne[4])*pow(10, 19)
        te_fit = self.fmc.tanh(shift_psipro, popt_te[0], popt_te[1], popt_te[2], 
                         popt_te[3], popt_te[4])*pow(10, 3)
        
        
        ne_dat = [x*pow(10, -20) for x in ne_fit]
        te_dat = [x*pow(10, -3) for x in te_fit]
        
        # exp_dic = {'exp_psi': exp_psi, 'exp_ne': ne_cut, 'exp_te': te_cut}
        
        
        fitprofile = {'exp_psi': exp_psi, 'fit_ne': ne_fit, 
                      'fit_te': te_fit}
        
        self.data['ExpDict']['exp_dic']['exp_psi'] = exp_psi
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
            filename = 'fit_{}.dat'.format(self.data['dircomp']['Shot'])
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
            plt.scatter(rad_cut, ne_cut, label='ne', color='b', linewidth=2)

            # Add labels and title
            plt.xlabel("radius (m)", fontsize=12)
            plt.title("ne [$m^{-3}$] vs. radius", fontsize=14)
            plt.legend()

            # Create the plot
            plt.figure()
            plt.errorbar(new_psi, ne_cut, yerr= ne_er_cut, label='ne', color='b', fmt = 'o')
            plt.plot(shift_psipro, ne_fit, label='ne_fit', color='g', linewidth=2)

            # Add labels and title
            plt.xlabel("psiN", fontsize=12)
            plt.title("ne [$m^{-3}$] vs. psiN", fontsize=14)
            plt.legend()
            
            
            # Create the plot
            plt.figure()
            plt.scatter(rad_cut, te_cut, label='te', color='b', linewidth=2)

            # Add labels and title
            plt.xlabel("radius (m)", fontsize=12)
            plt.title("te [eV] vs. radius", fontsize=14)
            plt.legend()


            # Create the plot
            plt.figure()
            plt.errorbar(new_psi, te_cut, yerr= te_er_cut, label='te', color='b', fmt = 'o')
            plt.plot(shift_psipro, te_fit, label='te_fit', color='g', linewidth=2)

            # Add labels and title
            plt.xlabel("psiN", fontsize=12)
            plt.title("te [eV] vs. psiN", fontsize=14)


            # Add legend
            plt.legend()

            # Show the plot
            plt.show()
        
    
"""   
    
    def Amal_load_mastu_TS(self, plot_OD, plot_P, writefile):
        
        # Replace 'your_file.pkl' with the path to your pickle file
        gbase = self.data['dirdata']['gbase']
        file_path = '{}/TS_pedestal_49404.pkl'.format(gbase)


        radius = np.loadtxt('filename.txt')  # Replace with your actual file path
        print(data)


        
        
        
        
        
        
        
        
        
        

"""
