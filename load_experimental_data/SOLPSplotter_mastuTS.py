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
from scipy.interpolate import UnivariateSpline



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
            
            if it >= 1.33 and np.isnan(ne[kk]) == False and np.isnan(te[kk]) == False:
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
    
    
    
    def merge_sorted_lists(self, x_lists, y_lists):
        # Combine all x and y pairs
        combined = []
        for x_list, y_list in zip(x_lists, y_lists):
            combined.extend(zip(x_list, y_list))
        
        # Sort the combined list based on x values
        combined.sort(key=lambda pair: pair[0])
        
        # Unzip into final x and y lists
        final_x, final_y = zip(*combined)
        
        return list(final_x), list(final_y)
    
    
    
    
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


        
        dlcfs_flat = []
        ne_flat = []
        te_flat = []
        len_log = []

        for tm in time_list:
                         
            dlcfs_flat.append(alldata_dic[tm]['dlcfs'])
            ne_flat.append(alldata_dic[tm]['ne'])
            te_flat.append(alldata_dic[tm]['te'])
            len_log.append(len(alldata_dic[tm]['dlcfs']))
            
        
        radius_sort, ne_sort = self.merge_sorted_lists(x_lists = dlcfs_flat, y_lists = ne_flat)
        radius_sort, te_sort = self.merge_sorted_lists(x_lists = dlcfs_flat, y_lists = te_flat)
        
        # print(dlcfs_flat)
        # print(radius_sort)
        # print(ne_sort)
        print(len_log)
        
      
        
        "try to turn into array"


        # print(dlcfs_flat)

        plt.figure()

        plt.scatter(radius_sort, ne_sort, label= 'ne_flatter')

        # Add labels and title
        plt.xlabel("dlcfs (m)", fontsize = 12)
        plt.title("ne vs. dlcfs", fontsize = 14)

        # Add legend
        plt.legend()

        # Show the plot
        plt.show()
        
        
        func = self.data['midplane_calc']['midR_psi_func']
        
        RBS_func = self.data['gfile']['gcomp']['interp_dic']['RBS']
        
        
        input_len = len(radius_sort)
        exp_psi = np.zeros(input_len)
        z_flat = 0.015*np.ones(input_len)
        
        print('test np ones function')
        print(input_len)
        # print(z_flat)
        
        for i in range(input_len):
            
            exp_psi[i] = RBS_func(radius_sort[i], z_flat[i])
        
        
        psi_mid = RBS_func(1.389, 0.015)
        # print(psi_mid)
        # print(exp_psi)
        # print('separatrix is {:.2f}'.format(psi_mid))
        # exp_psi = func(dlcfs_flat)
        
        n_tot = 500
        psi_pro = np.linspace(exp_psi.min(), exp_psi.max(), num= n_tot, dtype= float)
        
        shift_psi = -0.03
        
        new_psi = exp_psi + shift_psi*np.ones(input_len)
        
        # print(new_psi)
        
        

        Ne = [x *pow(10, -19) for x in ne_sort]
        Te = [x*pow(10, -3) for x in te_sort]
        
        p0 = [1, 0.6, 0.01, 0.01, 3/14]
        p1 = [1, 0.6, 0.01, 0.01, 3/14]

        popt_ne, pcov_ne = curve_fit(self.fmc.tanh, new_psi, Ne, p0)      
        popt_te, pcov_te = curve_fit(self.fmc.tanh, new_psi, Te, p1)
        
        shift_psipro = np.linspace(new_psi.min(), new_psi.max(), num= n_tot, dtype= float)
        
          
        ne_fit = self.fmc.tanh(psi_pro, popt_ne[0], popt_ne[1], popt_ne[2], 
                         popt_ne[3], popt_ne[4])*pow(10, 19)
        te_fit = self.fmc.tanh(psi_pro, popt_te[0], popt_te[1], popt_te[2], 
                         popt_te[3], popt_te[4])*pow(10, 3)
        
        
        ne_dat = [x*pow(10, -20) for x in ne_fit]
        te_dat = [x*pow(10, -3) for x in te_fit]
        
        exp_dic = {'exp_psi': exp_psi, 'exp_ne': ne_sort, 'exp_te': te_sort}
        
        
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
            plt.scatter(exp_psi, ne_sort, label='ne', color='b', linewidth=2)
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
            plt.scatter(exp_psi, te_sort, label='te', color='b', linewidth=2)
            plt.plot(psi_pro, te_fit, label='te_fit', color='g', linewidth=2)

            # Add labels and title
            plt.xlabel("psiN", fontsize=12)
            plt.title("te [eV] vs. psiN", fontsize=14)


            # Add legend
            plt.legend()

            # Show the plot
            plt.show()
            
        
        
"""


ln = len(time_list)

an_list = []

for ti in time_list:
    
    an_list.append(len(alldata_dic[ti]['dlcfs']))
    

an = max(an_list)

print(an)




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
 
print(dlcfs_pro)
print(ne_pro)




spline = UnivariateSpline(radius_sort, te_sort, s= 100)  # s controls smoothing
y_smooth = spline(radius_sort)


# 2. Create a dense x-axis for smooth plotting
x_dense = np.linspace(min(radius_sort), max(radius_sort), 500)

# 3. Evaluate the spline
y_smooth = spline(x_dense)

# 4. Plot original data and spline
plt.figure()
plt.plot(radius_sort, te_sort, 'o', label='Original data')
plt.plot(x_dense, y_smooth, '-', label='Spline fit')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Spline Fit to Data')
plt.show()




"""          
        

    
    
    
    




