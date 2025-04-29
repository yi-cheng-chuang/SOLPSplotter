# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:07:41 2024

@author: ychuang
"""




import matplotlib.pyplot as plt
import numpy as np


class aspect_ratio_correlate:
    
    
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    
    
    def plot_data_reorder(self, pol_list):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        if withshift == True and withseries == False:
            
            result = self.data['opacity_poloidal']
            ll = len(self.data['dircomp']['Attempt'].keys())
            mm = len(pol_list)
            
            series_list = self.data['dircomp']['Attempt'].keys()
            
            # data_collect_opq = xr.DataArray(np.zeros((ll, mm)), 
            #                       coords=[series_list, pol_list], 
            #                 dims=['different_density','Poloidal_Location'], 
            #                  name = r'dimensionless opaqueness $m$')
            
            
            data_collect_opq = np.zeros((mm, ll))
            i = 0
            for ia, s_item in enumerate(self.data['dircomp']['multi_shift']):
                lb = np.asarray(result[s_item]['dimensionless_opaqueness'])
                data_collect_opq[:, ia] = lb
                # i = i + 1
            
            self.data['data_collect'] = data_collect_opq
            
            shift_list = np.zeros(ll)
            ka = 0
            for k, item in enumerate(self.data['dircomp']['multi_shift']):
                shift_list[k] = float(self.data['dircomp']['shift_dic'][item])
                # ka = ka + 1

            self.paper_data_reorder(iter_list = pol_list, change_var = shift_list,
                             data_collect = data_collect_opq)

        else:
            print('data reorder function is not there yet!')
            
        
    
    
    def data_reorder(self, iter_list, change_var, data_collect):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == True:
            plt.figure()
            color_list = ['red', 'salmon', 'orange', 'lime', 'green', 'darkgreen', 
                          'cyan', 'deepskyblue', 'navy', 'purple']
            
            
            series_flag = self.DF.series_flag
            
            
            if series_flag == 'change_den':
                for p in range(len(change_var)):
                    x_cor = np.ones(len(iter_list))*change_var[p]
                    plt.scatter(x_cor, data_collect[:, p], color= color_list[p], 
                                    label= 'efold length for different core electron density {}'.format(change_var[p]))
                plt.xlabel('Electron density: ${n_e}$ (m$^{-3}$)')
                plt.title('Experimental opaqueness verses shift core electron density')
                plt.legend()
                
            elif series_flag == 'eireneN':
                for p in range(len(change_var)):
                    x_cor = np.ones(len(iter_list))*change_var[p]
                    plt.scatter(x_cor, data_collect[:, p], color= color_list[p], 
                                    label= '{}'.format(change_var[p]))
                plt.xlabel('Eirene SOL particle number')
                # plt.yscale('log')
                plt.xscale('log')
                # plt.title('dimensionless opaqueness verses different SOL eirene particle number')
                plt.title('Experimental opaqueness')
                # plt.legend()
            
        elif withshift == True and withseries == False:
            
            
            plt.figure()
            
            A_list = [1.4, 2.0, 2.4, 2.8, 3.4]
            color_list = ['red', 'orange', 'green', 'blue', 'purple']
            
            
            for p in range(len(change_var)):
                x_cor = np.ones(len(iter_list))*A_list[p]
                plt.scatter(x_cor, data_collect[:, p], color= color_list[p], 
                            label= 'aspect ratio= {}'.format(str(A_list[p])))
            plt.axvline(x= 0.7/0.5, color='salmon',lw=3, ls='--', 
                        label= 'MAST aspect ratio')
            plt.axvline(x= 1.67/0.67, color='darkgreen',lw=3, ls='--', label= 'D3D aspect ratio')
            plt.axvline(x= 0.68/0.22, color='cyan',lw=3, ls='--', label= 'C-Mod/ ITER aspect ratio')
            # plt.axvline(x= 3.4, color='black',lw=3, ls='--', label= 'JT-60 aspect ratio')
            # plt.xlabel('aspect ratio')
            # plt.ylabel('dimensionless opaqueness')
            # plt.title('dimensionless opaqueness verses different modify distance')
            plt.title('Experimental opaqueness')
            plt.legend()
    
    
    
    def paper_data_reorder(self, iter_list, change_var, data_collect):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        if withshift == False and withseries == True:
            
            plt.figure()
            
            color_list = ['red', 'salmon', 'orange', 'lime', 'green', 'darkgreen', 
                          'cyan', 'deepskyblue', 'navy', 'purple']
            
            series_flag = self.DF.series_flag
            
            if series_flag == 'change_den':
                for p in range(len(change_var)):
                    x_cor = np.ones(len(iter_list))*change_var[p]
                    plt.scatter(x_cor, data_collect[:, p], color= color_list[p], 
                                    label= 'efold length for different core electron density {}'.format(change_var[p]))
                plt.xlabel('Electron density: ${n_e}$ (m$^{-3}$)')
                plt.title('Experimental opaqueness verses shift core electron density')
                plt.legend()
                
            elif series_flag == 'eireneN':
                for p in range(len(change_var)):
                    x_cor = np.ones(len(iter_list))*change_var[p]
                    plt.scatter(x_cor, data_collect[:, p], color= color_list[p], 
                                    label= '{}'.format(change_var[p]))
                plt.xlabel('Eirene SOL particle number')
                # plt.yscale('log')
                plt.xscale('log')
                # plt.title('dimensionless opaqueness verses different SOL eirene particle number')
                plt.title('Neutral opaqueness')
                # plt.legend()
            
        elif withshift == True and withseries == False:
            
            
            plt.figure()
            A_list = [1.4, 2.0, 2.4, 2.8, 3.4]
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            color_list = ['red', 'orange', 'green', 'blue', 'purple']
            
            
            for p in range(len(change_var)):
                x_cor = np.ones(len(iter_list))*A_list[p]
                plt.scatter(x_cor, data_collect[:, p], color= color_list[p], label= 'A= {}'.format(A_list[p]))
            plt.axvline(x= 0.7/0.5, color='salmon',lw=3, ls='--', label= 'MAST')
            plt.axvline(x= 1.67/0.67, color='lime',lw=3, ls='--', label= 'D3D')
            plt.axvline(x= 6.2/2, color='cyan',lw=3, ls='--', label= 'ITER')
            # plt.axvline(x= 3.4, color='black',lw=3, ls='--', label= 'JT-60 aspect ratio')
            plt.xlabel('aspect ratio')
            # plt.ylabel('dimensionless opaqueness')
            # plt.title('dimensionless opaqueness verses different modify distance')
            plt.title('Neutral opaqueness')
            plt.legend(loc= 'upper left')
    