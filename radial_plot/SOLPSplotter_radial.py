# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 21:00:45 2024

@author: ychuang
"""

# from SOLPS_input.header import *

import matplotlib.pyplot as plt
import numpy as np


class radial_plot_opacity:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data

    
    
    def opacity_radial_method(self, result_dic, SEP, x_coord, Nd, Ne, Te, 
                                          P, log_flag):
        
        tanh_ne_fit = result_dic['tanh_ne_fit']
        tanh_te_fit = result_dic['tanh_te_fit']
        exp_an_fit = result_dic['exp_fit']
        dn = result_dic['pedestal_width_psiN'][0]
        dtn = result_dic['temperature_pedestal_width']
        efold = result_dic['efold_length_psiN'][0]
        opq = result_dic['dimensionless_opaqueness']
        xcoord_cut = result_dic['x_coord_cut']
        sym_pt = result_dic['ne_symmetry_point']
        te_sym_pt = result_dic['te_symmetry_point']
        
                      
        x = [-efold + max(xcoord_cut), max(xcoord_cut)]
        y = [min(exp_an_fit), min(exp_an_fit)]
        xd = [-dn + sym_pt, dn + sym_pt]
        yd = [tanh_ne_fit[SEP] , tanh_ne_fit[SEP]]
        xt = [-dtn + te_sym_pt, dtn + te_sym_pt]
        yt = [tanh_te_fit[SEP], tanh_te_fit[SEP]]
        
        
        
        
        plt.figure()
        if log_flag:
            plt.yscale('log')
        else:
            pass
        plt.plot(x_coord, Nd,'-', color = 'green', label= 'solps neutral density')
        # plt.plot(psi_RGI, Nd,'-', color = 'b', label= 'RGI_solps neutral density')
        plt.plot(xcoord_cut, exp_an_fit, color='r',lw= 5, ls='-', label= 'exponential fit')
        plt.axvline(x= max(xcoord_cut), color='orange',lw=3)
        plt.plot(x,y, color='orange', lw=3, ls='-', label= 'Neutral penetration length : $\lambda_{n_D}$')
        plt.axvline(x=-efold + max(xcoord_cut), color='orange',lw=3)
        plt.axvline(x= dn + sym_pt, color='black',lw=3, ls='--', 
                    label= 'fit range : $\Delta n_e$')
        plt.axvline(x= -dn + sym_pt, color='black',lw=3, ls='--')
        plt.xlabel('Normalized flux coordinate $\psi_N$')
        plt.title('Neutral density with fits')
        plt.legend()
            

        plt.figure()
        if log_flag:
            plt.yscale('log')
        else:
            pass
        # plt.plot(psi_xport, Ne,'-', color = 'r', label= 'solps_electron density')
        plt.plot(x_coord, Ne,'-', color = 'b', label= 'solps electron density')
        # plt.plot(fitdsa, cutNe,'-', color = 'g', label= 'experiment electron density')
        plt.plot(x_coord, tanh_ne_fit, ls='-', color='r',lw= 3, label= 'tanh fit')
        plt.plot(xd, yd, color='black', ls='-', lw=3, label= 'Pedestal width : $\Delta n_e$')
        plt.axvline(x=dn + sym_pt, color='black',lw=3)
        plt.axvline(x=-dn + sym_pt, color='black',lw=3)
        plt.axvline(x=max(xcoord_cut), color='orange',lw=3, ls='--', 
                    label = 'Neutral penetration length : $\lambda_{n_D}$')
        plt.axvline(x=-efold + max(xcoord_cut), color='orange',lw=3, ls='--')
        plt.xlabel('psiN')
        plt.title('Electron density with fits')
        plt.legend()
        
        
        plt.figure()
        if log_flag:
            plt.yscale('log')
        else:
            pass
        plt.plot(x_coord, Te, ls='-', color = 'b', label= 'solps electron tempurature')
        plt.plot(x_coord, tanh_te_fit, ls='-', color='r',lw= 3, label= 'tanh fit')
        plt.plot(xt, yt, color='black', ls='-', lw=3, label= 'temperature pedestal width [m]: $\Delta n_e$')
        plt.axvline(x=dtn + te_sym_pt, color='black', lw=3, ls='-')
        plt.axvline(x=-dtn + te_sym_pt, color='black', lw=3, ls='-')
        plt.xlabel('psiN')
        plt.title('Electron temperature with fits')
        plt.legend()

    

    
    def Opacity_study_radial_plot(self, pol_loc):
        
        
        dat_size = self.DF.data_size
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == False:
            
            pol_index = int(pol_loc[0])
            # Nd = self.data['radial_fit_data']['NeuDen'][:, pol_index]
            # Ne = self.data['radial_fit_data']['Ne'][:, pol_index]
            # Te = self.data['radial_fit_data']['Te'][:, pol_index]
            Nd = self.data['radial_fit_data']['NeuDen']
            Ne = self.data['radial_fit_data']['Ne']
            Te = self.data['radial_fit_data']['Te']
            SEP = int(self.data['DefaultSettings']['sep_index_dsa'])
            
            
            nx = self.data['b2fgeo']['nx']
            ny = self.data['b2fgeo']['ny']
            
            
            if dat_size == 'small':
                
                psi = self.data['psi']['psi_{}_val'.format(pol_loc[0])][:, 1][1:ny+1]
            
            else:
                
                print('I do not use full data size for my study')
                
            
            result_dic = self.data['radial_fit_data'] | self.data['opacity_poloidal']
            
            
            P = self.data['Parameter']
            self.opacity_radial_method(result_dic = result_dic, SEP = SEP, 
            x_coord = psi, Nd = Nd, Ne = Ne, Te = Te, P = P, log_flag = True)
        
        
        elif withshift == True and withseries == False:
            
            mix_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
            
                result_dic = self.data['radial_fit_data'][aa] | self.data['opacity_poloidal'][aa]
                mix_dic[aa] = result_dic
            
            self.data['mix_dic'] = mix_dic
            
            for aa in self.data['dircomp']['multi_shift']:
            
                pol_index = int(pol_loc[0])
                
                
                nx = self.data['b2fgeo'][aa]['nx']
                ny = self.data['b2fgeo'][aa]['ny']
                
                
                fstate = self.data['b2fstate'][aa]
                
                if dat_size == 'small':
                    
                    psi = self.data['psi']['psi_{}_val'.format(pol_loc[0])][aa][:, 1][1:ny+1]
                    data = self.data['ft44'][aa]['dab2']
                    Nd = np.transpose(data[:, :, 0])
                    Ne = fstate['ne'][1:nx+1, 1:ny+1].transpose()
                    Te_J = fstate['te'][1:nx+1, 1:ny+1].transpose()
                
                else:
                    
                    print('I do not use full data size for my study')
                
                
                ev = 1.6021766339999999 * pow(10, -19)
                Te = Te_J / ev
                
                Nd = Nd[:, pol_index]
                Ne = Ne[:, pol_index]
                Te = Te[:, pol_index]
                SEP = int(self.data['DefaultSettings']['sep_index_dsa'][aa])
                
                
                # result_dic = self.data['radial_fit_data'][aa] | self.data['opacity_poloidal'][aa]
                # mix_dic[aa] = result_dic
                
                P = self.data['Parameter']
                self.opacity_radial_method(result_dic = mix_dic[aa], SEP = SEP, 
                x_coord = psi, Nd = Nd, Ne = Ne, Te = Te, P = P, log_flag = False)
           
            
        
        elif withshift == False and withseries == True:
            
            mix_dic = {}
            
            for aa in self.data['dircomp']['Attempt'].keys():
            
                result_dic = self.data['radial_fit_data'][aa] | self.data['opacity_poloidal'][aa]
                mix_dic[aa] = result_dic
            
            self.data['mix_dic'] = mix_dic
            
            
            
            for aa in self.data['dircomp']['Attempt'].keys():
            
                pol_index = int(pol_loc[0])
                
                
                nx = self.data['b2fgeo'][aa]['nx']
                ny = self.data['b2fgeo'][aa]['ny']
                
                
                fstate = self.data['b2fstate'][aa]
                
                if dat_size == 'small':
                    
                    psi = self.data['psi']['psi_{}_val'.format(pol_loc[0])][aa][:, 2][1:ny+1]
                    data = self.data['ft44'][aa]['dab2']
                    Nd = np.transpose(data[:, :, 0])                    
                    Ne = fstate['ne'][1:nx+1, 1:ny+1].transpose()
                    Te_J = fstate['te'][1:nx+1, 1:ny+1].transpose()
                
                else:
                    
                    print('I do not use full data size for my study')
                
                
                ev = 1.6021766339999999 * pow(10, -19)
                Te = Te_J / ev
                
                Nd = Nd[:, pol_index]
                Ne = Ne[:, pol_index]
                Te = Te[:, pol_index]
                SEP = int(self.data['DefaultSettings']['sep_index_dsa'][aa])
                
                
                # result_dic = self.data['radial_fit_data'][aa] | self.data['opacity_poloidal'][aa]
                # mix_dic[aa] = result_dic
                
                P = self.data['Parameter']
                self.opacity_radial_method(result_dic = mix_dic[aa], SEP = SEP, 
                x_coord = psi, Nd = Nd, Ne = Ne, Te = Te, P = P, log_flag = False)
        

        
        
        elif withshift == True and withseries == True:
            print('Opacity_study_radial_plot_psi is not there yet, to be continue...')    
            
        else:
            print('Opacity_study_radial_plot_psi has a bug')      
        
        
"""      
        
    

            
            
"""




    
    
    