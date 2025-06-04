# -*- coding: utf-8 -*-
"""
Created on Thu May 29 17:22:33 2025

@author: ychuang
"""

import matplotlib.pyplot as plt
from SOLPS_input.input_setting import loadDS_dic
from load_experimental_data.load_expdata_method import read_expdata_method 
from fit_data.fitting_method import fit_method_collection
from load_coordinate.SOLPSplotter_geo import load_geometry
from scipy.optimize import curve_fit
import numpy as np



    
class load_piecewise_expdata:
    
    
    def __init__(self, DF, data, fmc: fit_method_collection, rem: read_expdata_method, 
                                     lg: load_geometry):
        
        self.DF = DF
        self.data = data
        self.fmc = fmc
        self.rem = rem
        self.lg = lg
        
        loadDS = loadDS_dic(DEV = self.DF.DEV)
        
        self.loadDS = loadDS
    
    
    
    def loadmastdata(self, EXP, fit):
        if EXP:
            # mastloc = '{}/{}/{}'.format(self.data['dirdata']['basedrt'], 
            #                        self.DEV, self.loadDS['expfilename'])
            
            
            terminal = self.DF.terminal
            
            if terminal == True:
                mastloc = '{}/{}'.format(self.data['dirdata']['topdrt'], 
                                        self.loadDS['expfilename'])
            
            
            elif terminal == False:
                mastloc = '{}/{}'.format(self.data['dirdata']['gbase'], 
                                        self.loadDS['expfilename'])
            
            else:
                print('loadmastdata has a bug')
                
            expdic = self.rem.read_mastfile(mastloc)
            self.data['ExpDict'] = expdic
            self.data['dirdata']['mastloc'] = mastloc
        
        if fit:
            fitloc = '{}/{}/{}'.format(self.data['dirdata']['basedrt'], 
                                    self.DEV, self.loadDS['fitfname'])
            fitdic = self.rem.read_fitfile(fitloc)
            self.data['fitprofile'] = fitdic
            self.data['dirdata']['fitloc'] = fitloc
            
    
    
    
        
    def solpsgrid_data_store(self, x_coord, ne_fit_coe, te_fit_coe, plot_solps_fit):
        
        
        
        
        ne_fit_solps = self.fmc.tanh(x_coord, ne_fit_coe[0], ne_fit_coe[1], ne_fit_coe[2], ne_fit_coe[3], ne_fit_coe[4])
        te_fit_solps = self.fmc.tanh(x_coord, te_fit_coe[0], te_fit_coe[1], te_fit_coe[2], te_fit_coe[3], te_fit_coe[4])
        
        exp_fit_dic = {'psiN': x_coord, 'ne': ne_fit_solps, 'te': te_fit_solps,
                       'ne_coe': ne_fit_coe, 'te_coe': te_fit_coe}
        
        self.data['experimental_fit'] = exp_fit_dic
        
        if plot_solps_fit:
            "fit electron density in 38 grids"
            plt.figure()
            plt.plot(x_coord, ne_fit_solps,'-o', color='r', label= 'electron density fit with shift')
             
            plt.xlabel('Magnetic flux coordinate: ${\psi_N}$')
            plt.ylabel('Electron density: ${n_e}$ (10$^{20}$*m$^{-3}$)')
            plt.title('Electron density')
            plt.legend()
            
            "fit electron temperature in 38 grids"
            plt.figure()
            plt.plot(x_coord, te_fit_solps,'-o', color='r', label= 'electron temperature fit with shift')
             
            plt.xlabel('Magnetic flux coordinate: ${\psi_N}$')
            plt.ylabel('Electron temperature: ${T_e}$ (KeV)')
            plt.title('Electron temperature')
            plt.legend()
            
            plt.show()
        elif plot_solps_fit == False:
            pass
        else:
            print('plot_solps_fit has a bug')
    
    
    
    
    def fitmastexp(self, plot_setting_dic):
        
        writefile = plot_setting_dic['writefile']
        plot_solps_fit = plot_setting_dic['plot_solps_fit']
        plot_exp_and_fit = plot_setting_dic['plot_exp_and_fit']
        plot_shift_compare = plot_setting_dic['plot_shift_compare']
        data_print = plot_setting_dic['data_print']
        piecewise = plot_setting_dic['piecewise']
        
        n_tot = 200
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == False:
            jxa = self.data['b2mn']['jxa']
            psi_solps = self.data['psi']['psival'][:, jxa]
            
        elif withshift == True and withseries == False:
            # jxa = self.data['b2mn']['org']['jxa']
            jxa = 56
            psi_solps = self.data['psi']['psival']['org'][:, jxa]
            
            core_psi = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                # jxa = self.data['b2mn'][aa]['jxa']
                jxa = 56
                psi_solps = self.data['psi']['psival'][aa][:, jxa]
                core_psi[aa] = round(psi_solps[0], 3)
        
        elif withshift == False and withseries == True:
            series_rap = list(self.data['dircomp']['Attempt'].keys())[0]
            jxa = self.data['b2mn'][series_rap]['jxa']
            psi_solps = self.data['psi']['psival'][series_rap][:, jxa]
        
        else:
            print('fitmastexp function has a bug checking b2mn')
            
        print('find psi_solps shape:')
        print(psi_solps.shape)
        
        
        p0 = [0.97, 0.6, 0.01, 0.01, 3/14]
        p0_1 = [1.05, 0.1, 0.05, 0.01]
        p1 = [0.95, 0.2, 0.02, 0.01, 6/7]
        self.loadmastdata(EXP= True, fit= False)
        mast_dat_dict = self.data['ExpDict']
        psi_org = mast_dat_dict['psi_normal']
        # Define your split point
        split_point = 0.55
        
        # Masks for the two regions
        mask = psi_org >= split_point
        
        psi = mast_dat_dict['psi_normal'][mask]
        ne = mast_dat_dict['electron_density(10^20/m^3)'][mask]
        ne_er = mast_dat_dict['density error(10^20/m^3)'][mask]
        te = mast_dat_dict['electron_temperature(KeV)'][mask]
        te_er = mast_dat_dict['temperature error(10^20/m^3)'][mask]
        
        
        if piecewise:
            # Define your split point
            split_point = 1
            
            # Masks for the two regions
            mask1 = psi < split_point
            mask2 = psi >= split_point
            
            
            popt_ne_ped, pcov_ne_ped = curve_fit(self.fmc.tanh, psi, ne, p0)
            
            popt_ne_SOL, pcov_ne_SOL = curve_fit(self.fmc.tanh0, psi[mask2], ne[mask2], p0_1)
        
        else:
            popt_ne, pcov_ne = curve_fit(self.fmc.tanh, psi, ne, p0)
        
              
        popt_te, pcov_te = curve_fit(self.fmc.tanh, psi, te, p1)

        if piecewise:
            
            # Define your split point
            split_point = 1
            # split_point_1 = 1.001
            
            # Masks for the two regions
            
            x_model = np.linspace(min(psi), 1.1, n_tot)
            
            
            mask3 = x_model < split_point
            mask4 = x_model >= split_point
            
            
            tanh_ne_fit = self.fmc.tanh(x_model, *popt_ne_ped)
            tanh_ne_fit_SOL = self.fmc.tanh0(x_model[mask4], *popt_ne_SOL)
        
        else:
            
            x_model = np.linspace(min(psi), 1.1, n_tot)
            tanh_ne_fit = self.fmc.tanh(x_model, popt_ne[0], popt_ne[1], popt_ne[2], popt_ne[3], popt_ne[4])
            
            

        tanh_te_fit = self.fmc.tanh(x_model, popt_te[0], popt_te[1], popt_te[2], popt_te[3], popt_te[4])
        
        shift = 0
        
        
        if piecewise:
            
            # Define your split point
            split_point = 1
            
            mask3 = x_model < split_point
            mask4 = x_model >= split_point
            
            sh_ne_fit_ped = self.fmc.tanh(x_model[mask3], popt_ne_ped[0] + shift, 
                            popt_ne_ped[1], popt_ne_ped[2], popt_ne_ped[3], popt_ne_ped[4])
            sh_ne_fit_SOL = self.fmc.tanh0(x_model[mask4], popt_ne_SOL[0] + shift, 
                            popt_ne_SOL[1], popt_ne_SOL[2], popt_ne_SOL[3])
            
        else:
            
            sh_ne_fit = self.fmc.tanh(x_model, popt_ne[0] + shift, popt_ne[1], popt_ne[2], popt_ne[3], popt_ne[4])
            
        
        # print(type(sh_ne_fit_ped))
        
            

        sh_te_fit = self.fmc.tanh(x_model, popt_te[0] + shift, popt_te[1], popt_te[2], popt_te[3], popt_te[4])
        
        
        if piecewise:
            pass
        
        else:
            
            coe_len = len(popt_ne) 
            
            sh_popt_ne = np.zeros(coe_len)
            sh_popt_te = np.zeros(coe_len)
            for i in range(len(popt_ne)):
                if i == 0:
                    sh_popt_ne[i] = popt_ne[i] + shift
                    sh_popt_te[i] = popt_te[i] + shift
                else:
                    sh_popt_ne[i] = popt_ne[i]
                    sh_popt_te[i] = popt_te[i]
            
        
        
        
        if piecewise:
            # dn = (popt_ne_ped[2] + popt_ne_SOL[2])/2
            # sym_pt = (popt_ne_ped[0] + popt_ne_SOL[2])/2
            dn = popt_ne_ped[2]
            sym_pt = popt_ne_ped[0]
        
        else:
            gnexp = np.gradient(tanh_ne_fit)
            dn = popt_ne[2]
            sym_pt = popt_ne[0]
            h = popt_ne[1]*pow(10, 20)
            
            

            
        dtn = popt_te[2]
        te_sym_pt = popt_te[0]
        ro_popt_te = np.round(popt_te, 2)        
        sep_pos = ro_popt_te[0] - 0.5*np.log(2 - np.sqrt(3))*ro_popt_te[2]
        
        if withshift:
            core_ne = {}
            core_te = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                
                
                if piecewise:
                    core_ne[aa] = '{:.3e}'.format(self.fmc.tanh(core_psi[aa], *popt_ne_ped)*pow(10, 20))
                
                else:
                    
                    core_ne[aa] = '{:.3e}'.format(self.fmc.tanh(core_psi[aa], *popt_ne)*pow(10, 20))
                    

                core_te[aa] = '{:.3e}'.format(self.fmc.tanh(core_psi[aa], *popt_te)*pow(10, 3))
        
        else:
            pass
        
        
        if plot_exp_and_fit:
            "experimental data and tanh fit"
            "electron density"
            
            plt.figure()
            if piecewise:
                
                # Define your split point
                split_point = 1
                
                mask3 = x_model < split_point
                mask4 = x_model >= split_point
                
                
                plt.plot(x_model[mask3], tanh_ne_fit[mask3], color='r', label= 'electron density ped fit')
                plt.plot(x_model[mask4], tanh_ne_fit[mask4], color='green', linestyle= '--')
                plt.plot(x_model[mask4], tanh_ne_fit_SOL, color='r', label= 'electron density SOL fit')
                
            else:                
                plt.plot(x_model, tanh_ne_fit, color='r', label= 'electron density fit')
                

            plt.errorbar(psi, ne, ne_er,fmt= "o", label= 'electron density experiment data')
            plt.axvline(x=dn + sym_pt, color= 'black',lw= 3, ls= '--', 
                        label= 'Pedestal width : $\Delta n_e$ from all range fit' )
            plt.axvline(x=-dn + sym_pt, color= 'black',lw= 3, ls= '--')
            plt.axvline(x=dn + sym_pt, color= 'black',lw= 3, ls= '--', 
                        label= 'Pedestal width : $\Delta n_e$ from all range fit')
            plt.axvline(x=-dn + sym_pt, color= 'black',lw= 3, ls= '--')
            plt.xlabel('Magnetic flux coordinate: ${\psi_N}$')
            # plt.ylabel('Electron density: ${n_e}$ (10$^{20}$*m$^{-3}$)')
            plt.title('Electron density: ${n_e}$ (10$^{20}$m$^{-3}$)')
            plt.legend()
            
            
            "electron temperature"
            
            plt.figure()
            plt.plot(x_model, tanh_te_fit, color='r', label= 'electron temperature fit')
            plt.errorbar(psi, te, te_er, fmt= "o", label= 'electron temperature experiment data')
            plt.axvline(x=dtn + te_sym_pt, color='black',lw=3, ls='--', 
                        label= 'temperature pedestal width : $\Delta T_e$')
            plt.axvline(x=-dtn + te_sym_pt, color='black',lw=3, ls= '--')
            plt.xlabel('Magnetic flux coordinate: ${\psi_N}$')
            # plt.ylabel('Electron temperature: ${T_e}$ (KeV)')
            plt.title('Electron temperature: ${T_e}$ (KeV)')
            plt.legend()
            
            plt.show()
        elif plot_exp_and_fit == False:
            pass
        else:
            print('plot_exp_and_fit has a bug')
        
        
        if plot_shift_compare:
            "fit profile and shift profile"
            "electron density"
            
            plt.figure()
            plt.plot(x_model, sh_ne_fit,'-o', color='r', label= 'electron density fit with shift')
            plt.plot(x_model, tanh_ne_fit,'-o', color='b', label= 'electron density fit')
            
            plt.xlabel('Magnetic flux coordinate: ${\psi_N}$')
            plt.ylabel('Electron density: ${n_e}$ (10$^{20}$*m$^{-3}$)')
            plt.title('Electron density')
            plt.legend()
            
            "electron tempurature"
            
            plt.figure()
            plt.plot(x_model, sh_te_fit,'-o', color='r', label= 'electron temperature fit with shift')
            plt.plot(x_model, tanh_te_fit,'-o', color='b', label= 'electron temperature fit')
            
            plt.xlabel('Magnetic flux coordinate: ${\psi_N}$')
            plt.ylabel('Electron temperature: ${T_e}$ (KeV)')
            plt.title('Electron temperature')
            plt.legend()
                   
            plt.show()
        elif plot_shift_compare == False:
            pass
        
        else:
            print('plot_shift_compare has a bug')
            

        try:
            self.solpsgrid_data_store(x_coord = psi_solps, ne_fit_coe = sh_popt_ne, 
                                      te_fit_coe = sh_popt_te, plot_solps_fit = plot_solps_fit)
        except:
            print('solpsgrid_data_store function has a bug')
        
            
        if writefile == True:
            w_datalist = []
            filename = 'fit_027205_275.dat'
            
            terminal = self.DF.terminal
            
            
            if terminal == True:
                fdir = '{}/{}'.format(self.data['dirdata']['topdrt'], self.loadDS['fitfname'])
                
            elif terminal == False:
                
                if piecewise:
                    
                    fdir = '{}/{}/{}'.format(self.data['dirdata']['basedrt'], 
                                            self.DF.DEV, 'pft_027205_275.dat')
                
                else:
                    
                    fdir = '{}/{}/{}'.format(self.data['dirdata']['basedrt'], 
                                            self.DF.DEV, self.loadDS['fitfname'])
                    
                

            else:
                print('exp fit file writing has a bug')
            
            
            
            if piecewise:
                
                import pandas as pd

                # # Full x array
                # x = x_model
                
                # y1 and y2 are piecewise results
                ne_ped = sh_ne_fit_ped  # corresponds to first 3 x values
                ne_SOL = sh_ne_fit_SOL  # corresponds to last 3 x values
                
                # Split x to match y1 and y2
                x1 = x_model[:len(ne_ped)]
                x2 = x_model[len(ne_ped):]
                
                # Combine using np.concatenate
                combined_x = np.concatenate([x1, x2])
                combined_ne = np.concatenate([ne_ped, ne_SOL])
                
                
                
                # Stack columns together
                data = np.column_stack((x_model, combined_ne, sh_te_fit))
                
                

                # Save to text file
                np.savetxt(fdir, data, fmt='%.6f', delimiter=' ')
                
                
                """
                
                segments = np.array(['part1'] * len(x1) + ['part2'] * len(x2))
                
                # Save to CSV
                df = pd.DataFrame({
                    'x': combined_x,
                    'y': combined_y,
                    'segment': segments
                })
                df.to_csv('piecewise_combined.csv', index=False)
                
                """

                
                
                
            else:
                
                for j in range(n_tot):
                    w_list =[]
                    w_list.append("{: .6f}".format(x_model[j]))
                    w_list.append("{: .6f}".format(sh_ne_fit[j]))
                    w_list.append("{: .6f}".format(sh_te_fit[j]))
                    w_writelist = ' '.join(str(y)+ "\t" for y in w_list)
                    w_datalist.append(w_writelist)
               
                with open(fdir, 'w') as f:
                    for l,w_line in enumerate(w_datalist):   
                        f.writelines(w_line + "\n")
            
            

            
            
            
            
            
        
        if data_print:
            if piecewise:
                
                print('the next line is popt_ne')
                print(popt_ne_ped)
            else:
                print('the next line is popt_ne')
                print(popt_ne)
                
            
            print('the next line is popt_te')
            print(popt_te)
            print('the next line is rounded popt_te')
            print(ro_popt_te)
            print('the next line is separatrix position')
            print(sep_pos)
            print('the next line is rounded separatrix position')
            print(round(sep_pos, 2))
            print('the next line is the temparature separatrix position calculation result')
            print(te_sym_pt + 0.5*np.log(2 + np.sqrt(3))*dtn + shift)
            
            
            if withshift:                
                print('ne core fixed value:')
                print(core_ne)
                
                print('te core fixed value:')
                print(core_te)
                
                print('core psi_N')
                print(core_psi)
                
           
        elif data_print == False:
            pass
        else:
            print('data_print has a bug')


    
    
    
        
        
        
"""





"""  
    
    
    
    
    
    