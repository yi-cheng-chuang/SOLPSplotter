# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 21:54:40 2023

@author: user
"""

import matplotlib.pyplot as plt
from SOLPS_input.input_setting import loadDS_dic
from load_experimental_data.load_expdata_method import read_expdata_method
# from load_simulation_data.load_B2_data_method import 
from fit_data.fitting_method import fit_method_collection
from load_coordinate.SOLPSplotter_geo import load_geometry
from scipy.optimize import curve_fit
import numpy as np



    
class load_expdata:
    
    
    def __init__(self, DF, data, fmc: fit_method_collection, rem: read_expdata_method, lg: load_geometry):
        
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
            
    
    def check_and_loadpsi1D(self, itername):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == False:
            
            self.lg.check_b2mn(itername = None)
            jxa = self.data['b2mn']['jxa']
            self.lg.calcpsi_1D(pol_loc= str(jxa), no_coord_avg_check = False)
            psi_solps = self.data['psi']['psi_{}_val'.format(str(jxa))]
            
            return psi_solps
        
        elif withshift == True and withseries == False:
            self.lg.check_b2mn(itername = itername)
            jxa = self.data['b2mn'][itername]['jxa']
            self.lg.calcpsi_1D(pol_loc= str(jxa), no_coord_avg_check = False)
            psi_solps = self.data['psi']['psi_{}_val'.format(str(jxa))][itername]
            
            return psi_solps
        
        elif withshift == False and withseries == True:
            
            self.lg.check_b2mn(itername = itername)
            jxa = self.data['b2mn']['jxa']
            self.lg.calcpsi_1D(pol_loc= str(jxa), no_coord_avg_check = False)
            psi_solps = self.data['psi']['psi_{}_val'.format(str(jxa))]
            
            return psi_solps
        
        else:
            print('check_and_loadpsi1D function has a bug')
    
        
    def solpsgrid_data_store(self, x_coord, ne_fit_coe, te_fit_coe, plot_solps_fit):
        ne_fit_solps = self.fmc.tanh(x_coord, ne_fit_coe[0], ne_fit_coe[1], ne_fit_coe[2], ne_fit_coe[3], ne_fit_coe[4])
        te_fit_solps = self.fmc.tanh(x_coord, te_fit_coe[0], te_fit_coe[1], te_fit_coe[2], te_fit_coe[3], te_fit_coe[4])
        
        exp_fit_dic = {'psiN': x_coord, 'ne': ne_fit_solps, 'te': te_fit_solps,
                       'ne_coe': ne_fit_coe, 'te_coe': te_fit_coe}
        
        self.data['experimental_fit'] = exp_fit_dic
        
        if plot_solps_fit:
            "fit electron density in 38 grids"
            plt.figure(figsize=(7,7))
            plt.plot(x_coord, ne_fit_solps,'-o', color='r', label= 'electron density fit with shift')
             
            plt.xlabel('Magnetic flux coordinate: ${\psi_N}$')
            plt.ylabel('Electron density: ${n_e}$ (10$^{20}$*m$^{-3}$)')
            plt.title('Electron density')
            plt.legend()
            
            "fit electron temperature in 38 grids"
            plt.figure(figsize=(7,7))
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
        
        n_tot = 200
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == False:
            psi_solps = self.check_and_loadpsi1D(itername = None)
            
        elif withshift == True and withseries == False:
            psi_solps = self.check_and_loadpsi1D(itername = 'org')
        
        elif withshift == False and withseries == True:
            series_rap = list(self.data['dircomp']['Attempt'].keys())[0]
            psi_solps = self.check_and_loadpsi1D(itername = series_rap)
        
        else:
            print('fitmastexp function has a bug checking b2mn')
            
                  
        p0 = [0.97, 0.6, 0.01, 0.01, 3/14]
        p1 = [0.95, 0.2, 0.02, 0.01, 6/7]
        self.loadmastdata(EXP= True, fit= False)
        mast_dat_dict = self.data['ExpDict']
        psi = mast_dat_dict['psi_normal']
        ne = mast_dat_dict['electron_density(10^20/m^3)']
        ne_er = mast_dat_dict['density error(10^20/m^3)']
        te = mast_dat_dict['electron_temperature(KeV)']
        te_er = mast_dat_dict['temperature error(10^20/m^3)']
        
        popt_ne, pcov_ne = curve_fit(self.fmc.tanh, psi, ne, p0)      
        popt_te, pcov_te = curve_fit(self.fmc.tanh, psi, te, p1)

          
        x_model = np.linspace(min(psi), 1.1, n_tot)
        tanh_ne_fit = self.fmc.tanh(x_model, popt_ne[0], popt_ne[1], popt_ne[2], popt_ne[3], popt_ne[4])
        tanh_te_fit = self.fmc.tanh(x_model, popt_te[0], popt_te[1], popt_te[2], popt_te[3], popt_te[4])
        
        shift = 0
                
        sh_ne_fit = self.fmc.tanh(x_model, popt_ne[0] + shift, popt_ne[1], popt_ne[2], popt_ne[3], popt_ne[4])
        sh_te_fit = self.fmc.tanh(x_model, popt_te[0] + shift, popt_te[1], popt_te[2], popt_te[3], popt_te[4])
        
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
                
                      
        gnexp = np.gradient(tanh_ne_fit)
        dn = popt_ne[2]
        sym_pt = popt_ne[0]
        dtn = popt_te[2]
        te_sym_pt = popt_te[0]
        h = popt_ne[1]*pow(10, 20)
        
        ro_popt_te = np.round(popt_te, 2)
        
        sep_pos = ro_popt_te[0] - 0.5*np.log(2 - np.sqrt(3))*ro_popt_te[2]
        
        
        if plot_exp_and_fit:
            "experimental data and tanh fit"
            "electron density"
            plt.figure(figsize=(7,7))
            plt.plot(x_model, tanh_ne_fit, color='r', label= 'electron density fit')
            plt.errorbar(psi, ne, ne_er,fmt= "o", label= 'electron density experiment data')
            plt.axvline(x=dn + sym_pt, color= 'black',lw= 3, ls= '--', 
                        label= 'Pedestal width : $\Delta n_e$')
            plt.axvline(x=-dn + sym_pt, color= 'black',lw= 3, ls= '--')
            plt.axvline(x=dn + sym_pt, color= 'black',lw= 3, ls= '--', 
                        label= 'Pedestal width : $\Delta n_e$')
            plt.axvline(x=-dn + sym_pt, color= 'black',lw= 3, ls= '--')
            plt.xlabel('Magnetic flux coordinate: ${\psi_N}$')
            # plt.ylabel('Electron density: ${n_e}$ (10$^{20}$*m$^{-3}$)')
            plt.title('Electron density: ${n_e}$ (10$^{20}$m$^{-3}$)')
            plt.legend()
            
            
            "electron temperature"
            
            plt.figure(figsize=(7,7))
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
            
            plt.figure(figsize=(7,7))
            plt.plot(x_model, sh_ne_fit,'-o', color='r', label= 'electron density fit with shift')
            plt.plot(x_model, tanh_ne_fit,'-o', color='b', label= 'electron density fit')
            
            plt.xlabel('Magnetic flux coordinate: ${\psi_N}$')
            plt.ylabel('Electron density: ${n_e}$ (10$^{20}$*m$^{-3}$)')
            plt.title('Electron density')
            plt.legend()
            
            "electron tempurature"
            
            plt.figure(figsize=(7,7))
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
            self.solpsgrid_data_store(x_coord = psi_solps[:, 1], ne_fit_coe = sh_popt_ne, 
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
                fdir = '{}/{}/{}'.format(self.data['dirdata']['basedrt'], 
                                        self.DF.DEV, self.loadDS['fitfname'])
            else:
                print('exp fit file writing has a bug')
            
            for j in range(n_tot):
                w_list =[]
                w_list.append("{: .6f}".format(x_model[j]))
                w_list.append("{: .6f}".format(sh_ne_fit[j]))
                w_list.append("{: .6f}".format(sh_te_fit[j]))
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
        
        if data_print:
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
        elif data_print == False:
            pass
        else:
            print('data_print has a bug')
        
        
        
            
            
            
"""





"""  
        
        
        
        
        
        