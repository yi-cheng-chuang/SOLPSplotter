# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 21:54:40 2023

@author: user
"""
from SOLPSplotter_geo import load_geometry
import matplotlib.pyplot as plt
import load_mast_expdata_method as lmem
import load_B2_data_method as lbdm
import fitting_method as fm 
from scipy.optimize import curve_fit
import numpy as np


class load_expdata(load_geometry):
    
    def __init__(self, DEV, withshift, withseries, DefaultSettings, loadDS):
        load_geometry.__init__(self, DEV, withshift, 
                                        withseries, DefaultSettings)
        # Employee.__init__(self, first, last, pay)
        self.loadDS = loadDS
        
    def loadmastdata(self, EXP, fit):
        if EXP:
            # mastloc = '{}/{}/{}'.format(self.data['dirdata']['basedrt'], 
            #                         self.DEV, self.loadDS['expfilename'])
            mastloc = '{}/{}'.format(self.data['dirdata']['gbase'], 
                                    self.loadDS['expfilename'])
            expdic = lmem.read_mastfile(mastloc)
            self.data['ExpDict'] = expdic
            self.data['dirdata']['mastloc'] = mastloc
        
        if fit:
            fitloc = '{}/{}/{}'.format(self.data['dirdata']['basedrt'], 
                                    self.DEV, self.loadDS['fitfname'])
            fitdic = lmem.read_fitfile(fitloc)
            self.data['fitprofile'] = fitdic
            self.data['dirdata']['fitloc'] = fitloc
            
    
    def check_and_loadpsi1D(self, itername):
        if itername == None:
            self.check_b2mn(itername = None)
            jxa = self.data['b2mn']['jxa']
            self.calcpsi_1D(pol_loc= str(jxa), no_coord_avg_check = False)
            psi_solps = self.data['psi']['psi_{}_val'.format(str(jxa))]
            
            return psi_solps
        
        elif itername != None:
            self.check_b2mn(itername = itername)
            jxa = self.data['b2mn'][itername]['jxa']
            self.calcpsi_1D(pol_loc= str(jxa), no_coord_avg_check = False)
            psi_solps = self.data['psi']['psi_{}_val'.format(str(jxa))][itername]
            
            return psi_solps
        else:
            print('check_and_loadpsi1D function has a bug')
    
        
    def solpsgrid_data_store(self, x_coord, ne_fit_coe, te_fit_coe, plot_solps_fit):
        ne_fit_solps = fm.tanh(x_coord, ne_fit_coe[0], ne_fit_coe[1], ne_fit_coe[2], ne_fit_coe[3], ne_fit_coe[4])
        te_fit_solps = fm.tanh(x_coord, te_fit_coe[0], te_fit_coe[1], te_fit_coe[2], te_fit_coe[3], te_fit_coe[4])
        
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
        
        
        if self.withshift == False and self.withseries == False:
            psi_solps = self.check_and_loadpsi1D(itername = None)
            
        elif self.withshift == True and self.withseries == False:
            psi_solps = self.check_and_loadpsi1D(itername = 'org')
        
        elif self.withshift == False and self.withseries == True:
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
        
        popt_ne, pcov_ne = curve_fit(fm.tanh, psi, ne, p0)      
        popt_te, pcov_te = curve_fit(fm.tanh, psi, te, p1)

          
        x_model = np.linspace(min(psi), 1.1, n_tot)
        tanh_ne_fit = fm.tanh(x_model, popt_ne[0], popt_ne[1], popt_ne[2], popt_ne[3], popt_ne[4])
        tanh_te_fit = fm.tanh(x_model, popt_te[0], popt_te[1], popt_te[2], popt_te[3], popt_te[4])
        
        shift = 0
                
        sh_ne_fit = fm.tanh(x_model, popt_ne[0] + shift, popt_ne[1], popt_ne[2], popt_ne[3], popt_ne[4])
        sh_te_fit = fm.tanh(x_model, popt_te[0] + shift, popt_te[1], popt_te[2], popt_te[3], popt_te[4])
        
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
        
        ro_popt_te = np.round_(popt_te, 2)
        
        sep_pos = ro_popt_te[0] - 0.5*np.log(2 - np.sqrt(3))*ro_popt_te[2]
        
        
        if plot_exp_and_fit:
            "experimental data and tanh fit"
            "electron density"
            plt.figure(figsize=(7,7))
            plt.plot(x_model, tanh_ne_fit, color='r', label= 'electron density fit')
            plt.errorbar(psi, ne, ne_er,fmt= "o", label= 'electron density experiment data')
            plt.axvline(x=dn + sym_pt, color= 'black',lw= 3, ls= '--', 
                        label= 'Pedestal width [m]: $\Delta n_e$')
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
                        label= 'temperature pedestal width [m]: $\Delta T_e$')
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
            self.solpsgrid_data_store(x_coord = psi_solps[:, 2], ne_fit_coe = sh_popt_ne, 
                                      te_fit_coe = sh_popt_te, plot_solps_fit = plot_solps_fit)
        except:
            print('solpsgrid_data_store function has a bug')
        
            
        if writefile == True:
            w_datalist = []
            filename = 'wsh_027205_275.dat'
            fdir = '{}/{}/{}'.format(self.data['dirdata']['basedrt'], 
                                    self.DEV, self.loadDS['fitfname'])
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
if self.data['b2mn']['jxa'] == None:
    b2mn = lcm.scrape_b2mn(self.data['dirdata']['simudir']
                          + '/b2mn.dat')
    self.data['b2mn'] = b2mn
else:
    pass
jxa = self.data['b2mn']['jxa']
self.calcpsi_1D(pol_loc= str(jxa), no_coord_avg_check = False)
psi_solps = self.data['psi']['psi_{}_val'.format(str(jxa))]

# psi_sh = psi + shift

# sh_opt_ne, sh_cov_ne = curve_fit(fm.tanh, psi_sh, ne, p0)
# print(sh_opt_ne)
# sh_opt_te, sh_cov_te = curve_fit(fm.tanh, psi_sh, te, p1)
# print(sh_opt_te)
       
# x_sh = np.linspace(min(psi_sh), max(psi_sh), n_tot)

# def chunks(xs, n):
#     n = max(1, n)
#     return (xs[i:i+n] for i in range(0, len(xs), n))

# sep_test = np.array_split(range(11), 3)
# self.data['sep_test'] = sep_test
# print(sep_test[0])


BASEDRT = self.data['dirdata']['outputdir']['Output']
Attempt = self.data['dircomp']['Attempt']
# Attempts = len([self.data['dircomp']['a_shift']])
# N = len(self.data['dircomp']['a_shift'])
XGrid = int(self.data['b2fgeo']['nx'])
# X = self.data['gridsettings']['X']
XMin= 1
XMax= XGrid
# print(XGrid)
XDIM = int(self.data['DefaultSettings']['XDIM'])
# YSurf = int(self.data['b2fgeo']['ny'])
# Y = self.data['gridsettings']['Y']
YDIM = int(self.data['DefaultSettings']['YDIM'])
n = 0


# DRT = '{}/Attempt{}'.format(BASEDRT, str(Attempt))   #Generate path
test = param in self.Parameters.keys()
self.data['outputdata'][param] = np.zeros([YDIM, XDIM], dtype= np.float32)
# self.data['outputdata'][param] = xr.DataArray(np.zeros((YSurf,XGrid,N)), 
#                                  coords=[Y,X,Attempts], 
# dims=['Radial_Location','Poloidal_Location','Attempt'], name = param)
if test:
    # print('yes, {} is in parameter'.format(param))
    RawData = np.loadtxt('{}/{}{}'.format(BASEDRT, param, str(Attempt)),usecols = (3))
    # temp_dic = {param: RawData}
    # self.data['temp'] = temp_dic
elif test == False:
    print('no, {} is not in parameter'.format(param))
else:
    print('there might be a bug')
    
    
# def chunks(xs, n):
#     n = max(1, n)
#     return (xs[i:i+n] for i in range(0, len(xs), n))

# sep_test = np.array_split(range(11), 3)
# self.data['sep_test'] = sep_test
# print(sep_test[0])

if len(RawData) > 0:        
    if RawData.size == XDIM*YDIM:
        # self.data['outputdata'][param].values[:,:,n] = RawData.reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
        self.data['outputdata'][param] = RawData.reshape((YDIM,XDIM))
    elif RawData.size == XDIM*YDIM*2:
        raw_split = np.array_split(RawData, 2)
        param_dic = {'D_0': raw_split[0].reshape((YDIM,XDIM)), 
                     'D_1': raw_split[1].reshape((YDIM,XDIM))}
        self.data['outputdata'][param] = param_dic
        # print('let work on it')
        
        
    elif RawData.size != XDIM*YDIM:
        print('rawdata size is not equal to {}, it is {}'.format(str(XDIM*YDIM), str(RawData.size)))
    # elif RawData.size == XDIM*YDIM*2:
    #     self.data['outputdata'][param].values[:,:,n] = RawData.reshape((2*YDIM,XDIM))[1+YDIM:2*YDIM-1,XMin:XMax+1]
else:
    print('we have a problem loading rawdata')



BASEDRT = self.data['dirdata']['infolderdir'][aa]['outputdir']['Output']
Attempt = self.data['dircomp']['Attempt'][aa]
XGrid = int(self.data['b2fgeo'][aa]['nx'])
# print(XGrid)
XDIM = int(self.data['DefaultSettings']['XDIM'][aa])
YDIM = int(self.data['DefaultSettings']['YDIM'][aa])


# DRT = '{}/Attempt{}'.format(BASEDRT, str(Attempt))   #Generate path
test = param in self.Parameters.keys()
param_data_dic[aa] = np.zeros([YDIM, XDIM], dtype= np.float32)
# self.data['outputdata'][param] = xr.DataArray(np.zeros((YSurf,XGrid,N)), 
#                                  coords=[Y,X,Attempts], 
# dims=['Radial_Location','Poloidal_Location','Attempt'], name = param)
if test:
    # print('yes, {} is in parameter'.format(param))
    RawData = np.loadtxt('{}/{}{}'.format(BASEDRT, param, str(Attempt)),usecols = (3))
elif test == False:
    print('no, {} is not in parameter'.format(param))
else:
    print('there might be a bug')

if len(RawData) > 0:        
    if RawData.size == XDIM*YDIM:
        # self.data['outputdata'][param].values[:,:,n] = RawData.reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
        param_data_dic[aa] = RawData.reshape((YDIM,XDIM))
    elif RawData.size != XDIM*YDIM:
        print('rawdata size is not equal to {}'.format(str(XDIM*YDIM)))
    # elif RawData.size == XDIM*YDIM*2:
    #     self.data['outputdata'][param].values[:,:,n] = RawData.reshape((2*YDIM,XDIM))[1+YDIM:2*YDIM-1,XMin:XMax+1]
else:
    print('we have a problem loading rawdata')


BASEDRT = self.data['dirdata']['outputdir'][aa]['Output']
Attempt = self.data['dircomp']['Attempt'][aa]
XGrid = int(self.data['b2fgeo']['nx'])
# print(XGrid)
XDIM = int(self.data['DefaultSettings']['XDIM'])
YDIM = int(self.data['DefaultSettings']['YDIM'])


# DRT = '{}/Attempt{}'.format(BASEDRT, str(Attempt))   #Generate path
test = param in self.Parameters.keys()
param_data_dic[aa] = np.zeros([YDIM, XDIM], dtype= np.float32)
# self.data['outputdata'][param] = xr.DataArray(np.zeros((YSurf,XGrid,N)), 
#                                  coords=[Y,X,Attempts], 
# dims=['Radial_Location','Poloidal_Location','Attempt'], name = param)
if test:
    # print('yes, {} is in parameter'.format(param))
    RawData = np.loadtxt('{}/{}{}'.format(BASEDRT, param, str(Attempt)),usecols = (3))
elif test == False:
    print('no, {} is not in parameter'.format(param))
else:
    print('there might be a bug')

if len(RawData) > 0:        
    if RawData.size == XDIM*YDIM:
        # self.data['outputdata'][param].values[:,:,n] = RawData.reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
        param_data_dic[aa] = RawData.reshape((YDIM,XDIM))
    elif RawData.size != XDIM*YDIM:
        print('rawdata size is not equal to {}'.format(str(XDIM*YDIM)))
    # elif RawData.size == XDIM*YDIM*2:
    #     self.data['outputdata'][param].values[:,:,n] = RawData.reshape((2*YDIM,XDIM))[1+YDIM:2*YDIM-1,XMin:XMax+1]
else:
    print('we have a problem loading rawdata')


"""  
        
        
        
        
        
        