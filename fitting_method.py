# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 13:28:56 2023

@author: user
"""

import numpy as np
from scipy.optimize import curve_fit
import load_mast_expdata_method as lmem
import matplotlib.pyplot as plt

def tanh(r,r0,h,d,b,m):
    return b+(h/2)*(np.tanh((r0-r)/d)+1) + m*(r0-r-d)*np.heaviside(r0-r-d, 1)

def expfit(x,A,l):  #Removed vertical displacement variable B; seemed to cause 'overfitting'
    return A*np.exp(l*x)

def tanh_dsa(r,h,d,b,m):
    return b+(h/2)*(np.tanh((-r)/d)+1) + m*(-r-d)*np.heaviside(-r-d, 1)


def tanh_dsa_fit(dsa, ne, te):
    p0 = [0.4, 0.005, 0.05, 0.5]
    p1 = [2.2, 0.005, 0.05, 3.5]
    # mast_dat_dict = lmem.read_mastfile(mastloc)
    # psi = mast_dat_dict['psi_normal']
    # ne = mast_dat_dict['electron_density(10^20/m^3)']
    # te = mast_dat_dict['electron_temperature(KeV)']
    
    Ne = ne/ max(ne)
    Te = te/ max(te)
       
    popt_ne, pcov_ne = curve_fit(tanh_dsa, dsa, Ne, p0)
    print(popt_ne)
    popt_te, pcov_te = curve_fit(tanh_dsa, dsa, Te, p1)
    print(popt_te) 
    tanh_ne_fit = tanh_dsa(dsa, popt_ne[0], popt_ne[1], popt_ne[2], popt_ne[3])*max(ne)
    tanh_te_fit = tanh_dsa(dsa, popt_te[0], popt_te[1], popt_te[2], popt_te[3])*max(te)
    
    fit_tanh_dic = {'tanh_ne_fit': tanh_ne_fit, 'tanh_te_fit': tanh_te_fit, 
               'popt_ne': popt_ne, 'popt_te': popt_te}
    
    return fit_tanh_dic


def tanh_psi_fit(psi, ne, te):
    p0 = [0.4, 0.005, 0.05, 2]
    p1 = [1, 0.5, 0.05, 0.05, 3.5]
    # mast_dat_dict = lmem.read_mastfile(mastloc)
    # psi = mast_dat_dict['psi_normal']
    # ne = mast_dat_dict['electron_density(10^20/m^3)']
    # te = mast_dat_dict['electron_temperature(KeV)']
    f_psi = psi -1
    Ne = ne*pow(10, -19)
    Te = te*pow(10, -3)
       
    popt_ne, pcov_ne = curve_fit(tanh_dsa, f_psi, Ne, p0)
    print(popt_ne)
    popt_te, pcov_te = curve_fit(tanh, psi, Te, p1)
    print(popt_te) 
    tanh_ne_fit = tanh_dsa(f_psi, popt_ne[0], popt_ne[1], popt_ne[2], 
                           popt_ne[3])*pow(10, 19)
    tanh_te_fit = tanh(psi, popt_te[0], popt_te[1], popt_te[2], 
                           popt_te[3], popt_te[4])*pow(10, 3)
    
    # plt.figure(1)
    # plt.plot(psi, ne,'o-', color = 'b', label= 'solps electron density')
    # plt.plot(psi, tanh_ne_fit, color='r',lw= 3, label= 'tanh fit')
    # plt.xlabel('Radial coordinate: $R- R_{sep}$')
    # plt.ylabel('Electron Density $n_e\;(m^{-3})$')
    # plt.title('Electron density with fits')
    # plt.legend()
    
    fit_tanh_dic = {'tanh_ne_fit': tanh_ne_fit, 'tanh_te_fit': tanh_te_fit, 
               'popt_ne': popt_ne, 'popt_te': popt_te}
    
    return fit_tanh_dic



def exp_dsa_fit(dsa, neuden):
    pn = [1.015, 200.5]
    
    NeuDen = neuden/ max(neuden)
    
    popt_an, pcov_an = curve_fit(expfit, dsa, NeuDen, pn)
    print(popt_an)
    
    exp_an_fit = expfit(dsa, popt_an[0], popt_an[1])*max(neuden)
    
    fit_exp_dic = {'exp_an_fit': exp_an_fit, 'popt_an': popt_an}
    
    return fit_exp_dic



def exp_psi_fit(psi, neuden):
    pn = [1, 100]
    
    f_psi = psi - 1
    NeuDen = neuden*pow(10, -16)
    # print(NeuDen)
    
    
    popt_an, pcov_an = curve_fit(expfit, f_psi, NeuDen, pn)
    print(popt_an)
    
    exp_an_fit = expfit(f_psi, popt_an[0], popt_an[1])*pow(10, 16)
   
    fit_exp_dic = {'exp_an_fit': exp_an_fit, 'popt_an': popt_an}
    
    # plt.figure(2)
    # plt.plot(psi, neuden,'o-', color = 'green', label= 'solps neutral density')
    # plt.plot(psi, exp_an_fit, color='r',lw= 5, label= 'exponential fit')
    # plt.xlabel('psiN')
    # plt.ylabel('Neutral Atom (D) Density $(m^{-3})$')
    # plt.title('Neutral density with fits')
    # plt.legend()
    
    return fit_exp_dic


def Opacity_calculator(dsa, ne, te, neuden):
    fit_tanh_dic = tanh_dsa_fit(dsa, ne, te)
    tanh_ne_fit = fit_tanh_dic['tanh_ne_fit']
    dn = fit_tanh_dic['popt_ne'][1]
    
    dsa_cut = []
    an_cut = []
    for j in range(len(dsa)):
        if dsa[j] <= dn and dsa[j] >= -dn:
            dsa_cut.append(dsa[j])
            an_cut.append(neuden[j])
    
    dsa_cut = np.asarray(dsa_cut)
    an_cut = np.asarray(an_cut)
    
    fit_exp_dic = exp_dsa_fit(dsa_cut, an_cut)
    exp_an_fit = fit_exp_dic['exp_an_fit']
    efold = 1/fit_exp_dic['popt_an'][1]
    
    
    opq = 2*dn/efold
    
    result_dic = {'tanh_fit': tanh_ne_fit, 'exp_fit': exp_an_fit,
                      'pedestal_width': dn, 'efold_length': efold,
                      'dimensionless_opaqueness': opq}
    
    return result_dic


def Opacity_calculator_psi(psi, ne, te, neuden):
    fit_tanh_dic = tanh_psi_fit(psi, ne, te)
    tanh_ne_fit = fit_tanh_dic['tanh_ne_fit']
    dn = fit_tanh_dic['popt_ne'][1]
    tanh_te_fit = fit_tanh_dic['tanh_te_fit']
    dtn = fit_tanh_dic['popt_te'][2]
    
    ne_ped = (fit_tanh_dic['popt_ne'][0] + fit_tanh_dic['popt_ne'][2])*pow(10, 19)

    "first cut"
    psi_cut = []
    an_cut = []
    for j in range(len(psi)):
        if psi[j] <= dn + 1:
            psi_cut.append(psi[j])
            an_cut.append(neuden[j])
    
    psi_cut = np.asarray(psi_cut)
    an_cut = np.asarray(an_cut)
    
    fit_exp_dic = exp_psi_fit(psi_cut, an_cut)
    
    
    "fit data cut"
    psi_cut_dn = []
    an_fit_cut = []
    for j in range(len(psi)):
        if psi[j] <= dn + 1 and psi[j] >= -dn + 1:
            psi_cut_dn.append(psi[j])
            # an_fit_cut.append(fit_exp_dic['exp_an_fit'][j])

    
    psi_cut_dn = np.asarray(psi_cut_dn)
    an_fit_cut = expfit(psi_cut_dn - 1, fit_exp_dic['popt_an'][0], fit_exp_dic['popt_an'][1])*pow(10, 16)
    
    exp_an_fit = fit_exp_dic['exp_an_fit']
    efold = 1/fit_exp_dic['popt_an'][1]
    
    opq = 2*dn/efold
    
    result_dic = {'tanh_ne_fit': tanh_ne_fit, 'exp_fit': exp_an_fit,
                  'exp_fit_in_width': an_fit_cut, 'electron_density_pedestal': ne_ped,
                  'tanh_te_fit': tanh_te_fit, 'pedestal_width': dn, 
                  'temperature_pedestal_width': dtn, 'psi_cut': psi_cut_dn,
                  'efold_length': efold, 'dimensionless_opaqueness': opq}
    
    return result_dic


def expfit_B(x, A, l, B):  
    return A*np.exp(l*x)+ B


def exp_dsa_fit_B(dsa, neuden):
    
    
    NeuDen = neuden/ max(neuden)
    
    pn = [1.015, 500, min(NeuDen)]
    
    popt_an, pcov_an = curve_fit(expfit_B, dsa, NeuDen, pn)
    print(popt_an)
    
    exp_an_fit = expfit_B(dsa, popt_an[0], popt_an[1], popt_an[2])*max(neuden)
    
    fit_exp_dic = {'exp_an_fit': exp_an_fit, 'popt_an': popt_an}
    
    return fit_exp_dic
