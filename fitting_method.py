# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 13:28:56 2023

@author: user
"""

import numpy as np
from scipy.optimize import curve_fit
import load_mast_expdata_method as lmem

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

def exp_dsa_fit(dsa, neuden):
    pn = [1.015, 200.5]
    
    NeuDen = neuden/ max(neuden)
    
    popt_an, pcov_an = curve_fit(expfit, dsa, NeuDen, pn)
    print(popt_an)
    
    exp_an_fit = expfit(dsa, popt_an[0], popt_an[1])*max(neuden)
    
    fit_exp_dic = {'exp_an_fit': exp_an_fit, 'popt_an': popt_an}
    
    return fit_exp_dic
