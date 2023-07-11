# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 14:59:41 2019

@author: rmreksoatmodjo

Collection of general Tools to perform oft-repeated SOLPS data analyis and post-processing tasks
"""
import os
import numpy as np
from scipy import interpolate

def set_wdir(): #Function to set correct Working Directory Path depending on which machine is in use
    if os.environ['OS'] == 'Windows_NT':
        if os.environ['USERNAME'] == 'Yi-Cheng':
            basedrt = r"C:/Users/Yi-Cheng/Documents/SOLPS_Data/Simulation_Data"
            topdrt = r"C:/Users/Yi-Cheng/Documents/SOLPS_Data/Experimental_Data"
        elif os.environ['USERNAME'] == 'user':
            basedrt = r"C:/Users/user/Documents/SOLPS data/simulation data"
            topdrt = r"C:/Users/user/Documents/SOLPS data/experiment data"
    else:
        print('please add new directory in tools')
    
    return basedrt, topdrt

def a_number(text):
    name = text.split('\\',1)[1]
    nu = int(name.split('_')[2])

    return [nu, text]

def ex1_number(text):
    name = text.split('\\',1)[1]
    nu = int(name.split('_')[3])

    return [nu, text]

def unit_dic():
    unit = {
        'ne3da.last10': ['Electron density, outboard midplane', 'Electron density: ${n_e}$ (m$^{-3}$)'],
        'te3da.last10': ['Electron temperature, outboard midplane', 'Electron temperature: ${T_e}$ (eV)'],
        'an3da.last10':['Neutral density, outboard midplane', 'Neutral density: ${n_D}$ (m$^{-3}$)'],
        'ne3di.last10': ['Electron density, inboard midplane', 'Electron density: ${n_e}$ (m$^{-3}$)'],
        'te3di.last10': ['Electron temperature, inboard midplane', 'Electron temperature: ${T_e}$ (eV)'],
        'an3di.last10':['Neutral density, inboard midplane', 'Neutral density: ${n_D}$ (m$^{-3}$)'],
        'ne3dl.last10':['Electron density, western midplane', 'Electron density: ${n_e}$ (m$^{-3}$)'],
        'te3dl.last10': ['Electron temperature, western midplane', 'Electron temperature: ${T_e}$ (eV)'],
        'an3dl.last10':['Neutral density, western midplane', 'Neutral density: ${n_D}$ (m$^{-3}$)'],
        'ne3dr.last10':['Electron density, eastern midplane', 'Electron density: ${n_e}$ (m$^{-3}$)'],
        'te3dr.last10': ['Electron temperature, eastern midplane', 'Electron temperature: ${T_e}$ (eV)'],
        'an3dr.last10':['Neutral density, eastern midplane', 'Neutral density: ${n_D}$ (m$^{-3}$)'],
        '1':['Particle density-driven diffusivity','Density-driven diffusivity: D (m$^{2}$/s)'],
        '3':['Ion thermal anomalous diffusivity', 'Ion thermal diffusivity: ${\chi_i}$ (m$^{2}$/s)'],
        '4':['Electron thermal anomalous diffusivity', 'Electron thermal diffusivity: ${\chi_e}$ (m$^{2}$/s)'],
        
        }
    return unit

psi_solps =[0.56591402, 0.59635553, 0.65365526, 0.70507622, 0.75083496,
        0.79132874, 0.82698182, 0.85806206, 0.88490369, 0.90789432,
        0.92738632, 0.94367313, 0.95706941, 0.96795829, 0.97677538,
        0.9838775 , 0.98955578, 0.99415907, 0.99803803, 1.002408  ,
        1.00753157, 1.01263476, 1.01772166, 1.02279374, 1.02785249,
        1.03288158, 1.03794617, 1.04306613, 1.04817989, 1.05328886,
        1.05838546, 1.06347049, 1.06855367, 1.07363646, 1.07872032,
        1.08380671, 1.08889011, 1.09145489]

dsa = [-0.10817856015924884,
-0.09990042377111603,
-0.0846133267035689,
-0.0713664335482361,
-0.05987293296951107,
-0.049898929370704115,
-0.04124897557162459,
-0.03377954856336941,
-0.02737025239404635,
-0.02191068496806607,
-0.017292192441575246,
-0.013432437532237718,
-0.01025503971802598,
-0.00766848237493005,
-0.005577570403317772,
-0.0038936225714468614,
-0.002540280074319473,
-0.0014428354897131068,
-0.0005202548013382713,
0.0005202548013382852,
0.0017412877642430796,
0.0029587093458706415,
0.004173484315468989,
0.005385479979830576,
0.0065946993526905945,
0.007799836184615799,
0.009009678339736343,
0.010226892341373167,
0.011443929068697442,
0.012666850670040336,
0.013891570075914711,
0.015112086991656865,
0.016332974621876967,
0.0175550249181618,
0.018778594838134655,
0.020003969353135337,
0.021230335606628042,
0.021849669506929625]

def psi_to_dsa():
    psi_to_dsa_func = interpolate.interp1d(psi_solps, dsa, fill_value = 'extrapolate')
    
    return psi_to_dsa_func()

def dsa_to_psi(array):
    dsa_to_psi_func = interpolate.interp1d(dsa, psi_solps, fill_value = 'extrapolate')
    
    return dsa_to_psi_func(array)



def read_mastfile(mastfile_loc):
    with open(mastfile_loc, mode='r') as dfile:
        lines = dfile.readlines()
    
    profiles = {}
    nlines_tot = len(lines)
    psi_n = np.zeros(nlines_tot)
    ne = np.zeros(nlines_tot)
    te = np.zeros(nlines_tot)
    i = 0
    
    while i < nlines_tot:
        r_line = lines[i].split()
        psi_n[i] = float(r_line[0])
        ne[i] = float(r_line[1])*pow(10, -20)
        te[i] = float(r_line[3])/1000
        i += 1

    profiles['psi_normal'] = psi_n
    profiles['electron_density(10^20/m^3)'] = ne
    profiles['electron_temperature(KeV)'] = te
    return profiles

def read_fitfile(mastfile_loc):
    with open(mastfile_loc, mode='r') as dfile:
        lines = dfile.readlines()
    
    profiles = {}
    nlines_tot = len(lines)
    psi_n = np.zeros(nlines_tot)
    ne = np.zeros(nlines_tot)
    te = np.zeros(nlines_tot)
    i = 0
    
    while i < nlines_tot:
        r_line = lines[i].split()
        psi_n[i] = float(r_line[0])
        ne[i] = float(r_line[1])*pow(10, 20)
        te[i] = float(r_line[2])*1000
        i += 1

    profiles['psi_normal'] = psi_n
    profiles['electron_density(m^(-3))'] = ne
    profiles['electron_temperature(eV)'] = te
    return profiles

def tanh(r,r0,h,d,b,m):
    return b+(h/2)*(np.tanh((r0-r)/d)+1) + m*(r0-r-d)*np.heaviside(r0-r-d, 1)

def expfit(x,A,l):  #Removed vertical displacement variable B; seemed to cause 'overfitting'
    return A*np.exp(l*x)

def flat_tanh(x,b,h,d):
    return b+(h/2)*(np.tanh(-x/d)+1)

