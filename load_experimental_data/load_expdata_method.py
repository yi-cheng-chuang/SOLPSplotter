# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 21:24:52 2025

@author: ychuang
"""

import numpy as np


class read_expdata_method:
    
    def __init__(self, DF, data):
        self.DF = DF
        self.data = data
    
    
    def read_mastfile(self, mastfile_loc):
        with open(mastfile_loc, mode='r') as dfile:
            lines = dfile.readlines()
        
        profiles = {}
        nlines_tot = len(lines)
        psi_n = np.zeros(nlines_tot)
        ne = np.zeros(nlines_tot)
        ne_er = np.zeros(nlines_tot)
        te = np.zeros(nlines_tot)
        te_er = np.zeros(nlines_tot)
        
        i = 0
        
        while i < nlines_tot:
            r_line = lines[i].split()
            psi_n[i] = float(r_line[0])
            ne[i] = float(r_line[1])*pow(10, -20)
            ne_er[i] = float(r_line[2])*pow(10, -20)
            te[i] = float(r_line[3])/1000
            te_er[i] = float(r_line[4])/1000
            i += 1

        profiles['psi_normal'] = psi_n
        profiles['electron_density(10^20/m^3)'] = ne
        profiles['density error(10^20/m^3)'] = ne_er
        profiles['electron_temperature(KeV)'] = te
        profiles['temperature error(10^20/m^3)'] = te_er
        return profiles


    def read_fitfile(self, mastfile_loc):
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


