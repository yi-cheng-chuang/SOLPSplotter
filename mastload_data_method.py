# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 21:05:04 2023

@author: user
"""

import glob
import numpy as np
import B2plotter_tool as b2t
import B2plotter_set as b2s

d = b2s.Setting_dic()

dd = b2s.helper_dic()

def mast_b2_dir(a_shift):
    basedrt, topdrt, tpdrt= b2t.set_wdir()
    shift_list = dd['shiftlist']
    
    for s in shift_list:
        if a_shift == s:
            i = shift_list.index(a_shift)
            filename = d['series'][i]
            newbase = '{}/{}/{}/{}/{}'.format(basedrt, d['DEV'], d['Shot'], 
                                              d['shift'][i], filename)
            tbase = '{}/{}/{}/{}'.format(basedrt, d['DEV'], d['Shot'], 
                                         d['shift'][i])
            gbase = '{}/{}/{}'.format(topdrt, d['DEV'], d['Shot'])
            gdir = glob.glob('{}/g{}*'.format(gbase, d['Shot']))
            adrt = '{}/Output'.format(newbase)
            bdrt = '{}/Output2'.format(newbase)
            Attempt = str(b2t.s_number(adrt)[0])
            print(filename)
            print(Attempt)
    
    mastdic = {'basedrt': basedrt, 'topdrt': topdrt, 'tpdrt': tpdrt, 
               'simudir': newbase, 'tbase': tbase, 'expdir': gbase, 'outone': adrt,
               'outtwo': bdrt, 'gdir': gdir}
    
    return mastdic

def read_mastfile(mastfile_loc):
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