# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 21:05:04 2023

@author: user
"""

import glob
import numpy as np 
import B2plotter_set as b2s


d = b2s.mast_comp_dic()
od = b2s.Setting_dic()
# print(type(d))

def mast_base_dir():
    basedrt, topdrt, tpdrt= b2s.set_wdir()
    gbase = '{}/{}/{}'.format(topdrt, od['DEV'], d['Shot'])
    gdir = glob.glob('{}/g{}*'.format(gbase, d['Shot']))
    a_shift = d['a_shift']
    shift_list = list(d['shiftdic'].keys())
    for s in shift_list:
        if a_shift == s:
            i = shift_list.index(a_shift)
            filename = d['series'][i]
            newbase = '{}/{}/{}/{}/{}'.format(basedrt,od['DEV'], d['Shot'],
                                              d['shift'][i], filename)
            tbase = '{}/{}/{}/{}'.format(basedrt, od['DEV'], 
                                         d['Shot'], d['shift'][i])
            adir = {}
            for i in d['Output']:
                adir[i] = '{}/{}'.format(newbase, i)
             
            attempt = str(b2s.s_number(adir['Output'])[0])
            
            mastdic = {'simudir': newbase, 'expdir': tbase, 
                           'outputdir': adir}
    mast_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 'gbase': gbase, 
                    'gdir': gdir, 'adir': mastdic}

    
    return mast_basedir, attempt

def mast_shift_dir():
    basedrt, topdrt, tpdrt= b2s.set_wdir()
    shift_list = list(d['shiftdic'].keys())
    # print(type(shift_list))
    a_shift = d['multi_shift']
    mastdic = {}
    att_dic = {}
    for aa in a_shift:
        for s in shift_list:
            if aa == s:
                i = shift_list.index(aa)
                filename = d['series'][i]
                newbase = '{}/{}/{}/{}/{}'.format(basedrt,od['DEV'], d['Shot'],
                                                  d['shift'][i], filename)
                tbase = '{}/{}/{}/{}'.format(basedrt, od['DEV'], 
                                             d['Shot'], d['shift'][i])
                adir = {}
                for i in d['Output']:
                    adir[i] = '{}/{}'.format(newbase, i)
                 
                att_dic[aa] = str(b2s.s_number(adir['Output'])[0])
                
                mastdic[aa] = {'simudir': newbase, 'expdir': tbase, 
                               'outputdir': adir}
    shift_dir = mastdic
    
    return shift_dir, att_dic




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