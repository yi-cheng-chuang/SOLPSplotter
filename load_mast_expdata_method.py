# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 21:05:04 2023

@author: user
"""

import glob
import numpy as np 
import SOLPS_set as sps


d = sps.mast_comp_dic()
od = sps.Setting_dic()
# print(type(d))

def mast_base_dir():
    basedrt, topdrt = sps.set_wdir()
    gbase = '{}/{}/{}'.format(topdrt, od['DEV'], d['Shot'])
    gdir = glob.glob('{}/g{}*'.format(gbase, d['Shot']))
    
    shift_list = list(d['shift_dic'].keys())
    # print(type(shift_list))
    a_shift = d['a_shift']
    for aa in shift_list:
        if aa == a_shift:  
            filename = d['series_dic'][aa]
            shift = d['shift_file_dic'][aa]
            newbase = '{}/{}/{}/{}/{}'.format(basedrt, od['DEV'], 
                                            d['Shot'], shift, filename)
            tbase = '{}/{}/{}/{}'.format(basedrt, od['DEV'], d['Shot'], shift)
            adir = {}
            for i in d['Output']:
                adir[i] = '{}/{}'.format(newbase, i)
     
    attempt = str(sps.s_number(adir['Output'], series_flag= None)[0])
    shift_value = d['shift_dic'][a_shift]
    
    mast_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 'gbase': gbase, 
                    'gdir': gdir, 'simudir': newbase, 'simutop': tbase, 
                    'outputdir': adir}

    return mast_basedir, attempt, shift_value


def terminal_single_dir():
    basedrt, topdrt = sps.set_wdir()
    gdir = glob.glob('{}/g{}*'.format(topdrt, d['Shot']))
    
    aa = d['a_shift']
 
    filename = d['series_dic'][aa]
    shift = d['shift_file_dic'][aa]
    newbase = '{}/{}/{}'.format(basedrt, shift, filename)
    tbase = '{}/{}'.format(basedrt, shift)
    adir = {}
    for i in d['Output']:
        adir[i] = '{}/{}'.format(newbase, i)
     
    attempt = str(sps.s_number(adir['Output'], series_flag= None)[0])
    shift_value = d['shift_dic'][aa]
    
    mast_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 
                    'gdir': gdir, 'simudir': newbase, 'simutop': tbase, 
                    'outputdir': adir}

    return mast_basedir, attempt, shift_value





def mast_withshift_dir():
    
    mwd = sps.mast_comp_dic_withshift()
    # print(mwd)
     
    basedrt, topdrt = sps.set_wdir()
    gbase = '{}/{}/{}'.format(topdrt, od['DEV'], d['Shot'])
    gdir = glob.glob('{}/g{}*'.format(gbase, d['Shot']))
    shift_list = list(mwd['shift_dic'].keys())
    # print(type(shift_list))
    a_shift = mwd['multi_shift']
    simudir_dic = {}
    simutop_dic = {}
    outputdir_dic = {}
    att_dic = {}
    for aa in a_shift:
        for s in shift_list:
            if aa == s:
                i = shift_list.index(aa)
                filename = mwd['series'][i]
                
                # print('look for series :')
                # print(mwd['series'])
                
                newbase = '{}/{}/{}/{}/{}'.format(basedrt,od['DEV'], d['Shot'],
                                           mwd['shift_filelist'][i], filename)
                tbase = '{}/{}/{}/{}'.format(basedrt, od['DEV'], 
                                      d['Shot'], mwd['shift_filelist'][i])
                adir = {}
                for i in d['Output']:
                    adir[i] = '{}/{}'.format(newbase, i)
                
                print(adir['Output'])
                
                att_dic[aa] = str(sps.s_number(adir['Output'], series_flag= None)[0])
                
                simudir_dic[aa] = newbase
                simutop_dic[aa] = tbase
                outputdir_dic[aa] = adir
    
    mast_withshift_dir_dic = {'basedrt': basedrt, 'topdrt': topdrt, 
                              'gbase': gbase, 'gdir': gdir, 
                              'simudir': simudir_dic, 'simutop': simutop_dic, 
                              'outputdir': outputdir_dic}
    
    
    return mast_withshift_dir_dic, att_dic


mwd = sps.mastcomp_withshift_compare()

def mast_withshift_cp_dir():
    basedrt, topdrt = sps.set_wdir()
    gbase = '{}/{}/{}'.format(topdrt, od['DEV'], d['Shot'])
    gdir = glob.glob('{}/g{}*'.format(gbase, d['Shot']))
    shift_list = list(mwd['shift_dic'].keys())
    # print(type(shift_list))
    a_shift = mwd['multi_shift']
    simudir_dic = {}
    simutop_dic = {}
    outputdir_dic = {}
    att_dic = {}
    for aa in a_shift:
        for s in shift_list:
            if aa == s:
                i = shift_list.index(aa)
                filename = mwd['series'][i]
                newbase_cp = {'fixed': {}, 'flux': {}}
                
                
                newbase_cp['fixed'] = '{}/{}/{}/{}/{}'.format(basedrt,od['DEV'], d['Shot'],
                                                mwd['shift_filelist'][i], filename[0])
                
                newbase_cp['flux'] = '{}/{}/{}/{}/{}'.format(basedrt,od['DEV'], d['Shot'],
                                                mwd['shift_filelist'][i], filename[1])
                
                
                tbase = '{}/{}/{}/{}'.format(basedrt, od['DEV'], 
                                          d['Shot'], mwd['shift_filelist'][i])
                    
                
                adir = {'fixed': {}, 'flux': {}}
                for k in ['fixed', 'flux']:
                    for i in d['Output']:
                        adir[k][i] = '{}/{}'.format(newbase_cp[k], i)
                    
                    
                 
                att_dic[aa] = str(sps.s_number(adir['fixed']['Output'], series_flag= None)[0])
                
                simudir_dic[aa] = newbase_cp
                simutop_dic[aa] = tbase
                outputdir_dic[aa] = adir
    
    mast_withshift_dir_dic = {'basedrt': basedrt, 'topdrt': topdrt, 
                              'gbase': gbase, 'gdir': gdir, 
                              'simudir': simudir_dic, 'simutop': simutop_dic, 
                              'outputdir': outputdir_dic}
    
    
    return mast_withshift_dir_dic, att_dic











mcds = sps.mast_comp_dir_series()

def mast_series_dir(series_flag):
    if series_flag == 'change_den':
        mcds = sps.mast_comp_dir_series()
    elif series_flag == 'eireneN':
        mcds = sps.mast_comp_dir_eireneN()
    elif series_flag == 'change_temp':
        mcds = sps.mast_comp_dir_tempscan()
    else:
        print('please check series_flag')
       
    basedrt, topdrt = sps.set_wdir()
    gbase = '{}/{}/{}'.format(topdrt, od['DEV'], mcds['Shot'])
    gdir = glob.glob('{}/g{}*'.format(gbase, mcds['Shot']))
    newbase = glob.glob('{}/{}/{}/{}/*{}'.format(basedrt,od['DEV'], mcds['Shot'], 
                                       mcds['shift'], mcds['tail']))
    tbase = '{}/{}/{}/{}'.format(basedrt, od['DEV'], 
                                 mcds['Shot'], mcds['shift'])

    attempt_dic = {}
    new_dic = {}
    for i in newbase:
        if series_flag == 'change_den':
            attempt_dic[sps.s_number(i, series_flag)[0][0]] = sps.s_number(i, series_flag)[0][1]
            new_dic[sps.s_number(i, series_flag)[0][0]] = i
        elif series_flag == 'eireneN':
            attempt_dic[sps.s_number(i, series_flag)[0][1]] = sps.s_number(i, series_flag)[0][0]
            new_dic[sps.s_number(i, series_flag)[0][1]] = i
        elif series_flag == 'change_temp':
            attempt_dic[sps.s_number(i, series_flag)[0][0]] = sps.s_number(i, series_flag)[0][1]
            new_dic[sps.s_number(i, series_flag)[0][0]] = i
    # print(attempt_list)
    
    adir = {}
    for ii in attempt_dic.keys():
        adir[ii] = {}
        for j in mcds['Output']:
            adir[ii][j] = '{}/{}'.format(new_dic[ii], j)
    
    mast_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 'gbase': gbase, 
                    'gdir': gdir, 'simudir': new_dic, 'simutop': tbase, 
                    'outputdir': adir}

    return mast_basedir, attempt_dic


mcds = sps.mast_comp_dir_series()

def series_twocompare_dir(series_flag):
    if series_flag == 'two_compare':
        mcds = sps.mast_twocompare_dir()
    else:
        print('please check series_flag')
       
    basedrt, topdrt = sps.set_wdir()
    gbase = '{}/{}/{}'.format(topdrt, od['DEV'], mcds['Shot'])
    gdir = glob.glob('{}/g{}*'.format(gbase, mcds['Shot']))
    newbase = []
    for aa in mcds['fname_list']:
        newbase.append('{}/{}/{}/{}/{}'.format(basedrt,od['DEV'], mcds['Shot'], 
                                           mcds['series_name'], aa))
        
    
    tbase = '{}/{}/{}/{}'.format(basedrt, od['DEV'], mcds['Shot'], mcds['series_name'])

    attempt_dic = {}
    new_dic = {}
    
    if series_flag == 'two_compare':
        for ab in newbase:
            
            if 'std_a' in ab:
                
                shotnum = sps.s_number(ab, series_flag)[0][1]
                print('std shotnum is {}'.format(str(shotnum)))
                attempt_dic['std'] = shotnum
                new_dic['std'] = ab
            
            elif 'save_a' in ab:
                
                shotnum = sps.s_number(ab, series_flag)[0][1]
                print('save shotnum is {}'.format(str(shotnum)))
                attempt_dic['save'] = shotnum
                new_dic['save'] = ab
                
            
    # print(attempt_list)
    
    adir = {}
    for ii in attempt_dic.keys():
        adir[ii] = {}
        for j in mcds['Output']:
            adir[ii][j] = '{}/{}'.format(new_dic[ii], j)
    
    mast_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 'gbase': gbase, 
                    'gdir': gdir, 'simudir': new_dic, 'simutop': tbase, 
                    'outputdir': adir}

    return mast_basedir, attempt_dic




def two_layer_dic(key_a, key_b):
    
    # Initialize an empty dictionary
    twinscan_dic = {}
    
    # Populate the dictionary using nested loops
    for ka in key_a:
        
        twinscan_dic[ka] = {}
        
        for kb in key_b:
            
            twinscan_dic[ka][kb] = {}
    
    return twinscan_dic


def series_terminal_dir(series_flag, dir_comp_dic):

    mcds = dir_comp_dic
    
    # if series_flag == 'terminal_test':
    #     mcds = sps.terminal_series_comp_dir(tail = '_leakbsol_nts5_a', 
    #                         filename = 'org_densityscan_027205')
    # else:
    #     print('please check series_flag')
       
    basedrt, topdrt = sps.set_wdir()
    gdir = glob.glob('{}/g{}*'.format(topdrt, mcds['Shot']))
    newbase = glob.glob('{}/{}/*{}'.format(basedrt, mcds['filename'], 
                                                 mcds['tail']))
    tbase = '{}/{}'.format(basedrt, mcds['filename'])

    attempt_dic = {}
    new_dic = {}
    twinscan_dic = {}
    
    ds_key = []
    ts_key = []
    
    
    print('this is mcds')
    print(type(mcds['denscan_list'][3]))
    
    for x in mcds['denscan_list']:
        ds_key.append('{:.3f}'.format(x))
        
    print(ds_key)
    for x in mcds['tempscan_list']:
        ts_key.append('{:.3f}'.format(x))
        

    print(ts_key)
    
    
    # ds_key = [str(x) for x in mcds['denscan_list']]
    # ts_key = [str(x) for x in mcds['tempscan_list']]
    
    
    twinscan_dic = two_layer_dic(key_a = ds_key, key_b = ts_key)
    
    
    for i in newbase:
        if series_flag == 'terminal_test':
            attempt_dic[sps.s_number(i, series_flag)[0][0]] = sps.s_number(i, series_flag)[0][1]
            new_dic[sps.s_number(i, series_flag)[0][0]] = i
        
        elif series_flag == 'twin_scan':
            attempt_dic[sps.atp_number(i, series_flag)[0]] = sps.atp_number(i, series_flag)[1]
            st = sps.atp_number(i, series_flag)[0]
            twinscan_dic[str(st[0])][str(st[1])] = i
        
        else:
            print('please check series_flag')
    # print(attempt_list)
    
    adir = {}
    
    if series_flag != 'twin_scan':
        for ii in attempt_dic.keys():
            adir[ii] = {}
            for j in mcds['Output']:
                adir[ii][j] = '{}/{}'.format(new_dic[ii], j)
    elif series_flag == 'twin_scan':
        pass
    else:
        print('series_terminal_dir, please check series_flag!')
    
    if series_flag == 'twin_scan':
        mast_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 
                        'gdir': gdir, 'simudir': twinscan_dic, 'simutop': tbase}
    else:
        mast_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 
                        'gdir': gdir, 'simudir': new_dic, 'simutop': tbase, 
                        'outputdir': adir}
        

    return mast_basedir, attempt_dic


def twinscan_local_dir(series_flag, dir_comp_dic):

    mcds = dir_comp_dic
       
    basedrt, topdrt = sps.set_wdir()
    gbase = '{}/{}/{}'.format(topdrt, od['DEV'], mcds['Shot']) 
    gdir = glob.glob('{}/g{}*'.format(gbase, mcds['Shot']))
    simbase = '{}/{}/{}'.format(basedrt, od['DEV'], d['Shot'])
    newbase = glob.glob('{}/{}/*{}'.format(simbase, mcds['filename'], 
                                                 mcds['tail']))
    tbase = '{}/{}'.format(simbase, mcds['filename'])

    attempt_dic = {}
    new_dic = {}
    twinscan_dic = {}
    ds_key = []
    ts_key = []
    
    print('this is mcds')
    print(type(mcds['denscan_list'][3]))
    
    for x in mcds['denscan_list']:
        ds_key.append('{:.3f}'.format(x))
        
    print(ds_key)
    for x in mcds['tempscan_list']:
        ts_key.append('{:.3f}'.format(x))
        

    print(ts_key)
    
    
    twinscan_dic = two_layer_dic(key_a = ds_key, key_b = ts_key)
    print(twinscan_dic)
    
    for i in newbase:
        attempt_dic[sps.atp_number(i, series_flag)[0]] = sps.atp_number(i, series_flag)[1]
        # print(sps.atp_number(i, series_flag))
        st = sps.atp_number(i, series_flag)[0]
        # print(st)
        # print(i)
        twinscan_dic[str(st[0])][str(st[1])] = i
    
    
    mast_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 'gbase': gbase, 
                    'gdir': gdir, 'simudir': twinscan_dic, 'simutop': tbase}
        

    return mast_basedir, attempt_dic







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