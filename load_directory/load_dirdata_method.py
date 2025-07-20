# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 21:05:04 2023

@author: user
"""

import glob
from SOLPS_input.directory_input import directory_input
from load_directory.grab_attempt_number import grab_aptn_method
from pathlib import Path




# print(type(d))


class load_dir_method:
    
    def __init__(self, DF, data, di: directory_input, gam: grab_aptn_method):
        # self.init = 'init_load_dir_method'
        
        self.DF = DF
        self.data = data
        self.di = di
        self.gam = gam


    def mastu_base_dir(self):
        nd = self.di.mastu_comp_dic()
        basedrt, topdrt = self.di.set_wdir()
        
        if self.DF.DEV == 'cross_machine':
            
            gbase = '{}/{}/{}'.format(topdrt, 'mastu', nd['Shot'])
            g_path = glob.glob('{}/*'.format(gbase))[0]
            gdir = Path(g_path).as_posix()
        
        else:
            
            gbase = '{}/{}/{}'.format(topdrt, self.DF.DEV, nd['Shot'])
            g_path = glob.glob('{}/*'.format(gbase))[0]
            gdir = Path(g_path).as_posix()
            
        
        
        # all_list = glob.glob('{}/*'.format(gbase))
        # print(all_list)
        

        series = nd['series']
        filename = nd['filename']
        
        if self.DF.DEV == 'cross_machine':
            
            newbase = '{}/{}/{}/{}/{}'.format(basedrt, 'mastu', 
                                            nd['Shot'], series, filename)
            tbase = '{}/{}/{}/{}'.format(basedrt, 'mastu', nd['Shot'], series)
        
        else:
            newbase = '{}/{}/{}/{}/{}'.format(basedrt, self.DF.DEV, 
                                            nd['Shot'], series, filename)
            tbase = '{}/{}/{}/{}'.format(basedrt, self.DF.DEV, nd['Shot'], series)
            
        

         
        attempt = str(self.gam.mastu_atp_number(newbase, usage = 'load_dir')[0])
        
        print('mastu attempt number is {}'.format(attempt))
        
        
        
        mastu_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 'gbase': gbase, 
                        'gdir': gdir, 'simudir': newbase, 'simutop': tbase}

        return mastu_basedir, attempt
    
    
    
    def mast_base_dir(self):
        d = self.di.mast_comp_dic()
        basedrt, topdrt = self.di.set_wdir()
        
        
        if self.DF.Dnames == 'mastu_mast':
            
            gbase = '{}/{}/{}'.format(topdrt, 'mast', d['Shot'])
            g_path = glob.glob('{}/g{}*'.format(gbase, d['Shot']))[0]
            print(g_path)
            gdir = Path(g_path).as_posix()
        
        else:
            
            gbase = '{}/{}/{}'.format(topdrt, self.DF.DEV, d['Shot'])
            g_path = glob.glob('{}/g{}*'.format(gbase, d['Shot']))[0]
            print(g_path)
            gdir = Path(g_path).as_posix()
            
        
        print(gdir)
        
        shift_list = list(d['shift_dic'].keys())
        # print(type(shift_list))
        a_shift = d['a_shift']
        for aa in shift_list:
            if aa == a_shift:  
                filename = d['series_dic'][aa]
                shift = d['shift_file_dic'][aa]
                
                if self.DF.DEV == 'cross_machine':
                    
                    newbase = '{}/{}/{}/{}/{}'.format(basedrt, 'mast', 
                                                    d['Shot'], shift, filename)
                    tbase = '{}/{}/{}/{}'.format(basedrt, 'mast', d['Shot'], shift)
                
                else:
                    
                    newbase = '{}/{}/{}/{}/{}'.format(basedrt, self.DF.DEV, 
                                                    d['Shot'], shift, filename)
                    tbase = '{}/{}/{}/{}'.format(basedrt, self.DF.DEV, d['Shot'], shift)
                
                
                
                
                adir = {}
                for i in d['Output']:
                    adir[i] = '{}/{}'.format(newbase, i)
         
        attempt = str(self.gam.s_number(adir['Output'])[0])
        shift_value = d['shift_dic'][a_shift]
        
        mast_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 'gbase': gbase, 
                        'gdir': gdir, 'simudir': newbase, 'simutop': tbase, 
                        'outputdir': adir}

        return mast_basedir, attempt, shift_value





    def terminal_single_dir(self):
        
        d = self.di.mast_comp_dic()
        basedrt, topdrt = self.di.set_wdir()
        gdir = glob.glob('{}/g{}*'.format(topdrt, d['Shot']))
        
        aa = d['a_shift']
     
        filename = d['series_dic'][aa]
        shift = d['shift_file_dic'][aa]
        newbase = '{}/{}/{}'.format(basedrt, shift, filename)
        tbase = '{}/{}'.format(basedrt, shift)
        adir = {}
        for i in d['Output']:
            adir[i] = '{}/{}'.format(newbase, i)
         
        attempt = str(self.gam.s_number(adir['Output'])[0])
        shift_value = d['shift_dic'][aa]
        
        mast_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 
                        'gdir': gdir, 'simudir': newbase, 'simutop': tbase, 
                        'outputdir': adir}

        return mast_basedir, attempt, shift_value





    def mast_withshift_dir(self):
        
        mwd = self.di.mast_comp_dic_withshift()
        d = self.di.mast_comp_dic()
        # print(mwd)
         
        basedrt, topdrt = self.di.set_wdir()
        gbase = '{}/{}/{}'.format(topdrt, self.DF.DEV, d['Shot'])
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
                    
                    newbase = '{}/{}/{}/{}/{}'.format(basedrt, self.DF.DEV, d['Shot'],
                                               mwd['shift_filelist'][i], filename)
                    tbase = '{}/{}/{}/{}'.format(basedrt, self.DF.DEV, 
                                          d['Shot'], mwd['shift_filelist'][i])
                    adir = {}
                    for i in d['Output']:
                        adir[i] = '{}/{}'.format(newbase, i)
                    
                    print(adir['Output'])
                    
                    att_dic[aa] = str(self.gam.s_number(adir['Output'])[0])
                    
                    simudir_dic[aa] = newbase
                    simutop_dic[aa] = tbase
                    outputdir_dic[aa] = adir
        
        mast_withshift_dir_dic = {'basedrt': basedrt, 'topdrt': topdrt, 
                                  'gbase': gbase, 'gdir': gdir, 
                                  'simudir': simudir_dic, 'simutop': simutop_dic, 
                                  'outputdir': outputdir_dic}
        
        
        return mast_withshift_dir_dic, att_dic


    

    def mast_withshift_cp_dir(self):
        
        mwd = self.di.mastcomp_withshift_compare()
        d = self.di.mast_comp_dic()
        basedrt, topdrt = self.di.set_wdir()
        
        gbase = '{}/{}/{}'.format(topdrt, self.DF.DEV, d['Shot'])
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
                    
                    
                    newbase_cp['fixed'] = '{}/{}/{}/{}/{}'.format(basedrt, self.DF.DEV, d['Shot'],
                                                    mwd['shift_filelist'][i], filename[0])
                    
                    newbase_cp['flux'] = '{}/{}/{}/{}/{}'.format(basedrt, self.DF.DEV, d['Shot'],
                                                    mwd['shift_filelist'][i], filename[1])
                    
                    
                    tbase = '{}/{}/{}/{}'.format(basedrt, self.DF.DEV, 
                                              d['Shot'], mwd['shift_filelist'][i])
                        
                    
                    adir = {'fixed': {}, 'flux': {}}
                    for k in ['fixed', 'flux']:
                        for i in d['Output']:
                            adir[k][i] = '{}/{}'.format(newbase_cp[k], i)
                        
                        
                     
                    att_dic[aa] = str(self.gam.s_number(adir['fixed']['Output'])[0])
                    
                    simudir_dic[aa] = newbase_cp
                    simutop_dic[aa] = tbase
                    outputdir_dic[aa] = adir
        
        mast_withshift_dir_dic = {'basedrt': basedrt, 'topdrt': topdrt, 
                                  'gbase': gbase, 'gdir': gdir, 
                                  'simudir': simudir_dic, 'simutop': simutop_dic, 
                                  'outputdir': outputdir_dic}
        
        
        return mast_withshift_dir_dic, att_dic



    

    def mast_series_dir(self):
        
        mcds = self.di.mast_comp_dir_series()
        
        
        if self.DF.series_flag == 'change_den':
            mcds = self.di.mast_comp_dir_series()
        elif self.DF.series_flag == 'eireneN':
            mcds = self.di.mast_comp_dir_eireneN()
        elif self.DF.series_flag == 'change_temp':
            mcds = self.di.mast_comp_dir_tempscan()
        else:
            print('please check series_flag')
           
        basedrt, topdrt = self.di.set_wdir()
        gbase = '{}/{}/{}'.format(topdrt, self.DF.DEV, mcds['Shot'])
        gdir = glob.glob('{}/g{}*'.format(gbase, mcds['Shot']))
        newbase = glob.glob('{}/{}/{}/{}/*{}'.format(basedrt, self.DF.DEV, mcds['Shot'], 
                                           mcds['shift'], mcds['tail']))
        tbase = '{}/{}/{}/{}'.format(basedrt, self.DF.DEV, 
                                     mcds['Shot'], mcds['shift'])

        attempt_dic = {}
        new_dic = {}
        for i in newbase:
            if self.DF.series_flag == 'change_den':
                attempt_dic[self.gam.s_number(i, self.DF.series_flag)[0][0]] = self.gam.s_number(i, self.DF.series_flag)[0][1]
                new_dic[self.gam.s_number(i, self.DF.series_flag)[0][0]] = i
            elif self.DF.series_flag == 'eireneN':
                attempt_dic[self.gam.s_number(i, self.DF.series_flag)[0][1]] = self.gam.s_number(i, self.DF.series_flag)[0][0]
                new_dic[self.gam.s_number(i, self.DF.series_flag)[0][1]] = i
            elif self.DF.series_flag == 'change_temp':
                attempt_dic[self.gam.s_number(i, self.DF.series_flag)[0][0]] = self.gam.s_number(i, self.DF.series_flag)[0][1]
                new_dic[self.gam.s_number(i, self.DF.series_flag)[0][0]] = i
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


    

    def series_twocompare_dir(self):
        
        mcds = self.di.mast_comp_dir_series()
        if self.DF.series_flag == 'two_compare':
            mcds = self.di.mast_twocompare_dir()
        else:
            print('please check series_flag')
           
        basedrt, topdrt = self.di.set_wdir()
        gbase = '{}/{}/{}'.format(topdrt, self.DF.DEV, mcds['Shot'])
        gdir = glob.glob('{}/g{}*'.format(gbase, mcds['Shot']))
        newbase = []
        for aa in mcds['fname_list']:
            newbase.append('{}/{}/{}/{}/{}'.format(basedrt, self.DF.DEV, mcds['Shot'], 
                                               mcds['series_name'], aa))
            
        
        tbase = '{}/{}/{}/{}'.format(basedrt, self.DF.DEV, mcds['Shot'], mcds['series_name'])

        attempt_dic = {}
        new_dic = {}
        
        if self.DF.series_flag == 'two_compare':
            for ab in newbase:
                
                if 'std_a' in ab:
                    
                    shotnum = self.gam.s_number(ab, self.DF.series_flag)[0][1]
                    print('std shotnum is {}'.format(str(shotnum)))
                    attempt_dic['std'] = shotnum
                    new_dic['std'] = ab
                
                elif 'save_a' in ab:
                    
                    shotnum = self.gam.s_number(ab, self.DF.series_flag)[0][1]
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




    def two_layer_dic(self, key_a, key_b):
        
        # Initialize an empty dictionary
        twinscan_dic = {}
        
        # Populate the dictionary using nested loops
        for ka in key_a:
            
            twinscan_dic[ka] = {}
            
            for kb in key_b:
                
                twinscan_dic[ka][kb] = {}
        
        return twinscan_dic





    def series_terminal_dir(self, dir_comp_dic):

        mcds = dir_comp_dic
        
        # if series_flag == 'terminal_test':
        #     mcds = sps.terminal_series_comp_dir(tail = '_leakbsol_nts5_a', 
        #                         filename = 'org_densityscan_027205')
        # else:
        #     print('please check series_flag')
           
        basedrt, topdrt = self.di.set_wdir()
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
        
        
        twinscan_dic = self.two_layer_dic(key_a = ds_key, key_b = ts_key)
        
        
        for i in newbase:
            if self.DF.series_flag == 'terminal_test':
                attempt_dic[self.gam.s_number(i, self.DF.series_flag)[0][0]] = self.gam.s_number(i, self.DF.series_flag)[0][1]
                new_dic[self.gam.s_number(i, self.DF.series_flag)[0][0]] = i
            
            elif self.DF.series_flag == 'twin_scan':
                attempt_dic[self.gam.atp_number(i, self.DF.series_flag)[0]] = self.gam.atp_number(i, self.DF.series_flag)[1]
                st = self.gam.atp_number(i, self.DF.series_flag)[0]
                twinscan_dic[str(st[0])][str(st[1])] = i
            
            else:
                print('please check series_flag')
        # print(attempt_list)
        
        adir = {}
        
        if self.DF.series_flag != 'twin_scan':
            for ii in attempt_dic.keys():
                adir[ii] = {}
                for j in mcds['Output']:
                    adir[ii][j] = '{}/{}'.format(new_dic[ii], j)
        elif self.DF.series_flag == 'twin_scan':
            pass
        else:
            print('series_terminal_dir, please check series_flag!')
        
        if self.DF.series_flag == 'twin_scan':
            mast_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 
                            'gdir': gdir, 'simudir': twinscan_dic, 'simutop': tbase}
        else:
            mast_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 
                            'gdir': gdir, 'simudir': new_dic, 'simutop': tbase, 
                            'outputdir': adir}
            

        return mast_basedir, attempt_dic




    def twinscan_local_dir(self, dir_comp_dic):
        d = self.di.mast_comp_dic()
        mcds = dir_comp_dic
           
        basedrt, topdrt = self.di.set_wdir()
        gbase = '{}/{}/{}'.format(topdrt, self.DF.DEV, mcds['Shot']) 
        gdir = glob.glob('{}/g{}*'.format(gbase, mcds['Shot']))
        simbase = '{}/{}/{}'.format(basedrt, self.DF.DEV, d['Shot'])
        newbase = glob.glob('{}/{}/*{}'.format(simbase, mcds['filename'], 
                                                     mcds['tail']))
        tbase = '{}/{}'.format(simbase, mcds['filename'])

        attempt_dic = {}
        new_dic = {}
        twinscan_dic = {}
        ds_key = []
        ts_key = []
        
        # print('this is mcds')
        # print(type(mcds['denscan_list'][3]))
        
        for x in mcds['denscan_list']:
            ds_key.append('{:.3f}'.format(x))
            
        # print(ds_key)
        for x in mcds['tempscan_list']:
            ts_key.append('{:.3f}'.format(x))
            

        # print(ts_key)
        
        
        twinscan_dic = self.two_layer_dic(key_a = ds_key, key_b = ts_key)
        
        # print(twinscan_dic)
        
        
        
        
        for i in newbase:
            attempt_dic[self.gam.atp_number(i)[0]] = self.gam.atp_number(i)[1]
            # print(sps.atp_number(i, series_flag))
            st = self.gam.atp_number(i)[0]
            # print(st)
            # print(i)
            twinscan_dic[str(st[0])][str(st[1])] = i
        
        
        mast_basedir = {'basedrt': basedrt, 'topdrt': topdrt, 'gbase': gbase, 
                        'gdir': gdir, 'simudir': twinscan_dic, 'simutop': tbase}
            

        return mast_basedir, attempt_dic









