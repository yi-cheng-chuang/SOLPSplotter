# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 12:38:48 2023

@author: user
"""

import os
import re
import numpy as np



def Setting_dic():
    
    if os.environ['OS'] == 'Windows_NT':
        terminal = False
    elif os.environ['OS'] == '5.14.0-362.24.1.el9_3.0.1.x86_64':
        terminal = True
    else:
        print('there is a bug at Setting_dic function or unrecognized operating system')
    
    
    
    set_dic = {'DEV': 'mast', 'minor_rad': 0.5, 'withshift': False, 'withseries': True,
               'Parameters': P, 'series_flag': 'two_compare',
    'series_filename': 'org_25scan_027205', 'series_tail': '_leakbsol_nts5_a',
               'Publish': 'b2plottersetting', 'terminal': terminal}
    return set_dic

series_flag = ['eireneN','change_den','change_temp']


def mast_comp_dic():
    a_shift = 'org'
    shift_dic = {'org': 0, 'dot3': 0.3, 'dot5': 0.5, 'dot7': 0.7, 'one': 1}
    
    twinscan = True
    
    if twinscan:
        series_name = 'org_cfluxb_std'
        file_name = '77_nf5.52tf4.11_save_a'
    
    else:
        series_name = 'org_cfluxb_std'
        file_name = '76_nf5.52tf4.11_std_a'
        
    
    
    shift_file_dic = {'org': series_name,'dot3': 'dot3','dot5': 'dot5',
                          'dot7': 'dot7','one': 'one_LS'}
    
    
    series_dic = {'org': file_name, 'dot3': '16_n900000_leakbtarnsol_dot3_a', 
                  'dot5': '26_n100000_leakagebou_dot5_a', 'dot7': '14_n100000_leakagebou_dot7_a', 
                  'one': '33_n100000_leakagebou_one_a'}
    
    outputlist = ['Output', 'Output2', 'EirOutput']
    mast_dir_dic = {'Shot': '027205', 'shift_dic': shift_dic, 
                    'shift_file_dic': shift_file_dic, 'series_dic': series_dic, 
                    'a_shift': a_shift, 'Output': outputlist}
    
    return mast_dir_dic


def mast_twocompare_dir():
    a_shift = 'org'
    shift = 0
    series_name = 'org_cfluxb_std'
    fname_list = ['76_nf5.52tf4.11_std_a', '79_nf5.52tf4.11_save_a']
    outputlist = ['Output', 'Output2', 'EirOutput']
    mast_series_dir_dic = {'Shot': '027205', 'series_name': series_name, 'shift_value': shift,
                    'fname_list': fname_list, 'a_shift': a_shift, 'Output': outputlist}
    
    return mast_series_dir_dic




def mast_comp_dic_withshift():
    multi_shift = ['org', 'dot3', 'dot5', 'dot7']
    # multi_shift = ['org', 'dot3', 'dot5', 'dot7', 'one']
    shift_dic = {'org': 0, 'dot3': 0.3, 'dot5': 0.5, 'dot7': 0.7, 'one': 1}
    shift = ['org_new_series', 'dot3', 'dot5', 'dot7', 'one_LS']
    tail = {'org': 'nts_a', 'dot3': 'dot3_a', 'dot5': 'dot5_a', 'dot7': 'dot7_a',
            'one': 'one_a'}
    
    
    series = ['75_n900000_leakbtarnsol_nts5_a', '16_n900000_leakbtarnsol_dot3_a', '27_n900000_leakbtarnsol_dot5_a', 
              '15_n900000_leakbtarnsol_dot7_a', '32_n100000_m12_one_a']
    
    # series = ['73_n100000_n5e3et1e2_nts5_a', '15_n100000_leakagebou_dot3_a', '26_n100000_leakagebou_dot5_a', 
    #           '14_n100000_leakagebou_dot7_a', '33_n100000_leakagebou_one_a']
    
    # series = ['72_n100000_m12n8e3_nts5_a', '14_n100000_m12_dot3_a', '25_n100000_m12_dot5_a', 
    #           '13_n100000_m12_dot7_a', '32_n100000_m12_one_a']
    

    # series = ['46_n100000_5c_nts5_a', '13_n100000_2_dot3_a', '24_n100000_2_dot5_a', 
    #           '12_n100000_2_dot7_a', '31_n100000_2_one_a']
    
    
    outputlist = ['Output', 'Output2', 'EirOutput']
    
    mast_withshift_dic = {'Shot': '027205', 'multi_shift': multi_shift, 'shift_dic': shift_dic, 
                          'shift_filelist': shift, 'tail': tail, 'series': series,
                          'Output': outputlist}
    
    return mast_withshift_dic


def Ashift_dir_comp():
    multi_shift = ['MAST', 'D3D']
    # multi_shift = ['org', 'dot3', 'dot5', 'dot7', 'one']
    shift_dic = {'MAST': 0, 'D3D': 0.55}
    shift = ['org_new_series', 'AD3D']
    tail = {'org': 'nts5_a', 'D3D': 'd3d_a'}
    
    
    series = ['75_n900000_leakbtarnsol_nts5_a']
    
    outputlist = ['Output', 'Output2', 'EirOutput']
    
    mast_withshift_dic = {'Shot': '027205', 'multi_shift': multi_shift, 'shift_dic': shift_dic, 
                          'shift_filelist': shift, 'tail': tail, 'series': series,
                          'Output': outputlist}
    
    return mast_withshift_dic




def mast_comp_dir_series():
    a_shift = 'org'
    shift = 0
    tail = '_leakbsol_nts5_a'
    outputlist = ['Output', 'Output2', 'EirOutput']
    shift_filename = 'org_denscan_fluxb_027205'
    mast_series_dir_dic = {'Shot': '027205', 'shift': shift_filename, 'shift_value': shift,
                    'tail': tail, 'a_shift': a_shift, 'Output': outputlist}
    
    return mast_series_dir_dic


def mast_comp_dir_tempscan():
    a_shift = 'org'
    shift = 0
    tail = '_leakbsol_nts5_a'
    outputlist = ['Output', 'Output2', 'EirOutput']
    shift_filename = 'org_tescan_fluxb_027205'
    mast_series_dir_dic = {'Shot': '027205', 'shift': shift_filename, 'shift_value': shift,
                    'tail': tail, 'a_shift': a_shift, 'Output': outputlist}
    
    return mast_series_dir_dic

    
def mast_comp_dir_eireneN():
    a_shift = 'org'
    shift = 0
    tail = '_nts5_a'
    outputlist = ['Output', 'Output2', 'EirOutput']
    shift_filename = 'org_change_particle_number'
    mast_eireneN_dir_dic = {'Shot': '027205', 'shift': shift_filename, 'shift_value': shift,
                    'tail': tail, 'a_shift': a_shift, 'Output': outputlist}
    
    return mast_eireneN_dir_dic


def terminal_series_comp_dir(tail, filename):
    a_shift = 'org'
    shift = 0
    outputlist = ['Output', 'Output2', 'EirOutput']
    
    ds_dic = {'start': 5.02, 'stop': 9.02, 'space': 5}
    ts_dic = {'start': 3.73, 'stop': 7.73, 'space': 5}
    
    ds_list, ts_list = scan_list(denscan_dic = ds_dic, tempscan_dic = ts_dic)
    
    
    mast_eireneN_dir_dic = {'Shot': '027205', 'filename': filename, 'shift_value': shift,
                    'tail': tail, 'a_shift': a_shift, 'Output': outputlist, 
            'denscan_list': ds_list, 'tempscan_list': ts_list}
    
    return mast_eireneN_dir_dic




def mast_comp_dir_compare():
    
    a_shift = 'org'
    shift = 0
    
    multi_shift = ['org', 'dot3', 'dot5', 'dot7', 'one']
    shift_dic = {'org': 0, 'dot3': 0.3, 'dot5': 0.5, 'dot7': 0.7, 'one': 1}
    shift = ['org_new_series', 'dot3', 'dot5', 'dot7', 'one_LS']
    tail = {'org': 'nts_a', 'dot3': 'dot3_a', 'dot5': 'dot5_a', 'dot7': 'dot7_a',
            'one': 'one_a'}
    series = ['72_n100000_n5e3et1e2_nts5_a', '14_n100000_leakagebou_dot3_a', '25_n100000_leakagebou_dot5_a', 
              '13_n100000_leakagebou_dot7_a', '32_n100000_leakagebou_one_a']
    
    series_2 = ['72_n100000_m12n8e3_nts5_a', '14_n100000_m12_dot3_a', '25_n100000_m12_dot5_a', 
              '13_n100000_m12_dot7_a', '32_n100000_m12_one_a']
    
    outputlist = ['Output', 'Output2', 'EirOutput']
    shift_filename = 'org_new_series'
    mast_compare_dir_dic = {'Shot': '027205', 'shift': shift, 'shift_value': shift_dic,
                    'tail': tail, 'a_shift': a_shift, 'Output': outputlist}
    
    return mast_compare_dir_dic



def set_wdir(): #Function to set correct Working Directory Path depending on which machine is in use
    if os.environ['OS'] == 'Windows_NT':
        if os.environ['USERNAME'] == 'Yi-Cheng':
            basedrt = r"C:/Users/Yi-Cheng/Documents/SOLPS_Data/Simulation_Data"
            topdrt = r"C:/Users/Yi-Cheng/Documents/SOLPS_Data/Experimental_Data"


        elif os.environ['USERNAME'] == 'user':
            basedrt = r"C:/Users/user/Documents/SOLPS data/simulation data"
            topdrt = r"C:/Users/user/Documents/SOLPS data/experiment data"
            
            
        elif os.environ['USERNAME'] == 'ychuang':
            basedrt = r"C:/Users/ychuang/Documents/SOLPS_data/simulation_data"
            topdrt = r"C:/Users/ychuang/Documents/SOLPS_data/experimental_data"
           

    elif os.environ['OS'] == '5.14.0-362.24.1.el9_3.0.1.x86_64':
        if os.environ['USER'] == 'ychuang':
            basedrt = r"/sciclone/data10/ychuang/solps-iter/runs/mast"
            topdrt = r"/sciclone/data10/ychuang/solps-iter/runs/mast/gnpfiles"
           
    else:
        print('please add new directory in tools')
    
    return basedrt, topdrt





def set_figdir(): #Function to set correct Working Directory Path depending on which machine is in use
    if os.environ['OS'] == 'Windows_NT':
            
        if os.environ['USERNAME'] == 'ychuang':
            
            fig_dir = r"C:\Users\ychuang\Documents\SOLPS_data\simulation_data\mast\027205\dataplot_fig"
    
    return fig_dir





def s_number(text, series_flag):
    sd = Setting_dic()
    if sd['withshift'] == False and sd['withseries'] == False:
        name = text.split("/",-1)[-2]
        nu = int(name.split('_')[0])
    elif sd['withshift'] == False and sd['withseries'] == True:
        if series_flag == 'change_den':
            name = text.split("\\",-1)[-1]
            nu = re.findall('\d+\.\d+', name)
            nu.append(name.split('_')[0])
            # print(nu)
        elif series_flag == 'eireneN':
            name = text.split("\\",-1)[-1]
            nu = re.findall('\d+', name)
            nu.append(name.split('_')[0])
            # print(nu)
        elif series_flag == 'change_temp':
            name = text.split("\\",-1)[-1]
            nu = re.findall('\d+\.\d+', name)
            nu.append(name.split('_')[0])
        
        elif series_flag == 'two_compare':
            name = text.split("\\",-1)[-1]
            nu = re.findall('\d+', name)
            nu.append(name.split('_')[0])
            # print(nu)
        
        elif series_flag == 'twin_scan':
            name = text.split("/",-1)[-1]
            nu = re.findall('\d+\.\d+', name)
            nu.append(name.split('_')[0])
        
        else:
            print('check the series flag')
        
    elif sd['withshift'] == True and sd['withseries'] == False:
        name = text.split("/",-1)[-2]
        nu = int(name.split('_')[0])
    elif sd['withshift'] == True and sd['withseries'] == True:
        print('unexpected situation, please check the parameter setting')
    else:
        print('There is a bug in s_number function')

    return [nu, name]
        



def atp_number(text, series_flag):
    sd = Setting_dic()
    if sd['withshift'] == False and sd['withseries'] == False:
        name = text.split("/",-1)[-2]
        nu = int(name.split('_')[0])
    elif sd['withshift'] == False and sd['withseries'] == True:
      
        if series_flag == 'twin_scan':
            divider = "\\"
            if divider in text:
                # print('{} is in directory'.format(divider))
                name = text.split("\\",-1)[-1]
            else:
                name = text.split("/",-1)[-1]
                # print(name)
            
            nu_list = re.findall('\d+\.\d+', name)
            nu_tuple = (nu_list[0], nu_list[1])
            nu = [nu_tuple, name.split('_')[0]]
        
        else:
            print('check the series flag')
        
    elif sd['withshift'] == True and sd['withseries'] == False:
        name = text.split("/",-1)[-2]
        nu = int(name.split('_')[0])
    elif sd['withshift'] == True and sd['withseries'] == True:
        print('unexpected situation, please check the parameter setting')
    else:
        print('There is a bug in s_number function')

    return nu




def scan_list(denscan_dic, tempscan_dic):
    
    ds_start_num = denscan_dic['start']
    ds_stop_num = denscan_dic['stop']
    ds_space_num = denscan_dic['space']

    ds_list = np.linspace(ds_start_num, ds_stop_num, ds_space_num)

    ts_start_num = tempscan_dic['start']
    ts_stop_num = tempscan_dic['stop']
    ts_space_num = tempscan_dic['space']

    ts_list = np.linspace(ts_start_num, ts_stop_num, ts_space_num)
    
    return ds_list, ts_list
    
    


def loadDS_dic(DEV):
    "New DefaultSettings for loading experimental data"
    
    bload = {'TimeRange' : [1.10,1.30], 'AVG': False, 'ROOTSHOT': ''}
    
    # bload = {'TimeRange' : [1.10,1.30], 'AVG': False, 'multishift': False,
    #          'EXP': False, 'fit': True, 'ROOTSHOT': ''}
    if DEV == 'mast':
        fndic = {'expfilename': 'yag_27205_275.dat', 'fitfname': 'fit_027205_275.dat'}
    else:
        print('please add the experimental file name')
        
    loadDS = {**bload, **fndic} 
    
    return loadDS
    

A = ['39']




"""
REQUIRED:
    
    ** Setup **
    
    Shot = Full discharge number, only 1 per object instance\n 
    Attempts = Simulation run attempt ID(s) [As many as desired, as a List!]\n
    Parameters = Parameters to be plotted [Specified by a list of strings] - Defaults are:  
        'Ne' - Electron Density\n
        'Te' - Electron Temperature\n
        'Ti' - Ion Temperature\n
        'DN' - Particle Diffusion Coefficient\n
        'KYE' - Electron Thermal Diffusion Coefficient\n
        'KYI' - Ion Thermal Diffusion Coefficient  \n
        'NeuDen' - Neutral Density\n
        'IonFlx' - Radial Ion Particle Flux\n
        'IonPol' - Poloidal Ion Particle Flux\n
        'RadPinch' - Radial Pinch Velocity
        
OPTIONAL:

    TimeRange = [1.10,1.30] > Time range (in sec) over which experimental data is averaged\n    
    DEV = 'cmod' > Tokamak device being investigated\n
    EXP = True > Plot Experimental Data Points\n
    LOG10 = 0 > Plot Base-10 Logarithm of Parameter Data (0, 1, or 2)
        0 : No log plot or log processing of data\n
        1 : Take base-10 log of data, then plot on linear y-scale (ticks will be 0, 1, 2,...)\n
        2 : Plot data as is, then change y-axis scale to logarithmic (ticks will be 0, 10, 100,...)
    GRAD = False > Calculate Gradient of Parameter Data w.r.t. Radial Coordinate specified by RADC\n
    ELEV = 75 > Elevation of viewing angle for Surface Plot\n
    AZIM = 270 > Azimuthal viewing angle for Surface Plot\n
    JXI = 37 > Poloidal position of Inner Midplane - default is 37\n
    JXA = 55 > Poloidal position of Outer Midplane - default is 55\n
    SEP = 20 > Radial position of Separatrix - default is 20\n
    XDIM = 98 > Dimension of computational grid in the x (poloidal) direction\n
    YDIM = 38 > Dimension of computational grid in the y (radial) direction\n
    CoreBound = [24,71] > X-coordinate grid cell numbers that define the [left, right] bounds of Core Region\n
    CMAP = cm.viridis > Colormap used for 2D Contour Maps\n
    Colors = None > Specific colors to color each contour level of a contour map\n
    Publish = [] > List of strings to use in legends of publication-quality plots; if not [], changes plotting rc.params\n
    Markers = True > Plot separatrix, outer midplane, and/or inner midplane marker line\n
    PlotScheme = [] > List of matplotlib plot style flags for line plots; must either be empty or a list of strings as long as the number of Attempts\n    
    PsinOffset = 0 > Psi_n offset of Experimental Data Points\n
    RadOffset = 0 > Radial (in meters) offset of Experimental Data Points\n
    RADC = 'psin' > Set radial coordinate convention - Either 'psin', 'radial', 'rrsep' or 'Y'\n    
    POLC = 'dXP' > Set poloidal coordinate convention - Either 'X', 'theta', 'dXP', or 'dXP_norm'\n
    RadSlc = None > Radial surface selection for poloidal plots - Can set specific radial index, 'all', or 'None' defaults to SEP\n
    PolSlc = None > Poloidal grid line selection for radial plots - Can set specific poloidal index, 'all', or 'None' defaults to JXA\n
    PHYS = True > Map Contour to PHYSical Geometry; if False, plots on rectangular grid\n     
    LVN = 100 > Number of colorbar levels for contour plots or List [] of specific contour levels\n
    LINTHRESH = 1.0 > +/- Value at which divergent-logarithmic Contour plots should switch to linear scale going from positive to negative values\n
    DIVREG = True > Include Divertor Region Data (May cause loss of logarithmic resolution)\n   
    SAVE = False > Save plot to a .png file\n 
    DIFFERENCE = False > For a series of Attempts, take the difference between the matrix data of each Attempt from the first declared Attempt\n
    AVG = False > Take Average of each parameter dataset over all supplied Attempts, append averaged data to end of matrix dataset\n
    TC_Flux = [] > For d3d cases only; allows specification of ONETWO-calculated Radial Flux values\n
    TC_Psin = [] > For d3d cases only; allows specification of Psin coordinates of ONETWO-calculated Radial Flux values\n
    GRID = False > Turn on plot grid\n
    AX = None > Pass the name of a matplotlib axis for the SOLPSPLOT object to plot on; by default SOLPSPLOT plots on a new axis\n
    BASEDRT = 'SOLPS_2D_prof/' > Local base directory for all stored SOLPS simulation run data\n
    TOPDRT = '' > Local home directory, parent of BASEDRT and 'gfileProcessing' directories\n
    ROOTSHOT = '' > Optional numerical prefix to designate experimental campaign that a series of shots may share in common\n 
"""

P = {'Ne': r'Electron Density $n_e\;(m^{-3})$',
        'Te': r'Electron Temperature $T_e\;(eV)$',
        'NI' : r'Ion (D+) Density $n_i\;(m^{-3})$',
        'Ti': r'Ion (D+) Temperature $T_i\;(eV)$',
        'DN': r'Particle Density Diffusivity $D\;(m^2/s)$',
        'KYE': r'Electron Thermal Diffusivity $\chi_e (m^2/s)$',
        'KYI': r'Ion Thermal Diffusivity $\chi_i (m^2/s)$',
        'NeuDen': r'Neutral Atom (D) Density $(m^{-3})$',
        'MolDen': r'Neutral Molecule (D2) Density $(m^{-3})$',
        'NeuTemp': r'Neutral Atom (D) Temperature (eV)',
        'MolTemp': r'Neutral Molecule (D2) Temperature (eV)',
        'IonFlx': r'Radial Ion (D+) Flux $s^{-1}$',
        'IonPol': r'Poloidal Ion (D+) Flux $s^{-1}$',
        'MolFlx': r'Radial Molecule Flux $m^{-2}s^{-1}$',
        'RadFlx': r'Radial Atomic Flux $m^{-2}s^{-1}$',
        'VLY': r'Radial Pinch $v_y (m^2/s)$',
        'SX': r'Poloidal Contact Area SX $(m^{2})$',
        'SY': r'Radial Contact Area SY $(m^{2})$',
        'SZ': r'Poloidal Cross-Sectional Area SZ $(m^{2})$',
        'VOL': r'Cell Volume VOL $(m^{3})$',
        'RadPinch': r'Anomalous Radial Pinch Velocity $v_y\;(m/s)$',
        'AHAL': r'Atomic $H_\alpha$ emissivity $(photons*m^{-3}s^{-1})$',
        'MHAL': r'Molecular $H_\alpha$ emissivity $(photons*m^{-3}s^{-1})$',
        'NTI' : r'Test Ion (D2+) Density $n_{ti}\;(m^{-3})$',
        'TestIonTemp' : r'Test Ion (D2+) Temperature $T_{testion}\;(eV)$',
        'LyaEmiss' : r'Lyman-alpha Emissivity $(photons*m^{-3}s^{-1}sr^{-1})$',
        'LyaEmissW' : r'Converted Lyman-alpha Emissivity $(W*m^{-3})$',
        'LyaEmissW_IE' : r'Ion-Electron Reaction Lyman-alpha Emissivity $(W*m^{-3})$'}

DP_backup = {'LOG10' : 0,
        'GRAD' : False,
        'ELEV' : 75,
        'AZIM' : 270,
        'JXI' : 37,
        'JXA' : 57,
        'SEP' : 20,
        'XDIM' : 98,
        'YDIM' : 38,
        'CoreBound' : [24,71],
        'Publish' : [],
        'Markers' : True,
        'PlotScheme' : [],
        'PsinOffset' : 0,
        'RadOffset' : 0,
        'RADC' : 'psin',
        'POLC' : 'dXP',
        'RadSlc' : None,
        'PolSlc' : None,
        'SURF' : 20,
        'GEO' : True,
        'LVN' : 100,
        'DIVREG' : True,
        'SAVE' : False,
        'SUBTRACT' : False,
        'AVG' : False,
        'TC_Flux' : [], 
        'TC_Psin' : [],
        'GRID': False,
        'AX' : None} #1160718

# for key, value in DP.items():
#     print(key)
