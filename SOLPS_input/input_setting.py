# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 12:38:48 2023

@author: user
"""

import os

def Setting_dic():
    
    if os.environ['OS'] == 'Windows_NT':
        terminal = False
    elif os.environ['OS'] == '5.14.0-362.24.1.el9_3.0.1.x86_64':
        terminal = True
    else:
        print('there is a bug at Setting_dic function or unrecognized operating system')
    
    
    class DefaultSettings:
        
        def __init__(self):
            
            
            self.DEV = 'mast'
            self.a = 0.5
            self.withshift = False
            self.withseries = True
            self.terminal = terminal
            if self.withseries:
                self.series_flag = 'twin_scan'
                self.series_filename = 'org_new25scan_fast_save'
                self.series_tail = 'fast_a'
                
            else:
                print('withseries is false!')
                
            
            self.series_compare = False
            self.plot_setting = 'mod_transcoe'
            self.data_size = 'small'
    
    DF = DefaultSettings()
    
    # set_dic = {'DEV': 'mastu', 'minor_rad': 0.5, 'withshift': False, 'withseries': False,
    #             'Parameters': P, 'series_flag': 'twin_scan', 'series_compare': False,
    # 'series_filename': 'org_new25scan_fast_save', 'series_tail': '_fast_a',
    #                 'terminal': terminal}
    
    return DF

"""
series_flag = ['eireneN','change_den','change_temp', 'two_compare', 'twin_scan']
plot_setting = ['mod_transcoe']

"""


def set_figdir(): #Function to set correct Working Directory Path depending on which machine is in use
    if os.environ['OS'] == 'Windows_NT':
            
        if os.environ['USERNAME'] == 'ychuang':
            
            fig_dir = r"C:\Users\ychuang\Documents\SOLPS_data\simulation_data\mast\027205\dataplot_fig"
    
    return fig_dir
    
    


def loadDS_dic(DEV):
    "New DefaultSettings for loading experimental data"
    
    bload = {'TimeRange' : [1.10,1.30], 'AVG': False, 'ROOTSHOT': ''}
    
    # bload = {'TimeRange' : [1.10,1.30], 'AVG': False, 'multishift': False,
    #          'EXP': False, 'fit': True, 'ROOTSHOT': ''}
    if DEV == 'mast':
        fndic = {'expfilename': 'yag_27205_275.dat', 'fitfname': 'fit_027205_275.dat'}
        loadDS = {**bload, **fndic}
        
        
    else:
        print('please add the experimental file name')
        loadDS = {'Notice': 'loadDS dic is empty, please add the experimental file name.'}
    
    
    
    return loadDS

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
