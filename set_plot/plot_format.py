# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 00:07:50 2025

@author: ychuang
"""


import matplotlib.pyplot as plt


class plot_setting():
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    
    def set_plot(self):
        
        if self.DF.plot_setting == 'mod_transcoe':
            
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 5, markersize= 9)
            plt.rcParams.update({'font.size': 16})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
        
        

        elif self.DF.plot_setting == 'simple_radial':
            
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 4, markersize= 6)
            plt.rcParams.update({'font.size': 10})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
            # plt.rcParams["text.usetex"] = True

        

        elif self.DF.plot_setting == 'aspect_ratio_correlate':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 5, markersize= 9)
            plt.rcParams.update({'font.size': 16})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
        
        
        elif self.DF.plot_setting == 'contour_plot':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 3, markersize= 6)
            plt.rcParams.update({'font.size': 10})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
        
        
        elif self.DF.plot_setting == 'pol_subplot':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 4, markersize= 7)
            plt.rcParams.update({'font.size': 14})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
            # plt.rcParams['figure.figsize'] = 40, 12
        
        

        elif self.Publish == 'b2plottersetting':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 5, markersize= 9)
            plt.rcParams.update({'font.size': 16})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
  

           
 
        else:
            print('Publish setting is incorrect or add another setting')



    
    
    def radial_import_Plabel(self):
        
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
        
        self.data['Parameter'] = P
    
    
    def opacity_study_unit(self):
        unit = {'efold_length_psiN': 'Neutral penetration length ($\psi_N$)',
                'pedestal_width_psiN': 'Pedestal width ($\psi_N$)',
                  'dimensionless_opaqueness': 'Modeled opaqueness', 
                  'neutral_density': 'Neutral density ${n_D}$ (m$^{-3}$)', 
                  'electron_pedestal_density': 'Electron pedestal density: $n_{ped}$ (m$^{-3}$)',
                  'temperature_pedestal_width': 'Temperature pedestal width: $\Delta T$',
                  'flux_expansion': 'Flux expansion',
                  'efold_length': '$\lambda_{n_D}$ [mm]',
                  'pedestal_width': '$\Delta n_e$ [mm]',
                  
                  }
        
        self.data['opacity_study_unit'] = unit
    
    
    
    

        


