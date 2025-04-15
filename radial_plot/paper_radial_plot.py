# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 19:11:41 2025

@author: ychuang
"""

from SOLPS_input.header import *


# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.offsetbox import AnchoredText


class paper_radial:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data




    def paper_neuden_radial_method(self, fit_dat, x_coord, Nd, log_flag, itername):
            
        exp_an_fit = fit_dat['exp_fit']
        xcoord_cut = fit_dat['x_coord_cut']
        
        fig, axs = plt.subplots()
        
        anchored_text = AnchoredText('{}'.format('Neutral density [$m^{-3}$] with fits'), loc='upper left')
        if log_flag:
            axs.set_yscale('log')
        else:
            pass
        axs.plot(x_coord, Nd,'-', color = 'green', label= 'solps neutral density [$m^{-3}$]')
        axs.plot(xcoord_cut, exp_an_fit, color='r',lw= 5, ls='-', label= 'exponential fit')
        axs.axvline(x= max(xcoord_cut), color='black',lw=3, ls='--', 
                    label= 'fit range : $\Delta n_e$')
        axs.axvline(x= min(xcoord_cut), color='black',lw=3, ls='--')
        axs.set_xlabel('$\psi_N$')
        axs.add_artist(anchored_text)
        axs.legend(loc= 'lower right')
        
        
        
        # fig.savefig('neuden_fit_{}.pdf'.format(itername))
        
    
    
    def paper_neuden_radial_plot(self, pol_loc, dat_size):
        
        
        
        if self.withshift == True and self.withseries == False:
            
            
            for aa in self.data['dircomp']['multi_shift']:
                
                
                nx = self.data['b2fgeo'][aa]['nx']
                ny = self.data['b2fgeo'][aa]['ny']
                
                if dat_size == 'full':
            
                    dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
                
                elif dat_size == 'small':
                    dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
                
                
                psi_list, nd_list, fit_dat = self.neuden_midprof(iter_index = aa, 
                                                        data_struc = dat_struc)
                
    
                self.paper_neuden_radial_method(fit_dat= fit_dat, itername = aa,
                    x_coord = psi_list, Nd = nd_list, log_flag = True)
        
        else:
            print('this is the modify major radius version')