# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 19:11:04 2025

@author: ychuang
"""


from SOLPS_input.header import *


class fits_on_midplane:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    def neuden_midprof(self, iter_index, data_struc):
        
        
        
        cut_st = result_dic['59']['x_coord_cut']
        w_cut = []
        
        
        for ic, coord in enumerate(psi_val):
            
            if coord >= cut_st.min() and coord <= cut_st.max():
                w_cut.append(weight[ic])
        
        wcut_B = np.ones(len(w_cut))- w_cut
        
        fit_dat = {}
        
        expfit_a = result_dic['58']['exp_fit']
        expfit_b = result_dic['60']['exp_fit']
        
        fit_dat['exp_fit'] = np.multiply(expfit_a, w_cut) + np.multiply(expfit_b, wcut_B)
        
        cut_a = result_dic['58']['x_coord_cut']
        cut_b = result_dic['60']['x_coord_cut']
        
        fit_dat['x_coord_cut'] = np.multiply(cut_a, w_cut) + np.multiply(cut_b, wcut_B)
        
        return psi_list, nd_list, fit_dat


