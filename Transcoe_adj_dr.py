# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 20:13:25 2024

@author: user
"""



import SOLPSplotter_set as sps
import B2plotter_contour as b2c
import SOLPSplotter_transcoe_adj as sta

d = sps.Setting_dic()

xt = sta.transport_coefficient_adjustment(DEV = d['DEV'], withshift= d['withshift'], 
        withseries= d['withseries'], DefaultSettings = d['DefaultSettings'])
xt.load_mast_dir()
xt.load_solpsgeo()
xt.calcpsi()
xt.mod_transco(withmod = True, de_SOL = 24, ki_SOL = 31, ke_SOL = 23, log_flag = False)