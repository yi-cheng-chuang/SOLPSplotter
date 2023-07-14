# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 18:39:23 2023

@author: user
"""

if 'd3d' in Shot or DEV=='d3d':
    
    BASEDRT = '{}d3d'.format(BASEDRT)
    
    JXI = 39
    JXA = 55
    SEP = 17
    
    GFILE = '{}gfileProcessing/d3d_files/g175060.02512'.format(TOPDRT)
    GF = eq.equilibrium(gfile=GFILE)
    
    ExpData = loadmat('{}gfileProcessing/d3d_files/175060_data_SOLPS.mat'.format(TOPDRT))
    ONETWO = loadmat('{}gfileProcessing/d3d_files/flow_transport.mat'.format(TOPDRT))
    
    ii = 0
    jj = 0
    kk = 0
    
    Z_mid = 0
    R_OMcore = 2.176
    psin_core = GF.psiN(R_OMcore,Z_mid)
    
    while ExpData['psin_ne'][ii] < psin_core:
        ii += 1
    while ExpData['psin_te'][jj] < psin_core:
        jj += 1 
    while ExpData['psin_ti'][kk] < psin_core:
        kk += 1
    
    PsinNe = ExpData['psin_ne'][ii:]
    PsinTe = ExpData['psin_te'][jj:]
    PsinTi = ExpData['psin_ti'][kk:]
    PsinAvg = [PsinNe,PsinTe,PsinTi]
    
    self.ExpDict['Ned3d'] = ExpData['Ne'][ii:]
    self.ExpDict['Ted3d'] = ExpData['Te'][jj:]
    self.ExpDict['Tid3d'] = ExpData['Ti'][kk:]