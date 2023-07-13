# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 12:38:48 2023

@author: user
"""

def Setting_dic():
    setting_dic = {'Shot': '27205', 'Attempt': A, 'Parameters': P, 
                   'DefaultSettings': DP}
    return setting_dic



A = ['39']

P = ['Ne','Te','Ti','DN','KYE','KYI','NeuDen','IonFlx','IonPol','RadPinch']

DP = {'TimeRange' : [1.10,1.30],  
                 'DEV': 'cmod',
                 'EXP' : True,
                 'LOG10' : 0,
                 'GRAD' : False,
                 'ELEV' : 75,
                 'AZIM' : 270,
                 'JXI' : 37,
                 'JXA' : 55,
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
                 'AX' : None,
                 'BASEDRT': 'solps-iter/runs/',
                 'TOPDRT' : '',
                 'ROOTSHOT' : ''} #1160718

for key, value in DP.items():
    print(key)
