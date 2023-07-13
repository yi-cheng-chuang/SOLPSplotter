# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 12:45:19 2023

@author: user
"""
import B2plotter_set as b2s
import New_plotter as b2p

d = b2s.Setting_dic()


xp = b2p.B2plotter(Shot= d['Shot'],Attempt= d['Attempt'], 
                   Parameters= d['Parameters'], 
                   DefaultSettings= d['DefaultSettings']).DefaultSettings
print(xp)