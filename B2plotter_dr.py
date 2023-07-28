# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 12:45:19 2023

@author: user
"""
import B2plotter_set as b2s
import B2plotter_class as b2c

d = b2s.Setting_dic()
lex = b2s.loadDS_dic(d['DEV'])

xl = b2c.load_data(DEV = d['DEV'], DefaultSettings = d['DefaultSettings']
                   , loadDS = lex, Parameters= d['Parameters'])
xl.load_mast_dir()
xl.load_solpsgeo()
xl.Calcpsi()
xl.load_vessel()
gd = xl.creat_grid()
xl.load_output_geo(grid_dic= gd)
xl.load_output_data('NeuDen')
p = xl.data['dircomp']['Attempt']
q = xl.data['gfile']['gcomp']['check']
print(p)
print(q)



