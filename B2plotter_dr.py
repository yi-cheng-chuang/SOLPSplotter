# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 12:45:19 2023

@author: user
"""
import B2plotter_set as b2s
import B2plotter_class as b2c
import B2plotter_plot as b2p

d = b2s.Setting_dic()
lex = b2s.loadDS_dic(d['DEV'])

xl = b2p.simple_plot(DEV = d['DEV'], withshift= d['withshift'], DefaultSettings = d['DefaultSettings'],
                   loadDS = lex, Parameters= d['Parameters'], 
                   Publish= d['Publish'])
xl.load_mast_dir()
xl.load_solpsgeo()
# xl.calcpsi()
xl.calcpsi_1D(pol_loc= '57')
# xl.load_vessel()
# gd = xl.creat_grid()
# xl.load_output_geo(grid_dic= gd)
xl.plot_Ne_NeuDen_withshift(pol_loc= '57')
# p = xl.data['dircomp']['Attempt']
q = xl.data['gfile']['gcomp']['check']
# print(p)
print(q)



