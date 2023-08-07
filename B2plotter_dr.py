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




xl = b2p.simple_plot(DEV = d['DEV'], withshift= d['withshift'], withseries= d['withseries'],
            DefaultSettings = d['DefaultSettings'], loadDS = lex, 
            Parameters= d['Parameters'], Publish= d['Publish'])
xl.load_mast_dir()
xl.load_solpsgeo()
# xl.calcpsi()
# xl.calcpsi_1D(pol_loc= '40')
# xl.load_vessel()
# gd = xl.creat_grid()
# xl.load_output_geo(grid_dic= gd)
xl.set_plot()
# xl.Opacity_study_radial_plot_psi(pol_loc= '40')
poloidal_index_list = []
for i in range(45):
    poloidal_index_list.append('{}'.format(26 + i))
xl.Opacity_study_poloidal_plot(pol_list= poloidal_index_list)
# p = xl.data['dircomp']['Attempt']
# q = xl.data['gfile']['gcomp']['check']
# print(p)
# print(q)



