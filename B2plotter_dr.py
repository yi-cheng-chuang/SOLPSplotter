# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 12:45:19 2023

@author: user
"""
import B2plotter_set as b2s
import SOLPSplotter_geoRR as sgr
import B2plotter_plot as b2p
import B2plotter_mid_calc as b2mc

d = b2s.Setting_dic()
lex = b2s.loadDS_dic(d['DEV'])




xl = b2p.Opacity_study(DEV = d['DEV'], withshift= d['withshift'], withseries= d['withseries'],
            DefaultSettings = d['DefaultSettings'], loadDS = lex, 
            Parameters= d['Parameters'], Publish= d['Publish'])
xl.load_mast_dir()
xl.load_solpsgeo()
xl.calc_RRsep()
# xl.calcpsi()
# xl.calc_RRsep()
xl.calcpsi_1D(pol_loc= '27')
# xl.calc_flux_expansion(polpos= '55')
# xl.load_vessel()
xl.set_plot()
# xl.Opacity_study_radial_plot_dsa(pol_loc= '57')
xl.Opacity_study_radial_plot(pol_loc= '27', x_choice= 'psiN')


# poloidal_index_list = []
# for i in range(38):
#     poloidal_index_list.append('{}'.format(27 + i))
# xl.calc_pol_angle(pol_list = poloidal_index_list)
# xl.Opacity_study_poloidal_plot(pol_list= poloidal_index_list, x_choice= 'psiN')


# p = xl.data['dircomp']['Attempt']
# q = xl.data['gfile']['gcomp']['check']
# print(p)
# print(q)




"Use functions Richard.R created"
# xd = sgr.geo_RR(DEV = d['DEV'], withshift= d['withshift'], withseries= d['withseries'],
#             DefaultSettings = d['DefaultSettings'])
# xd.load_mast_dir()
# xd.load_solpsgeo()
# gd = xd.creat_grid()
# xd.load_output_geo(grid_dic= gd)
# xd.dsa_1D(pol_loc= '46')

"Use functions for mid_plane calculation"
# xt = b2mc.psi_RRsep_calc(DEV = d['DEV'], withshift= d['withshift'], withseries= d['withseries'],
#             DefaultSettings = d['DefaultSettings'])

# xt.load_mast_dir()
# xt.load_solpsgeo()
# xt.calcpsi()
# # xt.plot_psi_surface()
# xt.plot_seperatrix()

# poloidal_index_list = []
# for i in range(30):
#     poloidal_index_list.append('{}'.format(30 + i))
# xt.plot_psi(pol_list= poloidal_index_list)
# xt.plot_RZ(pol_list= poloidal_index_list)





