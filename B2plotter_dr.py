# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 12:45:19 2023

@author: user
"""
import B2plotter_set as b2s
import SOLPSplotter_geoRR as sgr
import B2plotter_plot as b2p
import B2plotter_flux as b2f
import B2plotter_mid_calc as b2mc

d = b2s.Setting_dic()
lex = b2s.loadDS_dic(d['DEV'])




xl = b2f.flux_adjustment(DEV = d['DEV'], withshift= d['withshift'], withseries= d['withseries'],
            DefaultSettings = d['DefaultSettings'], loadDS = lex, 
            Parameters= d['Parameters'], Publish= d['Publish'])
xl.set_plot()
xl.load_mast_dir()
xl.load_solpsgeo()
xl.calc_RRsep(plotRR= False)
xl.fitmastexp(writefile= False)
xl.transport_coe_align_plot()
xl.load_vessel()
# xl.flux_comparison_plot()

# xl.calcpsi()

plot_flag = 'radial'
if plot_flag == 'radial':
    # xl.set_plot()
    PL = '59'
    xl.calcpsi_1D(pol_loc= PL)
    xl.calc_dsa(pol_loc= PL)
    xl.Opacity_study_radial_plot(pol_loc= PL, x_choice= 'psiN')
elif plot_flag == 'poloidal':
    # xl.set_plot()
    poloidal_index_list = []
    for i in range(47):
        poloidal_index_list.append('{}'.format(25 + i))
    xl.calc_pol_angle(pol_list = poloidal_index_list)
    xl.Opacity_study_poloidal_plot(pol_list= poloidal_index_list, x_choice= 'psiN')
elif plot_flag == 'skip':
    pass
else:
    print('check plot_flag')


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





