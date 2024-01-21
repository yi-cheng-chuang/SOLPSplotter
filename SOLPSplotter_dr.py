# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 12:45:19 2023

@author: user
"""
import SOLPSplotter_set as sps
import B2plotter_contour as b2c
import SOLPSplotter_transcoe_adj as sta

d = sps.Setting_dic()
lex = sps.loadDS_dic(d['DEV'])



# xt = sta.transport_coefficient_adjustment(DEV = d['DEV'], withshift= d['withshift'], 
#         withseries= d['withseries'], DefaultSettings = d['DefaultSettings'])
# xt.load_mast_dir()
# xt.load_solpsgeo()
# xt.calcpsi()
# xt.mod_transco(withmod = False, de_SOL = 26, ki_SOL = 31, ke_SOL = 25, log_flag = False)



xl = b2c.PlotContour(DEV = d['DEV'], withshift= d['withshift'], withseries= d['withseries'],
            DefaultSettings = d['DefaultSettings'], loadDS = lex, 
            Parameters= d['Parameters'], Publish= d['Publish'])
xl.set_plot()
xl.load_mast_dir()
xl.load_solpsgeo()
xl.calcpsi()
xl.calc_RRsep(plotRR= False, plot_psi_dsa_align= False)
fitmastexp_setting_dic = {'writefile': True, 'plot_solps_fit': False, 
                          'plot_exp_and_fit': True, 'plot_shift_compare': False,
                          'data_print': True}
xl.fitmastexp(plot_setting_dic = fitmastexp_setting_dic)
xl.transport_coe_align_plot(plot_transcoe= False)
xl.load_b2fstate()
xl.load_vessel()
# xl.flux_comparison_plot()



# xl.calcpsi()

plot_flag = 'skip'
if plot_flag == 'radial':
    # xl.set_plot()
    PL = '59'
    xl.calcpsi_1D(pol_loc= PL)
    xl.calc_dsa(pol_loc= PL)
    xl.flux_expansion_map(pol_loc= PL, iter_index= None)
    xl.Opacity_study_radial_plot(pol_loc= PL)
elif plot_flag == 'poloidal':
    # xl.set_plot()
    poloidal_index_list = []
    for i in range(40):
        poloidal_index_list.append('{}'.format(28 + i))
    xl.calc_pol_angle(pol_list = poloidal_index_list, plot_angle= True)
    xl.Opacity_study_poloidal_plot(pol_list= poloidal_index_list)
elif plot_flag == 'skip':
    pass
else:
    print('check plot_flag')


radial_plot_flag = True
if radial_plot_flag:
    xl.plot_all_radial()
else:
    pass
        









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





