# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 20:13:25 2024

@author: user
"""



import SOLPS_set as sps
import SOLPS_transcoe_adj as sta

d = sps.Setting_dic()

xt = sta.transport_coefficient_adjustment(DefaultSettings = d)
xt.load_mast_dir()
xt.load_solpsgeo()
xt.calcpsi_avcr()
xt.set_plot()
xt.mod_transco(withmod = False, de_SOL = 24, ki_SOL = 23, ke_SOL = 23, log_flag = False)
xt.transport_coe_align_plot(plot_transcoe = True, paper_transcoe = True, save_eps = True)
xt.align_transco(plot_align = True, log_flag = False)


extra = False


if extra:
    
    basedrt, topdrt = sps.set_wdir()

    FL = basedrt + '/mast/027205/org_cfluxb_std/80_lfbcheck_1e3_a'
    FL2 = basedrt + '/mast/027205/org_cfluxb_std/80_nf5.512tf4.115_fast_a'
    # xt.calcpsi_block_method(file_loc = FL, shift = 0)
    flist = [FL, FL2]

    xt.transport_coe_compare_plot(file_loc_list = flist, plot_compare = True)





