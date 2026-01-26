# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 17:38:04 2026

@author: ychuang
"""


import numpy as np
import statistics
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib.pylab as pylab


class series_1para_HLcompare:

    def __init__(self, DF, data):

        self.DF = DF
        self.data = data

    def HLcompare_1para_plot_method(self, pol_list):

        midplane_psi = self.data['midplane_calc']['psi_solps_mid']
        r_rsep = self.data['midplane_calc']['R_Rsep']
        # xcoord_cut = self.data['radial_fit_data']['1.0']['x_coord_cut']

        psi_to_dsa_func = interpolate.interp1d(
            midplane_psi, r_rsep, fill_value='extrapolate')

        fig, axs = plt.subplots()

        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        dat_struc = {'nx': nx, 'ny': ny}

        HL_list = ['1.0', '4.5']

        marker_high = ['o', 's']
        lins_map_high = dict(zip(HL_list, marker_high))
        marker_low = ['v', '*']
        lins_map_low = dict(zip(HL_list, marker_low))

        colors_high = ['red', 'orange']
        color_map_high = dict(zip(HL_list, colors_high))
        colors_low = ['navy', 'purple']
        color_map_low = dict(zip(HL_list, colors_low))

        for aa in HL_list:

            """
            label= 'core density {} $10^{19}$'.format(aa)

            """
            axs.legend(loc='lower left', fontsize=10)

            def return_increase(high_dat, std_dat):

                norm_inc_dat = (high_dat - std_dat) * 100 / std_dat

                return norm_inc_dat

            fnay = self.data['b2wdat'][aa]['b2npc_fnays'][0][1:97, 1:37]
            hz = self.data['b2wdat'][aa]['hz'][1:nx+1, 1:ny+1]
            hy = self.data['b2wdat'][aa]['hy'][1:nx+1, 1:ny+1]
            hx = self.data['b2wdat'][aa]['hx'][1:nx+1, 1:ny+1]
            tor_area = np.multiply(hz, hy)
            pol_area = np.multiply(hz, hx)
            fnays = np.divide(fnay, pol_area)

            s_term = self.data['b2wdat'][aa]['b2npc_sna'][0][1:nx+1, 1:ny+1]
            vol = self.data['b2wdat'][aa]['vol'][1:nx+1, 1:ny+1]
            source = np.divide(s_term, vol)

            mid_fnays = fnays[59, :]
            highmid_fnays = fnays[40, :]

            arc_map = self.data['arc_map']

            print("inner midplane")
            print(pol_list[8])
            print("outer midplane")
            print(pol_list[30])

            xcoord_high = arc_map[8, :]
            xcoord_low = arc_map[30, :]

            dat_name = 'source'
            dat = source

            dat_high = dat[int(pol_list[8]), :]
            dat_low = dat[int(pol_list[30]), :]

            axs.plot(xcoord_high, dat_high, lins_map_high[aa] + '-', lw=1.5,
                     color=color_map_high[aa], label='{}HFS'.format(aa))

            axs.plot(xcoord_low, dat_low, lins_map_low[aa] + '-', lw=1.5,
                     color=color_map_low[aa], label='{}LFS'.format(aa))

        axs.set_xlabel(r'radial arc length $(r-r_{sep}) \; [m]$')

        if dat_name == 'source':

            axs.set_xlim(-0.02, 0.01)
            axs.set_yscale("log")
            # axs.set_ylim(4E+21, 3E+23)

        axs.grid(True)
        axs.legend()

    def HLcompare_1para_plot(self, pol_list):

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == False and withseries == True:

            self.HLcompare_1para_plot_method(pol_list=pol_list)
