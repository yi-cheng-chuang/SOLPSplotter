# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 05:56:01 2025

@author: ychuang
"""


import numpy as np
import statistics
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib.pylab as pylab


class series_midS:

    def __init__(self, DF, data):

        self.DF = DF
        self.data = data

    def Splot_method(self, iterlist, cl_dic, xcoord_type):

        midplane_psi = self.data['midplane_calc']['psi_solps_mid']
        r_rsep = self.data['midplane_calc']['R_Rsep']
        xcoord_cut = self.data['radial_fit_data']['1.0']['x_coord_cut']

        psi_to_dsa_func = interpolate.interp1d(
            midplane_psi, r_rsep, fill_value='extrapolate')
        marker = ['o', 's', 'x', '^', 'v', '*', '.', 'o', 's']
        lins_map = dict(zip(iterlist, marker))

        fig, axs = plt.subplots()

        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        dat_struc = {'nx': nx, 'ny': ny}

        avg_te_list = []

        for aa in iterlist:

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

            mid_fnays = fnays[59, :]

            mid_ne_pro = self.data['midplane_profile'][aa]['mid_ne']
            mid_te_pro = self.data['midplane_profile'][aa]['mid_te']
            mid_ti_pro = self.data['midplane_profile'][aa]['mid_ti']

            midte_trim = mid_te_pro[:33]
            midti_trim = mid_ti_pro[:33]

            mid_tei_pro = midte_trim / midti_trim
            mid_tie_pro = midti_trim / midte_trim

            ne_std = self.data['midplane_profile']['1.0']['mid_ne']

            ne_norm_inc_dat = return_increase(
                high_dat=mid_ne_pro, std_dat=ne_std)

            psi_coord = self.data['midplane_profile'][aa]['psiN']

            mid_neu_pro = self.data['midplane_profile'][aa]['mid_nd']
            insep_nd = self.data['midplane_profile'][aa]['mid_nd'][18]
            outsep_nd = self.data['midplane_profile'][aa]['mid_nd'][19]
            sep_nd = statistics.mean([insep_nd, outsep_nd])

            # max_nd = mid_neu_pro.max()
            nd_std = self.data['midplane_profile']['1.0']['mid_nd']

            nd_norm_inc_dat = return_increase(
                high_dat=mid_neu_pro, std_dat=nd_std)

            nd_change_percent = np.divide(
                mid_neu_pro, nd_std, out=np.zeros_like(mid_neu_pro, dtype=float), where=nd_std != 0) * 100

            norm_nd = mid_neu_pro / sep_nd

            midS = self.data['midplane_profile'][aa]['mid_S']
            max_S = max(self.data['midplane_profile'][aa]['mid_S'][:35])
            source_std = self.data['midplane_profile']['1.0']['mid_S']
            norm_S = midS / max_S
            source_change_percent = np.divide(
                midS, source_std, out=np.zeros_like(midS, dtype=float), where=source_std != 0) * 100

            S_norm_inc_dat = return_increase(
                high_dat=midS, std_dat=source_std)

            molsource = self.data['ft44'][aa]['srcml'][59, :, 0]

            rrsep_solps = psi_to_dsa_func(psi_coord)

            # avg_te =
            # avg_te_list.append()

            trim = False
            dat = S_norm_inc_dat

            if xcoord_type == 'psi':

                if trim:

                    axs.plot(psi_coord[:33], dat, lins_map[aa] + '-', lw=1.5,
                             color=cl_dic[aa], label='{}'.format(aa))

                else:

                    axs.plot(psi_coord, dat, lins_map[aa] + '-', lw=1.5,
                             color=cl_dic[aa], label='{}'.format(aa))

            elif xcoord_type == 'rrsep':

                if trim:

                    axs.plot(rrsep_solps[:33], dat, lins_map[aa] + '-', lw=1.5, color=cl_dic[aa],
                             label='{}'.format(aa))

                else:

                    axs.plot(rrsep_solps, dat, lins_map[aa] + '-', lw=1.5, color=cl_dic[aa],
                             label='{}'.format(aa))

        # axs.set_title('midS')
        # axs.set_yscale('log')
        # axs.set_ylim(50, 175)
        axs.set_xlim(-0.02, 0.01)
        # axs.set_xlim(left=-0.02)
        # axs.set_ylim(0.05, 1.2)
        # axs.set_yscale("log")
        # axs.set_ylim(4E+21, 3E+23)
        axs.axvline(x=max(psi_to_dsa_func(xcoord_cut)), color='black', lw=3, ls='--',
                    label='pedestal $\Delta n_e$')
        axs.axvline(x=min(psi_to_dsa_func(xcoord_cut)),
                    color='black', lw=3, ls='--')
        if xcoord_type == 'psi':

            axs.set_xlabel('$\psi_N$')
        elif xcoord_type == 'rrsep':

            axs.set_xlabel(r'$(R-R_{sep}) \; [m]$')

        axs.grid(True)
        axs.legend()

    def srmidS_plot(self, xcoord_type):

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == False and withseries == True:

            # series_flag = self.DefaultSettings['series_flag']

            labels = list(self.data['dircomp']['Attempt'].keys())
            colors = pylab.cm.gnuplot2(np.linspace(
                0, 1, len(labels) + 1))
            color_map = dict(zip(labels, colors))

            self.Splot_method(
                iterlist=labels, cl_dic=color_map, xcoord_type=xcoord_type)
