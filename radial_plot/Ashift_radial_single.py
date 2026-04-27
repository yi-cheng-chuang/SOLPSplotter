# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 17:54:25 2026

@author: ychuang
"""


import numpy as np
import statistics
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib.pylab as pylab
from scipy.optimize import curve_fit
from fit_data.fitting_method import fit_method_collection


class Ashift_midsingle:

    def __init__(self, DF, data, fmc: fit_method_collection):

        self.DF = DF
        self.data = data
        self.fmc = fmc

    def Ashiftplot_method(self, iterlist, cl_dic, xcoord_type):

        # midplane_psi = self.data['midplane_calc']['psi_solps_mid']
        # r_rsep = self.data['midplane_calc']['R_Rsep']
        xcoord_cut = self.data['radial_fit_data']['org']['x_coord_cut']

        # psi_to_dsa_func = interpolate.interp1d(
        #     midplane_psi, r_rsep, fill_value='extrapolate')
        marker = ['o', 's', 'x', '^', 'v', '*', '.', 'o', 's']
        lins_map = dict(zip(iterlist, marker))

        fig, axs = plt.subplots()

        nx = self.data['b2fgeo']['org']['nx']
        ny = self.data['b2fgeo']['org']['ny']
        dat_struc = {'nx': nx, 'ny': ny}

        avg_te_list = []
        fit_exp_dic = {}
        te_range_list = []
        # te_range_list.append(
        #     round(self.data['opacity_poloidal']['4.5']['electron_pedestal_temperature'][0], 2))
        # te_range_list.append(
        #     round(self.data['opacity_poloidal']['4.5']['electron_temperature_separatrix'][0], 2))
        # tie_range_list = []
        # tie_range_list.append(
        #     round(self.data['midplane_profile']['4.5']['mid_ti'][0], 2))
        # tie_range_list.append(
        #     round(self.data['midplane_profile']['4.5']['mid_te'][0], 2))

        for aa in iterlist:

            """
            label= 'core density {} $10^{19}$'.format(aa)

            """

            psi_coord = self.data['midplane_profile'][aa]['psiN']
            r_rsep = self.data['midplane_calc'][aa]['R_Rsep'][1:ny+1]
            psi_to_dsa_func = interpolate.interp1d(
                psi_coord, r_rsep, fill_value='extrapolate')

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

            ne_std = self.data['midplane_profile']['org']['mid_ne']

            ne_norm_inc_dat = return_increase(
                high_dat=mid_ne_pro, std_dat=ne_std)

            mid_neu_pro = self.data['midplane_profile'][aa]['mid_nd']
            insep_nd = self.data['midplane_profile'][aa]['mid_nd'][18]
            outsep_nd = self.data['midplane_profile'][aa]['mid_nd'][19]
            sep_nd = statistics.mean([insep_nd, outsep_nd])

            # max_nd = mid_neu_pro.max()
            nd_std = self.data['midplane_profile']['org']['mid_nd']

            nd_norm_inc_dat = return_increase(
                high_dat=mid_neu_pro, std_dat=nd_std)

            nd_change_percent = np.divide(
                mid_neu_pro, nd_std, out=np.zeros_like(mid_neu_pro, dtype=float), where=nd_std != 0) * 100

            norm_nd = mid_neu_pro / sep_nd

            midS = self.data['midplane_profile'][aa]['mid_S']
            max_S = max(self.data['midplane_profile'][aa]['mid_S'][:35])
            source_std = self.data['midplane_profile']['org']['mid_S']
            norm_S = midS / max_S
            source_change_percent = np.divide(
                midS, source_std, out=np.zeros_like(midS, dtype=float), where=source_std != 0) * 100

            S_norm_inc_dat = return_increase(
                high_dat=midS, std_dat=source_std)

            molsource = self.data['ft44'][aa]['srcml'][59, :, 0]

            exp_fit = False
            rrsep_solps = psi_to_dsa_func(psi_coord)

            if exp_fit:

                index_max = np.where(norm_S == 1)[0][0]
                print("This is index max:")
                print(index_max)

                inside_peak_list = norm_S[:index_max]

                pn = [1.0, 200.5]
                x_sh = rrsep_solps[:index_max]
                popt_an, pcov_an = curve_fit(
                    self.fmc.expfit, x_sh, inside_peak_list, pn)
                # print(popt_an)

                exp_an_fit = self.fmc.expfit(
                    x_sh, popt_an[0], popt_an[1])

                plt.figure()
                plt.yscale('log')
                plt.plot(x_sh, inside_peak_list, 'o-',
                         color='green', label='norm source')
                plt.plot(x_sh, exp_an_fit, color='r', lw=5,
                         label='norm source exponential fit')
                plt.xlabel('rrsep')
                plt.ylabel('Source rate $(m^{-3}/s)$')
                plt.title('Normalize source rate with fits')
                plt.legend()

                fit_exp_dic[aa] = {
                    'exp_an_fit': exp_an_fit, 'popt_an': popt_an}

            te_ped = self.data['opacity_poloidal'][aa]['electron_pedestal_temperature'][0]
            te_sep = self.data['opacity_poloidal'][aa]['electron_temperature_separatrix'][0]

            avg_te = 0.5 * (te_ped + te_sep)
            avg_te_list.append(round(avg_te, 2))

            trim = False
            dat = midS

            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                     'dot7': '2.8', 'one': '3.4'}

            if xcoord_type == 'psi':

                if trim:

                    axs.plot(psi_coord[:33], dat, lw=2,
                             color=cl_dic[aa], label='{}'.format(A_dic[aa]))

                else:

                    axs.plot(psi_coord, dat, lw=2,
                             color=cl_dic[aa], label='{}'.format(A_dic[aa]))

            elif xcoord_type == 'rrsep':

                if trim:

                    axs.plot(rrsep_solps[:33], dat, lw=3, color=cl_dic[aa],
                             label='A = {}'.format(A_dic[aa]))

                else:

                    axs.plot(rrsep_solps, dat, lw=3, color=cl_dic[aa],
                             label='A = {}'.format(A_dic[aa]))

        # axs.set_title('midS')
        axs.set_yscale('log')
        # axs.set_ylim(50, 175)
        axs.set_xlim(-0.02, 0.01)
        # axs.set_ylim(3E+15, 6E+16)
        axs.set_ylim(4E+21, 4E+22)
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

        # print("average te is:")
        # print(avg_te_list)
        # print('te range list:')
        # print(te_range_list)
        # print('tie range list:')
        # print(tie_range_list)

    def Ashiftmid_plot(self, xcoord_type):

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == True and withseries == False:

            # series_flag = self.DefaultSettings['series_flag']
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                     'dot7': '2.8', 'one': '3.4'}
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            labels = list(self.data['dircomp']['Attempt'].keys())
            colors = ['red', 'orange', 'green', 'blue']
            color_map = dict(zip(labels, colors))

            self.Ashiftplot_method(
                iterlist=labels, cl_dic=color_map, xcoord_type=xcoord_type)
