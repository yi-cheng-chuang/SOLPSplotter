# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 19:15:19 2026

@author: ychuang
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnchoredText
import matplotlib.pylab as pylab
from midplane_data.cal_sep_length import seplength_calculator


class series_1para_totargets_plot:

    def __init__(self, DF, data, sc: seplength_calculator):

        self.DF = DF
        self.data = data
        self.sc = sc

    def series_1para_method(self, iterlist, cl_dic, plot_option, ang_list, pol_list, parameter, rad_loc):

        fig, axs = plt.subplots()

        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']

        marker = ['o', 's', 'x', '^', 'v', '*', '.', 'o', 's']
        lins_map = dict(zip(iterlist, marker))

        plt.subplots_adjust(hspace=.0)
        anchored_text_1 = AnchoredText('{}'.format('Electron density [$m^{-3}$]'),
                                       loc='upper center')
        anchored_text_2 = AnchoredText('{}'.format('Electron temperature [eV]'),
                                       loc='lower right')
        anchored_text_3 = AnchoredText('{}'.format('Ion temperature [eV]'),
                                       loc='lower right')

        rad_grid = self.data['grid']['RadLoc']
        vert_grid = self.data['grid']['VertLoc']
        pfluxa_std = self.data['ft44']['1.0']['pfluxa'][:, :, 0]
        fnax_std = self.data['b2wdat']['1.0']['b2npc_fnaxs'][0][1:97, 1:37]

        for aa in iterlist:

            b2fstate = self.data['b2fstate'][aa]
            ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
            Te_J = b2fstate['te'][1:nx+1, 1:ny+1]
            Ti_J = b2fstate['ti'][1:nx+1, 1:ny+1]
            ev = 1.6021766339999999 * pow(10, -19)
            te_pro = Te_J / ev
            ti_pro = Ti_J / ev

            fnax = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
            hz = self.data['b2wdat'][aa]['hz'][1:nx+1, 1:ny+1]
            hy = self.data['b2wdat'][aa]['hy'][1:nx+1, 1:ny+1]
            tor_area = np.multiply(hz, hy)
            fnaxs = np.divide(fnax, tor_area)
            abs_fnaxs = abs(fnaxs)
            fnaxs_std = np.divide(fnax_std, tor_area)

            s_term = self.data['b2wdat'][aa]['b2npc_sna'][0][1:nx+1, 1:ny+1]
            vol = self.data['b2wdat'][aa]['vol'][1:nx+1, 1:ny+1]
            source = np.divide(s_term, vol)

            neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
            molden = self.data['ft44'][aa]['dmb2'][:, :, 0]
            rfluxa = self.data['ft44'][aa]['rfluxa'][:, :, 0]
            rfluxm = self.data['ft44'][aa]['rfluxm'][:, :, 0]
            pfluxa = self.data['ft44'][aa]['pfluxa'][:, :, 0]
            pfluxm = self.data['ft44'][aa]['pfluxm'][:, :, 0]

            # pfluxas = np.divide(pfluxa, tor_area)
            # pfluxms = np.divide(pfluxm, tor_area)

            pfluxas_std = np.multiply(pfluxa_std, tor_area)

            def return_increase(high_dat, std_dat):

                norm_inc_dat = (high_dat - std_dat) * 100 / std_dat

                return norm_inc_dat

            fnaxs_norm_inc_dat = return_increase(
                high_dat=fnaxs, std_dat=fnaxs_std)

            # pfluxa_norm_inc_dat = return_increase(
            #     high_dat=pfluxas, std_dat=pfluxas_std)

            data_dic = {'ne': ne_dat, 'te': te_pro, 'ti': ti_pro, 'polflux': fnaxs, 'source': source,
                        'nd': neuden, 'nm': molden, 'pfluxa': pfluxa, 'abspolflux': abs_fnaxs}

            data = data_dic[parameter]

            if plot_option == 'core':

                dat_list = []

                for ii in pol_list:

                    dat_list.append(data[int(ii), rad_loc])

            elif plot_option == 'inner leg':

                index_list = np.linspace(1, 27, 27)

                dat_list = []

                R_loc = rad_grid[rad_loc, 1:28]

                Z_loc = vert_grid[rad_loc, 1:28]

                arclength, interpolated_points = self.sc.calc_sepx_length(cr=R_loc,
                                                                          cz=Z_loc, plot=False)

                for ii in index_list:

                    dat_list.append(data[int(ii), rad_loc])

            elif plot_option == 'outer leg':

                index_list = np.linspace(70, 95, 26)

                dat_list = []

                R_loc = rad_grid[rad_loc, 70:96][::-1]

                Z_loc = vert_grid[rad_loc, 70:96][::-1]

                arclength, interpolated_points = self.sc.calc_sepx_length(cr=R_loc,
                                                                          cz=Z_loc, plot=False)

                arc_list = index_list[::-1]

                for ii in arc_list:

                    dat_list.append(data[int(ii), rad_loc])

            if plot_option == 'core':
                axs.plot(ang_list, dat_list, lins_map[aa] + '-',
                         lw=1.5, color=cl_dic[aa], label='{}'.format(aa))

                if parameter == 'nd':
                    axs.set_yscale('log')

                axs.set_xlabel('poloidal angle')
            elif plot_option == 'outer leg':
                axs.plot(arclength, dat_list, lins_map[aa] + '-', lw=1.5, color=cl_dic[aa],
                         label='{}'.format(aa))
                # axs.set_ylim(1E+22, 3E+23)

                if parameter == 'nd':
                    axs.set_yscale('log')

                axs.set_xlabel(
                    'separatrix arclength distance to the outer target [m]')
            elif plot_option == 'inner leg':
                axs.plot(arclength, dat_list, lins_map[aa] + '-', lw=1.5, color=cl_dic[aa],
                         label='{}'.format(aa))

                if parameter == 'nd':
                    axs.set_yscale('log')

                axs.set_xlabel(
                    'separatrix arclength distance to the inner target [m]')
                # axs[0].plot(index_list, ne_list, lins_map[aa] + '-', lw=1.5, color=cl_dic[aa],
                #             label='{}'.format(aa))
                # axs[1].plot(index_list, te_list, lins_map[aa] +
                #             '-', lw=1.5, color=cl_dic[aa])

            axs.legend(loc='best')

        # axs[0].add_artist(anchored_text_1)
        # axs[1].add_artist(anchored_text_3)

        if plot_option == 'core':
            axs.axvline(x=0, color='black', lw=3, ls='-', label='LFS')
            axs.axvline(x=180, color='brown', lw=3, ls='-', label='HFS')

        else:
            pass

        # if log_flag:
        #     axs.set_yscale('log')

        axs.legend(loc='best')
        axs.grid(True)
        axs.grid(True)

    def series_1para_totargets_polplot(self, pol_list, parameter, rad_loc):

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == False and withseries == True:

            labels = list(self.data['dircomp']['Attempt'].keys())
            colors = pylab.cm.gnuplot2(np.linspace(
                0, 1, len(labels) + 1))
            color_map = dict(zip(labels, colors))

            ang_list = self.data['angle']['angle_list']

            self.series_1para_method(iterlist=labels, cl_dic=color_map, ang_list=ang_list, rad_loc=rad_loc,
                                     pol_list=pol_list, plot_option='core', parameter=parameter)

            self.series_1para_method(iterlist=labels, cl_dic=color_map, ang_list=ang_list, rad_loc=rad_loc,
                                     pol_list=pol_list, plot_option='inner leg', parameter=parameter)

            self.series_1para_method(iterlist=labels, cl_dic=color_map, ang_list=ang_list, rad_loc=rad_loc,
                                     pol_list=pol_list, plot_option='outer leg', parameter=parameter)
