# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 05:34:20 2025

@author: ychuang
"""


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnchoredText
import matplotlib.pylab as pylab


class series_polfluxndS_polplot:

    def __init__(self, DF, data):

        self.DF = DF
        self.data = data

    def srpolfluxndS_method(self, iterlist, cl_dic, plot_option, ang_list, pol_list, log_flag, rad_loc):

        fig, axs = plt.subplots(3, 1)

        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']

        plt.subplots_adjust(hspace=.0)
        anchored_text_1 = AnchoredText('{}'.format('Poloidal flux [$m^{-2}*s^{-1}$]'),
                                       loc='lower center')
        anchored_text_2 = AnchoredText('{}'.format('Source [$m^{-3}*s^{-1}$]'),
                                       loc='upper center')
        anchored_text_3 = AnchoredText('{}'.format('Neutral density [$m^{-3}$]'),
                                       loc='upper right')

        for aa in iterlist:

            fnax = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
            hz = self.data['b2wdat'][aa]['hz'][1:nx+1, 1:ny+1]
            hy = self.data['b2wdat'][aa]['hy'][1:nx+1, 1:ny+1]
            tor_area = np.multiply(hz, hy)
            fnaxs = np.divide(fnax, tor_area)

            s_term = self.data['b2wdat'][aa]['b2npc_sna'][0][1:nx+1, 1:ny+1]
            vol = self.data['b2wdat'][aa]['vol'][1:nx+1, 1:ny+1]
            source = np.divide(s_term, vol)

            neuden = self.data['ft44'][aa]['dab2'][:, :, 0]

            # st = int(pol_list[0])
            # ed = int(pol_list[-1]) + 1

            if plot_option == 'core':

                fx_list = []
                s_list = []
                nd_list = []

                for ii in pol_list:

                    fx_list.append(fnaxs[int(ii), rad_loc])
                    s_list.append(source[int(ii), rad_loc])
                    nd_list.append(neuden[int(ii), rad_loc])

            elif plot_option == 'inner leg':

                index_list = np.linspace(1, 27, 27)

                fx_list = []
                s_list = []
                nd_list = []

                for ii in index_list:

                    fx_list.append(abs(fnaxs[int(ii), rad_loc]))
                    s_list.append(source[int(ii), rad_loc])
                    nd_list.append(neuden[int(ii), rad_loc])

            elif plot_option == 'outer leg':

                index_list = np.linspace(70, 95, 26)

                fx_list = []
                s_list = []
                nd_list = []

                for ii in index_list:

                    fx_list.append(fnaxs[int(ii), rad_loc])
                    s_list.append(source[int(ii), rad_loc])
                    nd_list.append(neuden[int(ii), rad_loc])

            if log_flag:
                axs[1].set_yscale('log')
                axs[2].set_yscale('log')
            else:
                pass

            axs[0].set_title('fluxndS')

            if plot_option == 'core':
                axs[0].plot(ang_list, fx_list, color=cl_dic[aa],
                            label='{} (1/s)'.format(aa))
                axs[1].plot(ang_list, s_list, color=cl_dic[aa])
                axs[2].plot(ang_list, nd_list, color=cl_dic[aa])
            else:
                axs[0].plot(index_list, fx_list, color=cl_dic[aa],
                            label='{}'.format(aa))
                axs[1].plot(index_list, s_list, color=cl_dic[aa])
                axs[2].plot(index_list, nd_list, color=cl_dic[aa])

            axs[2].set_xlabel('poloidal angle')
            axs[0].legend(loc='best')

            axs[0].add_artist(anchored_text_1)
            axs[1].add_artist(anchored_text_2)
            axs[2].add_artist(anchored_text_3)

            if plot_option == 'core':
                axs[0].axvline(x=0, color='black', lw=3, ls='-')
                axs[1].axvline(x=0, color='black', lw=3, ls='-', label='LFS')
                axs[2].axvline(x=0, color='black', lw=3, ls='-')
                axs[0].axvline(x=180, color='brown', lw=3, ls='-')
                axs[1].axvline(x=180, color='brown', lw=3, ls='-', label='HFS')
                axs[2].axvline(x=180, color='brown', lw=3, ls='-')

            else:
                pass

            axs[1].legend(loc='best')

    def srpolfluxndS_plot(self, pol_list, log_flag, rad_loc):

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == False and withseries == True:

            # series_flag = self.DefaultSettings['series_flag']

            labels = list(self.data['dircomp']['Attempt'].keys())
            colors = pylab.cm.gnuplot2(np.linspace(
                0, 1, len(labels) + 1))
            color_map = dict(zip(labels, colors))

            ang_list = self.data['angle']['angle_list']

            self.srpolfluxndS_method(iterlist=labels, cl_dic=color_map, ang_list=ang_list, rad_loc=rad_loc,
                                     pol_list=pol_list, plot_option='core', log_flag=log_flag)

            self.srpolfluxndS_method(iterlist=labels, cl_dic=color_map, ang_list=ang_list, rad_loc=rad_loc,
                                     pol_list=pol_list, plot_option='inner leg', log_flag=log_flag)

            self.srpolfluxndS_method(iterlist=labels, cl_dic=color_map, ang_list=ang_list, rad_loc=rad_loc,
                                     pol_list=pol_list, plot_option='outer leg', log_flag=log_flag)

        else:
            print('neteTS_plot, please check the series flag')
