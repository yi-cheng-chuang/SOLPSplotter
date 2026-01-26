# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 15:36:16 2026

@author: ychuang
"""


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnchoredText
import matplotlib.pylab as pylab
from midplane_data.cal_sep_length import seplength_calculator


class series_ndnmtotarget_poloidal_plot:

    def __init__(self, DF, data, sc: seplength_calculator):

        self.DF = DF
        self.data = data
        self.sc = sc

    def ndnm_poloidal_method(self, iterlist, cl_dic, plot_option, ang_list, pol_list, log_flag, rad_loc):

        fig, axs = plt.subplots(2, 1)

        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']

        marker = ['o', 's', 'x', '^', 'v', '*', '.', 'o', 's']
        lins_map = dict(zip(iterlist, marker))

        plt.subplots_adjust(hspace=.0)
        anchored_text_1 = AnchoredText('{}'.format('Atomic neutral density [$m^{-3}$]'),
                                       loc='upper center')
        anchored_text_2 = AnchoredText('{}'.format('Electron temperature [eV]'),
                                       loc='lower right')
        anchored_text_3 = AnchoredText('{}'.format('Molecular neutral density [$m^{-3}$]'),
                                       loc='lower right')

        rad_grid = self.data['grid']['RadLoc']
        vert_grid = self.data['grid']['VertLoc']

        for aa in iterlist:

            s_term = self.data['b2wdat'][aa]['b2npc_sna'][0][1:nx+1, 1:ny+1]
            vol = self.data['b2wdat'][aa]['vol'][1:nx+1, 1:ny+1]
            source = np.divide(s_term, vol)

            neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
            molden = self.data['ft44'][aa]['dmb2'][:, :, 0]

            if plot_option == 'core':

                nd_list = []
                s_list = []
                md_list = []

                for ii in pol_list:

                    nd_list.append(neuden[int(ii), rad_loc])
                    s_list.append(source[int(ii), rad_loc])
                    md_list.append(molden[int(ii), rad_loc])

            elif plot_option == 'inner leg':

                index_list = np.linspace(1, 27, 27)

                nd_list = []
                s_list = []
                md_list = []

                R_loc = rad_grid[rad_loc, 1:28]

                Z_loc = vert_grid[rad_loc, 1:28]

                arclength, interpolated_points = self.sc.calc_sepx_length(cr=R_loc,
                                                                          cz=Z_loc, plot=False)

                for ii in index_list:

                    nd_list.append(neuden[int(ii), rad_loc])
                    s_list.append(source[int(ii), rad_loc])
                    md_list.append(molden[int(ii), rad_loc])

            elif plot_option == 'outer leg':

                index_list = np.linspace(70, 95, 26)

                nd_list = []
                s_list = []
                md_list = []

                R_loc = rad_grid[rad_loc, 70:96][::-1]

                Z_loc = vert_grid[rad_loc, 70:96][::-1]

                arclength, interpolated_points = self.sc.calc_sepx_length(cr=R_loc,
                                                                          cz=Z_loc, plot=False)

                arc_list = index_list[::-1]

                for ii in arc_list:

                    nd_list.append(neuden[int(ii), rad_loc])
                    s_list.append(source[int(ii), rad_loc])
                    md_list.append(molden[int(ii), rad_loc])

            if plot_option == 'core':
                axs[0].plot(ang_list, nd_list, lins_map[aa] + '-',
                            lw=1.5, color=cl_dic[aa], label='{}'.format(aa))
                axs[1].plot(ang_list, md_list, lins_map[aa] +
                            '-', lw=1.5, color=cl_dic[aa])
                axs[1].set_xlabel('poloidal angle')
            elif plot_option == 'outer leg':
                axs[0].plot(arclength, nd_list, lins_map[aa] + '-', lw=1.5, color=cl_dic[aa],
                            label='{}'.format(aa))
                axs[1].plot(arclength, md_list, lins_map[aa] +
                            '-', lw=1.5, color=cl_dic[aa])
                axs[0].set_yscale('log')
                axs[1].set_yscale('log')
                axs[1].set_xlabel(
                    'separatrix arclength distance to the outer target [m]')
            elif plot_option == 'inner leg':
                axs[0].plot(arclength, nd_list, lins_map[aa] + '-', lw=1.5, color=cl_dic[aa],
                            label='{}'.format(aa))
                axs[1].plot(arclength, md_list, lins_map[aa] +
                            '-', lw=1.5, color=cl_dic[aa])
                axs[0].set_yscale('log')
                axs[1].set_yscale('log')
                axs[1].set_xlabel(
                    'separatrix arclength distance to the inner target [m]')
                # axs[0].plot(index_list, ne_list, lins_map[aa] + '-', lw=1.5, color=cl_dic[aa],
                #             label='{}'.format(aa))
                # axs[1].plot(index_list, te_list, lins_map[aa] +
                #             '-', lw=1.5, color=cl_dic[aa])

            axs[0].legend(loc='best')

        axs[0].add_artist(anchored_text_1)
        axs[1].add_artist(anchored_text_3)

        if plot_option == 'core':
            axs[0].axvline(x=0, color='black', lw=3, ls='-')
            axs[1].axvline(x=0, color='black', lw=3, ls='-', label='LFS')
            axs[0].axvline(x=180, color='brown', lw=3, ls='-')
            axs[1].axvline(x=180, color='brown', lw=3, ls='-', label='HFS')

        else:
            pass

        if log_flag:
            axs[0].set_yscale('log')

        axs[1].legend(loc='best')
        axs[0].grid(True)
        axs[1].grid(True)

    def serise_ndnm_polplot(self, pol_list, log_flag, rad_loc):

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == False and withseries == True:

            labels = list(self.data['dircomp']['Attempt'].keys())
            colors = pylab.cm.gnuplot2(np.linspace(
                0, 1, len(labels) + 1))
            color_map = dict(zip(labels, colors))

            ang_list = self.data['angle']['angle_list']

            self.ndnm_poloidal_method(iterlist=labels, cl_dic=color_map, ang_list=ang_list, rad_loc=rad_loc,
                                      pol_list=pol_list, plot_option='core', log_flag=log_flag)

            self.ndnm_poloidal_method(iterlist=labels, cl_dic=color_map, ang_list=ang_list, rad_loc=rad_loc,
                                      pol_list=pol_list, plot_option='inner leg', log_flag=log_flag)

            self.ndnm_poloidal_method(iterlist=labels, cl_dic=color_map, ang_list=ang_list, rad_loc=rad_loc,
                                      pol_list=pol_list, plot_option='outer leg', log_flag=log_flag)


"""

b2fstate = self.data['b2fstate'][aa]
ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
Te_J = b2fstate['te'][1:nx+1, 1:ny+1]
Ti_J = b2fstate['ti'][1:nx+1, 1:ny+1]
ev = 1.6021766339999999 * pow(10, -19)
te_pro = Te_J / ev
ti_pro = Ti_J / ev


ne_list = []
te_list = []
ti_list = []

ne_list.append(ne_dat[int(ii), rad_loc])
te_list.append(te_pro[int(ii), rad_loc])
ti_list.append(ti_pro[int(ii), rad_loc])

"""
