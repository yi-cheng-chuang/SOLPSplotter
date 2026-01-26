# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 21:31:15 2025

@author: ychuang
"""

from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np


class srtarget_radial:

    def __init__(self, DF, data):

        self.DF = DF
        self.data = data

    def srtargetNT_method(self, log_flag, cl_dic, iterlist, side):

        fig, axs = plt.subplots(1, 3)

        plt.subplots_adjust(hspace=.0)
        if side == 'inner target':

            anchored_text_1 = AnchoredText('{}'.format(
                'inner target electron density [$m^{-3}$]'), loc='upper left')
            anchored_text_2 = AnchoredText('{}'.format(
                'inner target electron temperature [eV]'), loc='upper left')
            anchored_text_3 = AnchoredText('{}'.format(
                'inner target source [$m^{-3}*s^{-1}$]'), loc='upper left')

        elif side == 'outer target':

            anchored_text_1 = AnchoredText('{}'.format(
                'outer target electron density [$m^{-3}$]'), loc='upper left')
            anchored_text_2 = AnchoredText('{}'.format(
                'outer target electron temperature [eV]'), loc='upper left')
            anchored_text_3 = AnchoredText('{}'.format(
                'outer target source [$m^{-3}*s^{-1}$]'), loc='upper left')

        # print('this is 201:')

        for aa in iterlist:

            psi_dic = self.data['target_profile'][aa]['psiN']
            ne_dic = self.data['target_profile'][aa]['ne']
            te_dic = self.data['target_profile'][aa]['te']
            sx_dic = self.data['target_profile'][aa]['source']
            neuden_dic = self.data['target_profile'][aa]['neuden']

            psi_list = psi_dic[side]
            ne_list = ne_dic[side]
            te_list = te_dic[side]
            sx_list = sx_dic[side]

            axs[0].plot(psi_list, ne_list, '-',
                        color=cl_dic[aa], label='{}'.format(aa))
            axs[1].plot(psi_list, te_list, '-',
                        color=cl_dic[aa], label='{}'.format(aa))
            axs[2].plot(psi_list, sx_list, '-',
                        color=cl_dic[aa], label='{}'.format(aa))
            axs[0].set_xlabel('$\psi_N$')
            axs[1].set_xlabel('$\psi_N$')
            axs[2].set_xlabel('$\psi_N$')
            axs[0].axvline(x=1, color='black', lw=3, ls='--')
            axs[1].axvline(x=1, color='black', lw=3, ls='--')
            axs[2].axvline(x=1, color='black', lw=3, ls='--')
            axs[0].add_artist(anchored_text_1)
            axs[1].add_artist(anchored_text_2)
            axs[2].add_artist(anchored_text_3)
            axs[0].set_title('target plot')

            axs[0].legend(loc='lower right')

    def srtargetNT_combine_method(self, log_flag, cl_dic, iterlist):

        fig, axs = plt.subplots(1, 2)

        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        dat_struc = {'nx': nx, 'ny': ny}

        plt.subplots_adjust(hspace=.0)

        anchored_text_1 = AnchoredText('{}'.format(
            'Electron density at target [$m^{-3}$]'), loc='upper left')
        anchored_text_2 = AnchoredText('{}'.format(
            'Electron temperature at target [eV]'), loc='upper left')
        # anchored_text_3 = AnchoredText('{}'.format('Source at target [$m^{-3}*s^{-1}$]'), loc='upper left')

        for aa in iterlist:

            psi_dic = self.data['target_profile'][aa]['psiN']
            ne_dic = self.data['target_profile'][aa]['ne']
            te_dic = self.data['target_profile'][aa]['te']

            psi_list_in = psi_dic['inner target']
            ne_list_in = ne_dic['inner target']
            te_list_in = te_dic['inner target']
            # sx_list_in = sx_dic['inner target']

            psi_list_out = psi_dic['outer target']
            ne_list_out = ne_dic['outer target']
            te_list_out = te_dic['outer target']
            # sx_list_out = sx_dic['outer target']

            axs[0].plot(psi_list_in, ne_list_in, '-', color=cl_dic[aa],
                        label='HFS {}'.format(aa))
            axs[1].plot(psi_list_in, te_list_in, '-', color=cl_dic[aa],
                        label='HFS {}'.format(aa))
            # axs[2].plot(psi_list_in, sx_list_in,'-', color = cl_dic[ad], label= 'inner {:.3E} (1/s)'.format(label_ad))
            axs[0].plot(psi_list_out, ne_list_out, '--', color=cl_dic[aa],
                        label='LFS {}'.format(aa))
            axs[1].plot(psi_list_out, te_list_out, '--', color=cl_dic[aa],
                        label='LFS {}'.format(aa))
            # axs[2].plot(psi_list_out, sx_list_out,'--', color = cl_dic[ad], label= 'outer {:.3E} (1/s)'.format(label_ad))
            axs[0].set_xlabel('$\psi_N$')
            axs[1].set_xlabel('$\psi_N$')
            # axs[2].set_xlabel('$\psi_N$')
            axs[0].axvline(x=1, color='black', lw=3, ls='--')
            axs[1].axvline(x=1, color='black', lw=3, ls='--')
            # axs[2].axvline(x= 1, color='black', lw=3, ls='--')
            axs[0].add_artist(anchored_text_1)
            axs[1].add_artist(anchored_text_2)
            # axs[2].add_artist(anchored_text_3)
            axs[0].set_title('target NT combine')

            axs[0].legend(loc='lower right')

    def srtargetnd_combine_method(self, log_flag, cl_dic, iterlist):

        fig, axs = plt.subplots(1, 2)

        plt.subplots_adjust(hspace=.0)

        anchored_text_1 = AnchoredText('{}'.format(
            'Source at target [$m^{-3}*s^{-1}$]'), loc='upper left')
        anchored_text_2 = AnchoredText('{}'.format(
            'Atomic neutral density at target [$m^{-3}$]'), loc='upper left')

        for aa in iterlist:

            psi_dic = self.data['target_profile'][aa]['psiN']
            sx_dic = self.data['target_profile'][aa]['source']
            neuden_dic = self.data['target_profile'][aa]['neuden']

            psi_list_in = psi_dic['inner target']

            sx_list_in = sx_dic['inner target']
            nd_list_in = neuden_dic['inner target']

            psi_list_out = psi_dic['outer target']

            sx_list_out = sx_dic['outer target']
            nd_list_out = neuden_dic['outer target']

            axs[0].plot(psi_list_in, sx_list_in, '-', color=cl_dic[aa],
                        label='HFS {}'.format(aa))
            axs[1].plot(psi_list_in, nd_list_in, '-', color=cl_dic[aa],
                        label='HFS {}'.format(aa))
            axs[0].plot(psi_list_out, sx_list_out, '--', color=cl_dic[aa],
                        label='LFS {}'.format(aa))
            axs[1].plot(psi_list_out, nd_list_out, '--', color=cl_dic[aa],
                        label='LFS {}'.format(aa))
            axs[0].set_xlabel('$\psi_N$')
            axs[1].set_xlabel('$\psi_N$')
            axs[0].axvline(x=1, color='black', lw=3, ls='--')
            axs[1].axvline(x=1, color='black', lw=3, ls='--')
            axs[0].add_artist(anchored_text_1)
            axs[1].add_artist(anchored_text_2)
            axs[0].set_title(
                'target nd combine')

            axs[0].legend(loc='lower right')

    def series_targetNT(self, log_flag, match):

        if self.DF.withshift == False and self.DF.withseries == True:

            labels = list(self.data['dircomp']['Attempt'].keys())
            colors = pylab.cm.gnuplot2(np.linspace(
                0, 1, len(labels) + 1))
            color_map = dict(zip(labels, colors))

            side_list = ['inner target', 'outer target']

            if match:

                self.srtargetNT_combine_method(
                    iterlist=labels, cl_dic=color_map, log_flag=log_flag)

                self.srtargetnd_combine_method(
                    iterlist=labels, cl_dic=color_map, log_flag=log_flag)

            else:

                for sk in side_list:

                    self.srtargetNT_method(
                        iterlist=labels, cl_dic=color_map, log_flag=log_flag, side=sk)
