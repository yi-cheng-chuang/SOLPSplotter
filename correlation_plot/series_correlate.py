# -*- coding: utf-8 -*-
"""
Created on Tue Jan  6 19:10:57 2026

@author: ychuang
"""

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np


class plot_series_correlate:

    def __init__(self, DF, data):

        self.DF = DF
        self.data = data

    def plot_series_opq(self, pol_list, dat_type):

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == False and withseries == True:

            result = self.data['opacity_poloidal']
            ll = len(self.data['dircomp']['Attempt'].keys())
            mm = len(pol_list)

            series_list = self.data['dircomp']['Attempt'].keys()

            color_list = ['red', 'salmon', 'orange', 'lime', 'green', 'darkgreen',
                          'cyan', 'deepskyblue', 'navy', 'purple']

            # data_collect_opq = xr.DataArray(np.zeros((ll, mm)),
            #                       coords=[series_list, pol_list],
            #                 dims=['different_density','Poloidal_Location'],
            #                  name = r'dimensionless opaqueness $m$')
            if dat_type == 'all':

                data_collect_opq = np.zeros((mm, ll))
                for ia, s_item in enumerate(self.data['dircomp']['Attempt'].keys()):
                    lb = np.asarray(result[s_item]['dimensionless_opaqueness'])
                    data_collect_opq[:, ia] = lb

            elif dat_type == 'std':

                data_collect_opq = np.zeros((ll, 4))
                for ia, s_item in enumerate(self.data['dircomp']['Attempt'].keys()):
                    lb = np.asarray(result[s_item]['dimensionless_opaqueness'])
                    den_ped = np.asarray(
                        result[s_item]['electron_pedestal_density'])
                    den_sep = np.asarray(
                        result[s_item]['electron_density_separatrix'])
                    std_den_ped = np.asarray(
                        result['1.0']['electron_pedestal_density'])
                    std_den_sep = np.asarray(
                        result['1.0']['electron_density_separatrix'])
                    den_avg = 0.5*(den_ped + den_sep)
                    std_den_avg = 0.5*(std_den_ped + std_den_sep)
                    norm_den_avg = den_avg / std_den_avg
                    # print(norm_den_avg)
                    opq_mean = np.mean(lb)
                    opq_std = np.mean(lb)
                    den_mean = np.mean(norm_den_avg)
                    den_std = np.mean(norm_den_avg)
                    data_collect_opq[ia, 0] = opq_mean
                    data_collect_opq[ia, 1] = opq_std
                    data_collect_opq[ia, 2] = den_mean
                    data_collect_opq[ia, 3] = den_std

            # data_collect_opq =

            self.data['data_collect'] = data_collect_opq

            series = np.array(list(self.data['dircomp']
                              ['Attempt'].keys()), dtype=float)
            shift_list = [round(x, 2) for x in series]

            # ka = 0
            # for k, item in enumerate(self.data['dircomp']['multi_shift']):
            #     shift_list[k] = float(self.data['dircomp']['shift_dic'][item])
            #     # ka = ka + 1

            colors = pylab.cm.gnuplot2(np.linspace(
                0, 1, len(shift_list) + 1))
            color_map = dict(zip(shift_list, colors))

            marker = ['o', 's', 'x', '^', 'v', '*', '.', 'o', 's']
            lins_map = dict(zip(shift_list, marker))

            fig, axs = plt.subplots()

            if dat_type == 'all':

                for p in range(len(shift_list)):
                    x_cor = np.ones(len(pol_list))*shift_list[p]
                    # plt.scatter(x_cor, data_collect_opq[:, p], lins_map[p] + '-', lw=1.5,
                    #             color=color_map[p], label='{}'.format(shift_list[p]))

                    axs.scatter(x_cor, data_collect_opq[:, p], color=color_list[p], label='{}'.format(
                        shift_list[p]))

            elif dat_type == 'std':

                axs.errorbar(
                    shift_list, data_collect_opq[:, 0], yerr=data_collect_opq[:, 1], fmt="o-", label='opaqueness')

                axs.errorbar(
                    shift_list, data_collect_opq[:, 2], yerr=data_collect_opq[:, 3], fmt="o-", label='average density')

            axs.set_xlabel('ratio')
            axs.set_title(
                'Experimental opaqueness verses shift core electron density')
            axs.legend()

        else:
            print('data reorder function is not there yet!')
