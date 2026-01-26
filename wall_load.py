# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 03:58:00 2026

@author: ychuang
"""


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnchoredText
import matplotlib.pylab as pylab
from midplane_data.cal_sep_length import seplength_calculator


class wall_load_plot:

    def __init__(self, DF, data):

        self.DF = DF
        self.data = data

    def wall_loads_method(self, ax, what, wld, cl):

        # check_pickle_files(where = where)

        # what="wldnek", vmin=None, vmax=None, scale="lin",
        #                     legend=True, labels=None, loc="best", grid=True, save=Fals

        # user to make sure the plotting dimensions are coheret (see SOLPS-ITER manual)

        poly = wld['wall_geometry']
        area = wld['wlarea']
        if what == 'wldnek':
            data = wld[what] + wld['wldnep']
        elif what == 'neudiff':

            data_diff = wld['wldna'][:, :, :51] - wld['wldra']
            data = data_diff.sum(axis=2)

        elif len(wld[what].shape) == 3:

            data = wld[what].sum(axis=2)

        for i in range(poly.shape[1]):
            ax[0].plot([poly[0, i], poly[2, i]], [
                       poly[1, i], poly[3, i]], 'ko-')
            ax[0].annotate(str(i), (poly[0, i], poly[1, i]),
                           textcoords="offset points",
                           xytext=(0, 8),
                           ha='left')

        ax[0].set_xlabel('$R \\; [m]$')
        ax[0].set_ylabel('$Z \\; [m]$')
        ax[0].set_title('Standard surfaces')
        ax[0].set_aspect('equal')

        ax[1].plot(np.linspace(0, poly.shape[1]-1, poly.shape[1]),
                   data[0, :poly.shape[1]] / area[:poly.shape[1]], color=cl)
        ax[1].set_xlabel('element number [-]')
        ax[1].set_yscale("log")
        if what == 'wldnek':
            ax[1].set_ylabel('$[W \\cdot m^{-2}]$')
        ax[1].set_title(what)

        plt.show()

    def plot_wall_loads(self, what):

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == False and withseries == True:

            labels = list(self.data['dircomp']['Attempt'].keys())
            colors = pylab.cm.gnuplot2(np.linspace(
                0, 1, len(labels) + 1))
            color_map = dict(zip(labels, colors))

            fig, ax = plt.subplots(1, 2)

            for aa in labels:

                col = color_map[aa]
                wld = self.data['ft44'][aa]

                self.wall_loads_method(ax=ax, what=what, wld=wld, cl=col)
