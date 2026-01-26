# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 01:48:56 2026

@author: ychuang
"""


import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np


class series_wall_loads_plot:

    def __init__(self, DF, data):

        self.DF = DF
        self.data = data

    def plot_wall_loads(where=".", what="wldnek", vmin=None, vmax=None, scale="lin",
                        legend=True, labels=None, loc="best", grid=True, save=False):

        # check_pickle_files(where = where)

        (neut, wld) = read_ft44(where=where, save=True)

        # user to make sure the plotting dimensions are coheret (see SOLPS-ITER manual)

        poly = wld['poly']
        area = wld['wlarea']
        if what == 'wldnek':
            data = wld[what] + wld['wldnep']
        elif len(wld[what].shape) == 3:
            data = wld[what].sum(axis=2)

        fig, ax = plt.subplots(1, 2)

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
                   data[:poly.shape[1], 0] / area[:poly.shape[1]], 'r-')
        ax[1].set_xlabel('element number [-]')
        if what == 'wldnek':
            ax[1].set_ylabel('$[W \\cdot m^{-2}]$')
        ax[1].set_title(what)

        plt.show()

        return
