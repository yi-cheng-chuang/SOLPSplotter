# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 01:12:00 2026

@author: ychuang
"""

from matplotlib.offsetbox import AnchoredText
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker, colors
import matplotlib.tri as tri
from matplotlib.colors import LogNorm
from numpy import ma
from contour_plot.contourplot_toolbox import contour_plot_method_collect
# from midplane_data.SOLPSplotter_PRmap import RP_mapping


class seriesrates_Rectangular_contour:

    def __init__(self, DF, data, cpmc: contour_plot_method_collect): v

    self.DF = DF
    self.data = data
    self.cpmc = cpmc

    def plot_rates(where=".", database=None, group=None, reaction=None,
                   vmin=None, vmax=None, cmap=None, scale="log", mask=True,
                   top=None, bottom=None, right=None, left=None, save=False):

        rates = compute_rates(where=where)

        b2fgmtry = read_b2fgmtry(where=where, verbose=False, save=True)

        if isinstance(where, list) is True:
            raise ValueError('No multiple \"where\" accepted!')

        if not isinstance(vmin, list):
            vmin = [vmin]
        if not isinstance(vmax, list):
            vmax = [vmax]

        if database is None:
            database = []
            group = []
            reaction = []
            for key in rates.keys():
                database.append(rates[key]['database'])
                group.append(rates[key]['group'])
                reaction.append(rates[key]['reaction'])
        else:
            if isinstance(database, list) is False:
                database = [database]
            if isinstance(group, list) is False:
                group = [group]
            if isinstance(reaction, list) is False:
                reaction = [reaction]

        for i in range(len(database)):

            try:
                low = vmin[i]
            except:
                low = vmin[0]
            try:
                up = vmax[i]
            except:
                up = vmax[0]

            for key in rates.keys():
                if (
                        rates[key]['database'] == database[i] and
                        rates[key]['group'] == group[i] and
                        rates[key]['reaction'] == reaction[i]
                ):
                    try:
                        value = rates[key]['reaction_rate']
                        this = ' '.join([database[i], group[i], reaction[i]])
                        break
                    # no data if database == ADAS or others
                    except:
                        pass

            if value != [] and value.max() > 0:

                if database[i] == 'AMMONX':
                    triangles = read_triangle_mesh(
                        where=where, verbose=False, save=False)
                    cells = triangles['cells']
                    nodes = triangles['nodes']
                else:
                    (value, nodes, cells) = to_triangles(
                        nsp=1, vs=value, b2fgmtry=b2fgmtry)

                fancy_title = this + ' ' + \
                    '(' + rates[key]['type'] + ')' '\n' + \
                    rates[key]['equation']

                triplot(where=where, what=fancy_title, value=value, sp=[''], nodes=nodes, cells=cells,
                        vmin=low, vmax=up, cmap=cmap, scale=scale, top=top, bottom=bottom,
                        right=right, left=left, mask=mask, save=save)
            else:
                print()
                print(this + ' missing...')

        plt.show()
