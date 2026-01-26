# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 10:37:43 2026

@author: ychuang
"""


import matplotlib.pyplot as plt


class poloidal_index_plotter:

    def __init__(self, DF, data):

        self.DF = DF
        self.data = data

    def poloidal_index_plot(self, pol_list):

        RadLoc = self.data['grid']['RadLoc']
        VertLoc = self.data['grid']['VertLoc']

        plt.figure()
        for in_pol in pol_list:
            crloc = RadLoc[:, int(in_pol)]
            czloc = VertLoc[:, int(in_pol)]
            plt.plot(crloc, czloc, color='g',
                     label='R&Zlocation_index{}'.format(in_pol))
        plt.xlabel('R: [m]')
        plt.ylabel('Z: [m]')
        plt.legend(loc='best')
