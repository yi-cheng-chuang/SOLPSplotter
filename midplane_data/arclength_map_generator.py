# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 19:54:35 2026

@author: ychuang
"""

from midplane_data.SOLPSplotter_PRmap import RP_mapping
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np


class arc_length_measure:

    def __init__(self, DF, data, rp: RP_mapping):

        self.DF = DF
        self.data = data
        self.rp = rp

    def arc_length_map(self, pol_list):

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == True and withseries == False:

            arc_map_dic = {}

            for aa in self.data['dircomp']['Attempt'].keys():

                nx = self.data['b2fgeo'][aa]['nx']
                ny = self.data['b2fgeo'][aa]['ny']

                RadLoc = np.transpose(self.data['grid']['RadLoc'][aa])[
                    1:nx + 1, 1:ny + 1]
                VertLoc = np.transpose(self.data['grid']['VertLoc'][aa])[
                    1:nx + 1, 1:ny + 1]

                pt = len(pol_list)
                arc_map = np.zeros((pt, 36))

                for it, bb in enumerate(pol_list):

                    rline = RadLoc[int(bb), :]
                    zline = VertLoc[int(bb), :]

                    arclength, inp = self.rp.RR_sep_calculator(
                        cr=rline, cz=zline)
                    arc_arcsep = arclength - \
                        np.mean([arclength[18], arclength[19]])

                    arc_map[it, :] = arc_arcsep

                arc_map_dic[aa] = arc_map

            self.data['arc_map'] = arc_map_dic

        else:

            nx = self.data['b2fgeo']['nx']
            ny = self.data['b2fgeo']['ny']

            RadLoc = np.transpose(self.data['grid']['RadLoc'])[
                1:nx + 1, 1:ny + 1]
            VertLoc = np.transpose(self.data['grid']['VertLoc'])[
                1:nx + 1, 1:ny + 1]

            pt = len(pol_list)
            arc_map = np.zeros((pt, 36))

            for it, aa in enumerate(pol_list):

                rline = RadLoc[int(aa), :]
                zline = VertLoc[int(aa), :]

                arclength, inp = self.rp.RR_sep_calculator(cr=rline, cz=zline)
                arc_arcsep = arclength - \
                    np.mean([arclength[18], arclength[19]])

                arc_map[it, :] = arc_arcsep

            self.data['arc_map'] = arc_map
