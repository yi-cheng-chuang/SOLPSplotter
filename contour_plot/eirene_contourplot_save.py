# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 18:54:47 2025

@author: ychuang
"""

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from numpy import ma
from matplotlib.offsetbox import AnchoredText
from matplotlib import colors, cm
from twscan_module.twinscan_prepare import twscan_assist


class series_Eirene_contour:

    def __init__(self, DF, data):

        self.DF = DF
        self.data = data

    def eirene_contourplot_tool(self, simudir, data, plot33, vmin=None, vmax=None):

        data_mask = ma.masked_where(data <= 0, data)

        datamap = np.abs(data_mask)
        CMAP = cm.viridis

        NORM_data = colors.LogNorm(1E+11, np.nanmax(datamap))
        # NORM_data = colors.LogNorm(np.nanmin(datamap), np.nanmax(datamap))
        # Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())

        # Alternatively use fort.33
        Nodes = np.fromfile('{}/fort.33'.format(simudir), sep=' ')
        NN = int(Nodes[0])
        XNodes = Nodes[1:NN+1]
        YNodes = Nodes[NN+1:]

        if plot33:

            numberlist = np.zeros(NN)
            for i in range(NN):
                numberlist[i] = i

            plt.figure()
            plt.scatter(XNodes[:500], YNodes[:500])

        Triangles = np.loadtxt('{}/fort.34'.format(simudir),
                               skiprows=1, usecols=(1, 2, 3))  # Alternatively use fort.34
        # print(Triangles -1)

        TP = tri.Triangulation(XNodes, YNodes, triangles=(Triangles - 1))

        # CMAP = cm.viridis

        fig, axs = plt.subplots()

        axs.tripcolor(TP, data, shading='flat', cmap='jet', norm=NORM_data)
        axs.set_title('Neutral density contour plot')
        axs.set_aspect('equal')

        axs.set_ylim(-200, 0)
        # plt.title('Neutral density contour plot outerleg')

        SM_data = cm.ScalarMappable(NORM_data, 'jet')
        cbar = plt.colorbar(SM_data)

        # Label
        cbar.set_label(r"[m$^{-3}$]")

    def eirene_percent_tripcolor_tool(self, simudir, percent_data, lognorm, plot33=False,
                                      vmin=None, vmax=None, nlevels=None,
                                      cmap='jet', center_zero=True,
                                      title="EIRENE value change (%)",
                                      ylim=(-200, 0), xlim=None,
                                      cbar_label="[%]",
                                      shading="flat"):
        """
        Linear (no LogNorm) tripcolor plot for EIRENE percentage-change data.

        Notes
        -----
        - tripcolor does NOT support `levels=` like tricontourf.
        - If you pass nlevels (int), this function will emulate discrete levels via BoundaryNorm.
          If nlevels is None, you get continuous colors (linear scaling).
        """

        # --- Read nodes (fort.33) ---
        Nodes = np.fromfile(f"{simudir}/fort.33", sep=" ")
        NN = int(Nodes[0])
        XNodes = Nodes[1:NN+1]*1E-2
        YNodes = Nodes[NN+1:]*1E-2

        # percent_data = np.asarray(percent_data)
        # if percent_data.shape[0] != NN:
        #     raise ValueError(
        #         f"percent_data length ({percent_data.shape[0]}) != NN from fort.33 ({NN}).")

        # Optional geometry preview
        if plot33:
            plt.figure()
            plt.scatter(XNodes[:min(500, NN)], YNodes[:min(500, NN)], s=8)
            plt.gca().set_aspect("equal")
            plt.title("fort.33 node preview")
            plt.show()

        # --- Read triangles (fort.34) and build triangulation ---
        Triangles = np.loadtxt(f"{simudir}/fort.34",
                               skiprows=1, usecols=(1, 2, 3))
        TP = tri.Triangulation(
            XNodes, YNodes, triangles=(Triangles - 1).astype(int))

        # Mask invalid values
        data = ma.masked_invalid(percent_data)
        finite = np.isfinite(percent_data)
        if not np.any(finite):
            raise ValueError("percent_data contains no finite values.")

        dmin = float(np.nanmin(percent_data[finite]))
        dmax = float(np.nanmax(percent_data[finite]))

        # Choose vmin/vmax (linear)
        if vmin is None or vmax is None:
            if center_zero:
                m = max(abs(dmin), abs(dmax))
                vmin = -m if vmin is None else vmin
                vmax = m if vmax is None else vmax
            else:
                vmin = dmin if vmin is None else vmin
                vmax = dmax if vmax is None else vmax

        # Optional: emulate "levels" using BoundaryNorm (DISCRETE bands)
        if lognorm:
            data_mask = ma.masked_where(data <= 0, data)
            datamap = np.abs(data_mask)
            norm = colors.LogNorm(vmin, vmax)

            # norm = None
            # if nlevels is not None:
            #     nlevels = int(nlevels)
            #     if nlevels < 2:
            #         raise ValueError("nlevels must be >= 2 if provided.")

            #     logmin = np.log10(vmin)   # = 21
            #     logmax = np.log10(vmax)   # ≈ 23.69897

            #     boundaries = np.logspace(logmin, logmax, nlevels)
            #     norm = colors.BoundaryNorm(boundaries, ncolors=256, clip=True)

        else:

            norm = None
            if nlevels is not None:
                nlevels = int(nlevels)
                if nlevels < 2:
                    raise ValueError("nlevels must be >= 2 if provided.")
                boundaries = np.linspace(vmin, vmax, nlevels)
                norm = colors.BoundaryNorm(boundaries, ncolors=256, clip=True)

        # --- Plot with tripcolor (NO LogNorm) ---
        fig, ax = plt.subplots()
        pc = ax.tripcolor(TP, data, cmap=cmap, norm=norm, shading=shading)

        ax.set_title(title)
        ax.set_aspect("equal")
        ax.set_ylim(*ylim)
        if xlim is not None:
            ax.set_xlim(*xlim)

        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label(cbar_label)

    def eirene_contour_plot(self, plot_type, eirene_param):

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == False and withseries == False:

            simu_dir = self.data['dirdata']['simudir']
            dat = self.data['ft46']['pdena'][:, 0]
            abs_dat = np.absolute(dat)

            self.eirene_contourplot_method(simudir=simu_dir, data=abs_dat)

        elif withshift == False and withseries == True:

            if plot_type == 'contour':

                for aa in self.data['dircomp']['Attempt'].keys():

                    simu_dir = self.data['dirdata']['simudir'][aa]
                    dat = self.data['ft46'][aa][eirene_param][:, 0]
                    # abs_dat = np.absolute(dat)
                    eirene_dic = self.Eirene_name_dic()

                    self.eirene_percent_tripcolor_tool(simudir=simu_dir, percent_data=dat,
                                                       lognorm=True,
                                                       plot33=False,
                                                       vmin=1E+11, vmax=2E+19, nlevels=30,
                                                       cmap='jet', center_zero=True,
                                                       title="{}".format(
                                                           eirene_dic[eirene_param]),
                                                       ylim=(-2, 0.5), xlim=None,
                                                       cbar_label="[$m^{-3}$]",
                                                       shading="flat")

            elif plot_type == 'change':

                simu_dir = self.data['dirdata']['simudir']['1.0']
                std_dat = self.data['ft46']['1.0'][eirene_param][:, 0]
                high_dat = self.data['ft46']['4.5'][eirene_param][:, 0]
                norm_inc_dat = (high_dat - std_dat) * 100 / std_dat

                eirene_dic = self.Eirene_name_dic()

                self.eirene_percent_tripcolor_tool(simudir=simu_dir, percent_data=norm_inc_dat, plot33=False,
                                                   vmin=-150, vmax=350, nlevels=30, lognorm=False,
                                                   cmap='jet', center_zero=True,
                                                   title="{} change (%)".format(
                                                       eirene_dic[eirene_param]),
                                                   ylim=(-2, 0.5), xlim=None,
                                                   cbar_label="[%]",
                                                   shading="flat")

        elif withshift == True and withseries == False:

            if plot_type == 'contour':

                for aa in self.data['dircomp']['Attempt'].keys():

                    simu_dir = self.data['dirdata']['simudir'][aa]
                    dat = self.data['ft46'][aa][eirene_param][:, 0]
                    # abs_dat = np.absolute(dat)
                    eirene_dic = self.Eirene_name_dic()

                    self.eirene_percent_tripcolor_tool(simudir=simu_dir, percent_data=dat,
                                                       lognorm=True,
                                                       plot33=False,
                                                       vmin=1E+11, vmax=2E+19, nlevels=30,
                                                       cmap='jet', center_zero=True,
                                                       title="{}".format(
                                                           eirene_dic[eirene_param]),
                                                       ylim=(-2, 0.5), xlim=None,
                                                       cbar_label="[$m^{-3}$]",
                                                       shading="flat")

        else:
            print('eirene_contour_plot function needs to be checked!')

    def Eirene_name_dic(self):

        EireneDict = {'pdena': 'Atomic neutral density $(m^{-3})$',
                      'pdenm': 'Molecular D2 Density $(m^{-3})$',
                      'edena': 'Neutral D Energy Density $(eV*m^{-3})$',
                      'edenm': 'Molecular D2 Energy Density $(eV*m^{-3})$'}

        return EireneDict
