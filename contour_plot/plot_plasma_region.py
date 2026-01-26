# -*- coding: utf-8 -*-
"""
Created on Sun Dec  7 15:52:24 2025

@author: ychuang
"""


import numpy as np
import matplotlib.pyplot as plt
from fit_data.fitting_method import fit_method_collection
from contour_plot.contourplot_toolbox import contour_plot_method_collect
from midplane_data.SOLPSplotter_PRmap import RP_mapping


class plasma_region_contour:

    def __init__(self, DF, data):

        self.DF = DF
        self.data = data
        # self.fmc = fmc
        # self.cpmc = cpmc
        # self.rp = rp

    def plasma_region_plot(self):

        # -------------------------------------------------------
        # 1. Load MAST equilibrium (your g-file)
        # -------------------------------------------------------

        g = self.data['gfile']['g']

        # Grid setup from the g-file
        nw = g['nw']
        nh = g['nh']
        rdim = g['rdim']
        zdim = g['zdim']
        rleft = g['rleft']
        zmid = g['zmid']

        R = np.linspace(rleft, rleft + rdim, nw)
        Z = np.linspace(zmid - 0.5 * zdim, zmid + 0.5 * zdim, nh)
        RR, ZZ = np.meshgrid(R, Z)

        # Poloidal flux and normalized flux
        psirz = g['psirz']          # shape (nh, nw)
        simag = g['simag']          # psi at magnetic axis
        sibry = g['sibry']          # psi at separatrix
        psiN = (psirz - simag) / (sibry - simag)  # normalized psi in [0, ~>1]

        # Magnetic axis and (optional) separatrix
        R_axis = g['rmaxis']
        Z_axis = g['zmaxis']

        # Vessel wall from g-file (or replace with your own arrays)
        Rwall = g.get('rlim')   # 1D array of wall R coordinates
        Zwall = g.get('zlim')   # 1D array of wall Z coordinates

        X = self.data['grid']['RadLoc']
        Y = self.data['grid']['VertLoc']

        # -------------------------------------------------------
        # 2. Plot magnetic coordinate (flux-surface) contours
        # -------------------------------------------------------
        fig, ax = plt.subplots(figsize=(5, 7))

        # Filled contours of normalized flux
        levels = np.linspace(psiN.min(), psiN.max(), 30)
        cs = ax.contourf(RR, ZZ, psiN, levels=levels)

        # Contour lines (psi_N flux surfaces)
        ax.contour(RR, ZZ, psiN, levels=levels, colors='k', linewidths=0.5)

        for k in range(X.shape[0]):  # horizontal lines
            if k <= 18:
                ax.plot(X[k, : 24], Y[k, : 24], 'g', linewidth=1)
                ax.plot(X[k, 73:], Y[k, 73:], 'g', linewidth=1)

            if k >= 8 and k <= 18:
                ax.plot(X[k, 25: 73], Y[k, 25: 73], 'cyan', linewidth=1)
            elif k > 18:
                ax.plot(X[k, :], Y[k, :], 'b', linewidth=1)

        # for k in range(X.shape[1]):  # vertical lines

        # Highlight separatrix (psi_N = 1) in red
        ax.contour(RR, ZZ, psiN, levels=[1.0], colors='red', linewidths=2)

        # ax.contour(RR, ZZ, psiN, levels=[0.85], colors='lime', linewidths=2)

        # Mark magnetic axis
        ax.plot(R_axis, Z_axis, 'r.', markersize=8)

        # Plot vessel wall
        if Rwall is not None and Zwall is not None:
            ax.plot(Rwall, Zwall, color='black', linewidth=3)

        ax.plot(X[:, 0], Y[:, 0], 'orange', linewidth=3)
        ax.plot(X[:, -1], Y[:, -1], 'orange', linewidth=3)

        # Colorbar
        cbar = plt.colorbar(cs, ax=ax)
        cbar.set_label(r'$\psi_N$')

        # Labels / formatting
        ax.set_xlabel('R (m)')
        ax.set_ylabel('Z (m)')
        ax.set_aspect('equal')
        ax.set_xlim(R.min(), R.max())
        ax.set_ylim(Z.min(), Z.max())
        # ax.set_title('Magnetic Reconstruction and Flux Surfaces')

        plt.tight_layout()
        plt.show()
