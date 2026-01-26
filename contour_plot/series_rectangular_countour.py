# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 23:22:34 2026

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


class series_Rectangular_contour:

    def __init__(self, DF, data, cpmc: contour_plot_method_collect):

        self.DF = DF
        self.data = data
        self.cpmc = cpmc

    def twrectangular_test(self):

        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']

        psi_map = np.transpose(self.data['psi']['psival'])[1:nx+1, 1:ny+1]

        Ra = np.arange(1, nx+1)
        extend_Ra = np.tile(Ra, (ny, 1))
        Ra_map = np.transpose(extend_Ra)

        self.data['save_rectangular'] = Ra_map

        # pol_list = self.DF.poloidal_loc_list

        # print(pol_list)

    def series_rectangular_contourplot(self, plot_name, norm_type, label_type, pol_loc_list, sted_index_list):

        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']

        psi_map = np.transpose(self.data['psi']['psival'])[1:nx+1, 1:ny+1]

        ang_label = self.data['angle']['angle_list']
        extend_ang = np.tile(ang_label, (ny, 1))
        ang_map = np.transpose(extend_ang)

        ins = sted_index_list[0]
        out = sted_index_list[1]

        Ra = np.arange(1, nx+1)
        extend_Ra = np.tile(Ra, (ny, 1))
        Ra_map = np.transpose(extend_Ra)

        arc_map = self.data['arc_map']

        rec_dic = {'angle': ang_map, 'index': Ra_map}

        self.data['save_rectangular'] = rec_dic

        RadLoc = np.transpose(self.data['grid']['RadLoc'])[1:nx + 1, 1:ny + 1]
        VertLoc = np.transpose(self.data['grid']['VertLoc'])[
            1:nx + 1, 1:ny + 1]

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == False and withseries == True:

            for aa in self.data['dircomp']['Attempt'].keys():

                fig, axs = plt.subplots()

                b2fstate = self.data['b2fstate'][aa]
                ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
                Te_J = b2fstate['te'][1:nx+1, 1:ny+1]

                ev = 1.6021766339999999 * pow(10, -19)
                te_dat = Te_J / ev

                source = self.data['b2wdat'][aa]['b2npc_sna'][0][1:nx+1, 1:ny+1]
                vol = self.data['b2wdat'][aa]['vol'][1:nx+1, 1:ny+1]
                sx = np.divide(source, vol)

                neuden_dat = self.data['ft44'][aa]['dab2'][:, :, 0]
                denm_dat = self.data['ft44'][aa]['dmb2'][:, :, 0]
                molsource_dat = self.data['ft44'][aa]['srcml'][:, :, 0]

                hz = self.data['b2wdat'][aa]['hz'][1:nx+1, 1:ny+1]
                hx = self.data['b2wdat'][aa]['hx'][1:nx+1, 1:ny+1]
                pol_area = np.multiply(hz, hx)
                rfluxm = self.data['ft44'][aa]['rfluxm'][:, :, 0]
                rfluxms = np.multiply(rfluxm, pol_area)
                rfluxm_dat = np.abs(rfluxms)

                rfluxa = self.data['ft44'][aa]['rfluxa'][:, :, 0]
                rfluxas = np.multiply(rfluxa, pol_area)
                rfluxa_dat = np.abs(rfluxas)

                nd_std_dat = self.data['ft44']['1.0']['dab2'][:, :, 0]
                nm_std_dat = self.data['ft44']['1.0']['dmb2'][:, :, 0]
                source_std = self.data['b2wdat']['1.0']['b2npc_sna'][0][1:nx+1, 1:ny+1]
                vol_std = self.data['b2wdat']['1.0']['vol'][1:nx+1, 1:ny+1]
                sx_std = np.divide(source_std, vol_std)

                def return_increase(high_dat, std_dat):

                    norm_inc_dat = (high_dat - std_dat) * 100 / std_dat

                    return norm_inc_dat

                if plot_name == 'sx':

                    datname = 'Source'
                    title_name = 'ion source rate [$m^{-3} s^{-1}$]'
                    dat = sx
                    vmin_0 = 1E+20
                    vmax_0 = 2E+23

                elif plot_name == 'neuden':

                    datname = 'Atomic neutral density'
                    title_name = 'Atomic neutral density [$m^{-3}$]'
                    dat = neuden_dat
                    vmin_0 = 1E+14
                    vmax_0 = 4E+17

                elif plot_name == 'denm':

                    datname = 'Molecular neutral density'
                    title_name = 'Molecular neutral density'
                    dat = denm_dat
                    vmin_0 = 1E+14
                    vmax_0 = 8E+17

                elif plot_name == 'ne':

                    datname = 'Electron density'
                    title_name = 'Electron density'
                    dat = ne_dat
                    vmin_0 = 1E+18
                    vmax_0 = 8E+19

                elif plot_name == 'molsource':

                    datname = 'Molecular neutral source'
                    title_name = 'Molecular neutral source'
                    dat = molsource_dat
                    vmin_0 = 1E+14
                    vmax_0 = 2E+17

                elif plot_name == 'mol_radial_flux':

                    datname = 'Molecular radial flux'
                    title_name = 'Molecular radial flux'
                    dat = rfluxm_dat
                    vmin_0 = 1E+19
                    vmax_0 = 1E+21

                elif plot_name == 'atomic_radial_flux':

                    datname = 'Atomic radial flux'
                    title_name = 'Atomic radial flux'

                    dat = rfluxa_dat
                    vmin_0 = 1E+19
                    vmax_0 = 1E+21

                elif plot_name == 'atomic_density_change':

                    datname = 'atomic_density_change'
                    title_name = 'Atomic neutral density change (%)'

                    norm_inc_dat = return_increase(
                        high_dat=neuden_dat, std_dat=nd_std_dat)

                    dat = norm_inc_dat
                    vmin_0 = -100
                    vmax_0 = 200

                elif plot_name == 'molecular_density_change':

                    datname = 'molecular_density_change'
                    title_name = 'molecular_density_change (%)'

                    norm_inc_dat = return_increase(
                        high_dat=denm_dat, std_dat=nm_std_dat)

                    dat = norm_inc_dat
                    vmin_0 = -100
                    vmax_0 = 150

                elif plot_name == 'ion_source_change':

                    datname = 'ion_source_change'
                    title_name = 'Ion source rate change (%)'

                    norm_inc_dat = return_increase(
                        high_dat=sx, std_dat=sx_std)

                    dat = norm_inc_dat
                    # vmin_0 = -150
                    # vmax_0 = 350

                    vmin_0 = -100
                    vmax_0 = 150

                if label_type == 'angle':

                    st = int(pol_loc_list[0])
                    ed = int(pol_loc_list[-1]) + 1

                    sx_limit = dat[st:ed, ins:out]
                    psi_limit = psi_map[st:ed, ins:out]
                    Ra_limit = ang_map[:, ins:out]

                    self.cpmc.twcontourp(plot_2dval=sx_limit, R_coord=Ra_limit,
                                         Z_coord=psi_limit, quantity=title_name, axs=axs,
                                         norm_type=norm_type, lv=40, vmin=vmin_0, vmax=vmax_0)

                    axs.set_ylabel("$\psi_N$")
                    axs.set_xlabel("poloidal angle")

                    axs.axvline(x=0, color='black', lw=3, ls='-', label='LFS')
                    axs.axvline(x=180, color='brown',
                                lw=3, ls='-', label='HFS')

                    axs.axhline(y=1, color='gray', lw=3,
                                ls='--', label='separatrix')
                    # axs.legend(loc='best')

                elif label_type == 'index':

                    sx_limit = dat[:, ins:out]
                    psi_limit = psi_map[:, ins:out]
                    Ra_limit = Ra_map[:, ins:out]

                    self.cpmc.twcontourp(plot_2dval=sx_limit, R_coord=Ra_limit,
                                         Z_coord=psi_limit, quantity=title_name, axs=axs,
                                         norm_type=norm_type, lv=40, vmin=vmin_0, vmax=vmax_0)

                    axs.set_ylabel("$\psi_N$")
                    axs.set_xlabel("poloidal index")

                elif label_type == 'arcmap':

                    st = int(pol_loc_list[0])
                    ed = int(pol_loc_list[-1]) + 1

                    sx_limit = dat[st:ed, ins:out]
                    psi_limit = arc_map[:, ins:out]
                    Ra_limit = ang_map[:, ins:out]

                    self.cpmc.twcontourp(plot_2dval=sx_limit, R_coord=Ra_limit,
                                         Z_coord=psi_limit, quantity=title_name, axs=axs,
                                         norm_type=norm_type, lv=40, vmin=vmin_0, vmax=vmax_0)

                    axs.set_ylabel("arc_length")
                    axs.set_xlabel("poloidal angle")
