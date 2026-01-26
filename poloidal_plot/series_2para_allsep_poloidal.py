# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 18:40:31 2026

@author: ychuang
"""


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnchoredText
import matplotlib.pylab as pylab
from midplane_data.cal_sep_length import seplength_calculator


class series_2para_allsep_polplot:

    def __init__(self, DF, data, sc: seplength_calculator):

        self.DF = DF
        self.data = data
        self.sc = sc

    def init_param_lists(self, parameter: str) -> dict:
        """Return a dict of empty lists for the requested parameter."""
        mapping = {
            "nete": {"ne": [], "te": []},
            "teti": {"te": [], "ti": []},
            "ndS":  {"nd": [], "s":  []},
            "ndnm": {"nd": [], "nm": []},
            "ndapf": {"nd": [], "abs_fnaxs": []},
            "rpa": {"rfluxa": [], "pfluxa": []},
            "rpm": {"rfluxm": [], "pfluxm": []},
            "pam": {"pfluxa": [], "pfluxm": []},
            "ppa": {"fnaxs": [], "pfluxa": []},
        }
        try:
            return mapping[parameter]
        except KeyError:
            raise ValueError(f"Unknown parameter: {parameter}")

    def append_param_values(self, parameter: str, lists: dict, ii: int, rad_loc: int, data: dict) -> None:
        """
        Append values for one poloidal index ii at rad_loc into the provided lists dict.
        data must contain the arrays referenced by the parameter.
        """
        config = {
            # parameter: ( (list_keys...), (array_keys...) )
            "nete": (("ne", "te"), ("ne_dat", "te_pro")),
            "teti": (("te", "ti"), ("te_pro", "ti_pro")),
            "ndS":  (("nd", "s"),  ("neuden", "source")),
            "ndnm": (("nd", "nm"), ("neuden", "molden")),
            "ndapf": (("nd", "abs_fnaxs"), ("neuden", "abs_fnaxs")),
            "rpa": (("rfluxa", "pfluxa"), ("rfluxa", "pfluxa")),
            "rpm": (("rfluxm", "pfluxm"), ("rfluxm", "pfluxm")),
            "pam": (("pfluxa", "pfluxm"), ("pfluxa", "pfluxm")),
            "ppa": (("fnaxs", "pfluxa"), ("fnaxs", "pfluxa")),
        }

        try:
            list_keys, array_keys = config[parameter]
        except KeyError:
            raise ValueError(f"Unknown parameter: {parameter}")

        for lk, ak in zip(list_keys, array_keys):
            lists[lk].append(data[ak][ii, rad_loc])

    def plot_parameter(self, axs, ang_list, parameter, lists, aa, lins_map, cl_dic):
        """
        axs: array-like of matplotlib axes
        ang_list: x values
        parameter: e.g. 'nete', 'teti', 'ndS', 'ndnm'
        lists: dict of y-lists returned from collect/init (e.g. {'ne': [...], 'te': [...]})
        aa: case label (used for color/linestyle and legend label)
        lins_map: dict mapping aa -> linestyle code
        cl_dic: dict mapping aa -> color
        """

        plot_cfg = {
            "nete": {
                "series": [("ne", 0, r"$n_e$"), ("te", 1, r"$T_e$")],
            },
            "teti": {
                "series": [("te", 0, r"$T_e$"), ("ti", 1, r"$T_i$")],
            },
            "ndS": {
                "series": [("nd", 0, r"$n_D$"), ("s", 1, r"Source")],
            },
            "ndnm": {
                "series": [("nd", 0, r"$n_D$"), ("nm", 1, r"$n_m$")],
            },
            "ndapf": {
                "series": [("nd", 0, r"$n_D$"), ("abs_fnaxs", 1, r"$\Gamma_\theta$")],
            },
            "rpa": {
                "series": [("rfluxa", 0, r"$\Gamma_{r, n_D}$"), ("pfluxa", 1, r"$\Gamma_{p, n_D}$")],
            },
            "rpm": {
                "series": [("rfluxm", 0, r"$\Gamma_{r, n_m}$"), ("pfluxm", 1, r"$\Gamma_{p, n_m}$")],
            },
            "pam": {
                "series": [("pfluxa", 0, r"$\Gamma_{p, n_D}$"), ("pfluxm", 1, r"$\Gamma_{p, n_m}$")],
            },
            "ppa": {
                "series": [("fnaxs", 0, r"$\Gamma_{p, ion}$"), ("pfluxa", 1, r"$\Gamma_{p, n_D}$")],
            },
        }

        if parameter not in plot_cfg:
            raise ValueError(f"Unknown parameter: {parameter}")

        ls = lins_map[aa] + "-"   # keep your convention
        color = cl_dic[aa]

        for key, ax_i, ylabel in plot_cfg[parameter]["series"]:
            axs[ax_i].plot(
                ang_list,
                lists[key],
                ls,
                lw=1.5,
                color=color,
                label=f"{aa}" if ax_i == 0 else None,  # legend only once
            )
            axs[ax_i].set_ylabel(ylabel)

    def series_2para_allsep_method(self, iterlist, cl_dic, parameter, ang_list, pol_list, rad_loc):

        fig, axs = plt.subplots(2, 1)

        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']

        marker = ['o', 's', 'x', '^', 'v', '*', '.', 'o', 's']
        lins_map = dict(zip(iterlist, marker))

        plt.subplots_adjust(hspace=.0)
        anchored_text_1 = AnchoredText('{}'.format('Electron density [$m^{-3}$]'),
                                       loc='upper center')
        anchored_text_2 = AnchoredText('{}'.format('Electron temperature [eV]'),
                                       loc='lower right')
        anchored_text_3 = AnchoredText('{}'.format('Ion temperature [eV]'),
                                       loc='lower right')

        rad_grid = self.data['grid']['RadLoc']
        vert_grid = self.data['grid']['VertLoc']

        for aa in iterlist:

            b2fstate = self.data['b2fstate'][aa]
            ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
            Te_J = b2fstate['te'][1:nx+1, 1:ny+1]
            Ti_J = b2fstate['ti'][1:nx+1, 1:ny+1]
            ev = 1.6021766339999999 * pow(10, -19)
            te_pro = Te_J / ev
            ti_pro = Ti_J / ev

            fnax = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
            hz = self.data['b2wdat'][aa]['hz'][1:nx+1, 1:ny+1]
            hy = self.data['b2wdat'][aa]['hy'][1:nx+1, 1:ny+1]
            tor_area = np.multiply(hz, hy)
            fnaxs = np.divide(fnax, tor_area)
            abs_fnaxs = abs(fnaxs)

            s_term = self.data['b2wdat'][aa]['b2npc_sna'][0][1:nx+1, 1:ny+1]
            vol = self.data['b2wdat'][aa]['vol'][1:nx+1, 1:ny+1]
            source = np.divide(s_term, vol)

            neuden = self.data['ft44'][aa]['dab2'][:, :, 0]
            molden = self.data['ft44'][aa]['dmb2'][:, :, 0]
            rfluxa = self.data['ft44'][aa]['rfluxa'][:, :, 0]
            rfluxm = self.data['ft44'][aa]['rfluxm'][:, :, 0]
            pfluxa = self.data['ft44'][aa]['pfluxa'][:, :, 0]
            pfluxm = self.data['ft44'][aa]['pfluxm'][:, :, 0]

            data = dict(ne_dat=ne_dat, te_pro=te_pro, ti_pro=ti_pro, neuden=neuden,
                        source=source, molden=molden, rfluxa=rfluxa, rfluxm=rfluxm, pfluxa=pfluxa,
                        pfluxm=pfluxm, fnaxs=fnaxs, abs_fnaxs=abs_fnaxs)

            R_loc = rad_grid[rad_loc, 1:97]

            Z_loc = vert_grid[rad_loc, 1:97]

            arclength, interpolated_points = self.sc.calc_sepx_length(cr=R_loc,
                                                                      cz=Z_loc, plot=False)

            lists = self.init_param_lists(parameter)

            index_list = np.linspace(0, 95, 96)

            for ii in index_list:
                ii = int(ii)
                self.append_param_values(
                    parameter, lists, ii, rad_loc, data)

            # self.plot_parameter(ang_list=arclength, axs=axs, parameter=parameter,
            #                     lists=lists, aa=aa, lins_map=lins_map, cl_dic=cl_dic)
            axs[0].legend()
            axs[1].set_xlabel(
                'separatrix arclength distance from inner to outer target[m]')
            axs[0].legend(loc='best')

        # axs[0].add_artist(anchored_text_1)
        # axs[1].add_artist(anchored_text_3)

        axs[0].axvline(x=0, color='gray', lw=3, ls='-')
        axs[1].axvline(x=0, color='gray', lw=3, ls='-', label='Inner target')
        axs[0].axvline(x=arclength[24], color='black', lw=3, ls='--')
        axs[1].axvline(x=arclength[24], color='black',
                       lw=3, ls='--', label='xpoint')
        axs[0].axvline(x=arclength[36], color='brown', lw=3, ls='-')
        axs[1].axvline(x=arclength[36], color='brown',
                       lw=3, ls='-', label='Inner midplane')
        axs[0].axvline(x=arclength[59], color='darkolivegreen', lw=3, ls='-')
        axs[1].axvline(x=arclength[59], color='darkolivegreen',
                       lw=3, ls='-', label='Outer midplane')
        axs[0].axvline(x=arclength[72], color='black', lw=3, ls='--')
        axs[1].axvline(x=arclength[72], color='black', lw=3, ls='--')
        axs[0].axvline(x=arclength[-1], color='steelblue', lw=3, ls='-')
        axs[1].axvline(x=arclength[-1], color='steelblue',
                       lw=3, ls='-', label='Outer target')

        if parameter == 'ndnm':

            axs[0].set_yscale('log')
            axs[1].set_yscale('log')

        elif parameter == 'ndS':

            axs[0].set_yscale('log')
            axs[1].set_yscale('log')

        elif parameter == 'ndapf':

            axs[0].set_yscale('log')
            axs[1].set_yscale('log')
            axs[0].set_ylim(7E+15, 1E+18)
            axs[1].set_ylim(1E+21, 1E+23)

        axs[1].legend(loc='best')
        axs[0].grid(True)
        axs[1].grid(True)

    def series_2para_allsepplot(self, pol_list, parameter, rad_loc):

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == False and withseries == True:

            labels = list(self.data['dircomp']['Attempt'].keys())
            colors = pylab.cm.gnuplot2(np.linspace(
                0, 1, len(labels) + 1))
            color_map = dict(zip(labels, colors))

            ang_list = self.data['angle']['angle_list']

            self.series_2para_allsep_method(iterlist=labels, cl_dic=color_map, ang_list=ang_list, rad_loc=rad_loc,
                                            pol_list=pol_list, parameter=parameter)
