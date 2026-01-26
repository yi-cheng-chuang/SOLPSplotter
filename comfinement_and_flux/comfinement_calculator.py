# -*- coding: utf-8 -*-
"""
Created on Sun Jan 11 21:48:47 2026

@author: ychuang
"""


import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np


class comfinement_calculation:

    def __init__(self, DF, data):

        self.DF = DF
        self.data = data

    def particle_balance_method(self, iterlist):

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == False and withseries == True:

            marker = ['o', 's', 'x', '^', 'v', '*', '.', 'o', 's']
            lins_map = dict(zip(iterlist, marker))

            color_list = ['red', 'salmon', 'orange', 'lime', 'green', 'darkgreen',
                          'cyan', 'deepskyblue', 'navy', 'purple']

            rfcore_list = np.zeros(len(iterlist))
            rfsep_list = np.zeros(len(iterlist))
            scomfine_list = np.zeros(len(iterlist))
            error_list = np.zeros(len(iterlist))
            mvalue_list = np.zeros(len(iterlist))
            error_percentage_list = np.zeros(len(iterlist))
            rhefcore_list = np.zeros(len(iterlist))
            rhifcore_list = np.zeros(len(iterlist))
            norm_p_list = np.zeros(len(iterlist))
            norm_neu_list = np.zeros(len(iterlist))
            fnay_sol_list = np.zeros(len(iterlist))
            fnax_tar_abs_list = np.zeros(len(iterlist))
            s_sol_list = np.zeros(len(iterlist))
            neu_sol_list = np.zeros(len(iterlist))
            tar_recycle_list = np.zeros(len(iterlist))
            tar_imping_list = np.zeros(len(iterlist))

            nx = self.data['b2fgeo']['nx']
            ny = self.data['b2fgeo']['ny']

            std_b2fstate = self.data['b2fstate']['1.0']
            std_ne_dat = std_b2fstate['ne'][1:nx+1, 1:ny+1]
            std_vol = self.data['b2wdat']['1.0']['vol'][1:nx+1, 1:ny+1]
            std_particle = np.multiply(std_ne_dat, std_vol)

            std_neu = self.data['ft44']['1.0']['dab2'][:, :, 0]
            std_allneu = np.multiply(std_neu, std_vol)

            for ik, aa in enumerate(iterlist):

                fnay = self.data['b2wdat'][aa]['b2npc_fnays'][0][1:97, 1:37]
                fnax = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                # fnax = self.data['b2wdat'][aa]['b2npc_fnaxs'][0][1:97, 1:37]
                pfluxa = self.data['ft44'][aa]['pfluxa'][:, :, 0]
                pfluxm = self.data['ft44'][aa]['pfluxm'][:, :, 0]
                rfluxa = self.data['ft44'][aa]['rfluxa'][:, :, 0]
                rfluxm = self.data['ft44'][aa]['rfluxm'][:, :, 0]

                hz = self.data['b2wdat'][aa]['hz'][1:nx+1, 1:ny+1]
                hy = self.data['b2wdat'][aa]['hy'][1:nx+1, 1:ny+1]
                hx = self.data['b2wdat'][aa]['hx'][1:nx+1, 1:ny+1]
                tor_area = np.multiply(hz, hy)
                pol_area = np.multiply(hz, hx)
                pfluxas = np.multiply(pfluxa, tor_area)
                pfluxms = np.multiply(pfluxm, tor_area)
                rfluxas = np.multiply(rfluxa, pol_area)
                rfluxms = np.multiply(rfluxm, pol_area)

                neu_radial_flux = sum(rfluxa[:, -1]) + sum(rfluxm[:, -1])*2
                # neu_radial_flux = sum(rfluxms[:, -1])*2

                fnay_core = sum(fnay[24:72, 0])
                fnay_sep = sum(fnay[24:72, 18])
                fnay_sol = sum(fnay[:, -1])
                fnax_tar_abs = (
                    sum(abs(fnax[0, :])) + sum(abs(fnax[-1, :])))*0.01

                fnax_tar = sum(abs(fnax[0, :])) + sum(abs(fnax[-1, :]))

                neu_tar_abs = (
                    sum(abs(pfluxas[0, :])) + sum(abs(pfluxas[-1, :])))*0.01/0.99

                neu_tar = sum(abs(pfluxas[0, :])) + sum(abs(pfluxas[-1, :]))

                mol_tar_abs = (
                    sum(abs(pfluxms[0, :])) + sum(abs(pfluxms[-1, :])))*0.02/0.99

                mol_tar = sum(abs(pfluxms[0, :])) + sum(abs(pfluxms[-1, :]))

                total_tar_abs = fnax_tar_abs + neu_tar_abs + mol_tar_abs
                # total_tar_abs = fnax_tar_abs + neu_tar_abs

                total_tar_recycle = mol_tar + neu_tar

                fhey = self.data['b2wdat'][aa]['b2nph9_fhey'][1:97, 1:37]

                fhey_core = sum(fhey[24:72, 0])

                fhiy = self.data['b2wdat'][aa]['b2nph9_fhiy'][1:97, 1:37]

                fhiy_core = sum(fhiy[24:72, 0])

                source = self.data['b2wdat'][aa]['b2npc_sna'][0][1:nx+1, 1:ny+1]

                s_comfine = np.sum(source[24:72, :18])
                s_sol = np.sum(source[24:72, 19:])

                b2fstate = self.data['b2fstate'][aa]
                ne_dat = b2fstate['ne'][1:nx+1, 1:ny+1]
                vol = self.data['b2wdat'][aa]['vol'][1:nx+1, 1:ny+1]
                particle = np.multiply(ne_dat, vol)
                comfine_particle = np.sum(particle[24:72, :18])

                comfine_std_particle = np.sum(std_particle[24:72, :18])

                norm_p = comfine_particle / comfine_std_particle

                neu = self.data['ft44'][aa]['dab2'][:, :, 0]
                all_neu = np.multiply(neu, vol)
                comfine_neu = np.sum(all_neu[24:72, :18])

                comfine_std_neu = np.sum(std_allneu[24:72, :18])

                norm_neu = comfine_neu / comfine_std_neu

                error = fnay_sep - fnay_core - s_comfine

                mvalue = max(fnay_core, fnay_sep, s_comfine)

                error_percentage = error * 100 / fnay_core

                rfcore_list[ik] = fnay_core
                rfsep_list[ik] = fnay_sep
                scomfine_list[ik] = s_comfine
                error_list[ik] = error
                mvalue_list[ik] = mvalue
                error_percentage_list[ik] = error_percentage
                rhefcore_list[ik] = fhey_core
                rhifcore_list[ik] = fhiy_core
                norm_p_list[ik] = norm_p
                norm_neu_list[ik] = norm_neu
                fnay_sol_list[ik] = fnay_sol
                fnax_tar_abs_list[ik] = total_tar_abs
                s_sol_list[ik] = s_sol
                neu_sol_list[ik] = abs(neu_radial_flux)
                tar_recycle_list[ik] = total_tar_recycle
                tar_imping_list[ik] = fnax_tar

            rounded = [round(float(x), 1) for x in iterlist]

            def round_dat(dat):

                dat_rounded = np.array([float(f"{x:.4g}") for x in dat])

                return dat_rounded

            rfcore_round = round_dat(dat=rfcore_list)
            rhefcore_round = round_dat(dat=rhefcore_list)
            rhifcore_round = round_dat(dat=rhifcore_list)

            pfcore = 5.512e+20
            hfcore = 4.115e+5

            print("print core particle flux")
            print(rfcore_round)
            print("print core electron heat flux")
            print(rhefcore_round)
            print("print core ion heat flux")
            print(rhifcore_round)

            print("print core particle flux ratio")
            print(rfcore_round / pfcore)
            print("print core electron heat flux ratio")
            print(rhefcore_round / hfcore)
            print("print core ion heat flux ratio")
            print(rhifcore_round / hfcore)

            fig, axs = plt.subplots()

            axs.plot(rounded, rfsep_list, '-o', color=color_list[3],
                     label='sep particle flux')

            axs.plot(rounded, scomfine_list, '-o', color=color_list[5],
                     label='particle source')

            axs.set_yscale("log")
            axs.set_title("flux element plot")
            axs.legend()

            fig, axs = plt.subplots()

            axs.plot(rounded, rfcore_list, '-o', color=color_list[0],
                     label='core particle flux')

            axs.plot(rounded, error_list, '-o', color=color_list[3],
                     label='error')

            axs.set_yscale("log")
            axs.set_title("compare error plot")
            axs.legend()

            fig, axs = plt.subplots()

            axs.plot(rounded, error_percentage_list, '-o', color=color_list[0],
                     label='error percentage')

            axs.set_title("error percentage plot")
            axs.legend()

            fig, axs = plt.subplots()

            axs.plot(rounded, rhefcore_list, '-o', color=color_list[0],
                     label='electron heat flux')

            axs.plot(rounded, rhifcore_list, '-o', color=color_list[3],
                     label='ion heat flux')

            axs.set_title("heat flux plot")
            axs.legend()

            fig, axs = plt.subplots()

            axs.plot(rounded, norm_p_list, '-o', color=color_list[0],
                     label='norm p')

            axs.plot(rounded, norm_neu_list, '-o', color=color_list[3],
                     label='norm neu')

            axs.set_title("particle compare")
            axs.legend()

            fig, axs = plt.subplots()

            axs.plot(rounded, fnay_sol_list, '-o', color=color_list[0],
                     label='SOL fnay')

            axs.plot(rounded, s_sol_list, '-o', color=color_list[3],
                     label='SOL source')

            axs.set_title("source and boundary fnay at SOL")
            axs.legend()

            fig, axs = plt.subplots()

            axs.plot(rounded, fnax_tar_abs_list, '-o', color=color_list[0],
                     label='absorbtion at the target')

            axs.plot(rounded, rfcore_list, '-o', color=color_list[3],
                     label='core particle flux')

            axs.set_yscale("log")
            axs.set_title("absorbtion at the target")
            axs.legend()

            fig, axs = plt.subplots()

            axs.plot(rounded, fnay_sol_list, '-o', color=color_list[0],
                     label='SOL fnay')

            axs.plot(rounded, neu_sol_list, '-o', color=color_list[3],
                     label='sol radial neutral flux')

            # axs.set_yscale("log")
            axs.set_title("compare radial flux")
            axs.legend()

            fig, axs = plt.subplots()

            axs.plot(rounded, tar_imping_list, '-o', color=color_list[0],
                     label='target fnax')

            axs.plot(rounded, tar_recycle_list, '-o', color=color_list[3],
                     label='target neutral flux')

            # axs.set_yscale("log")
            axs.set_title("compare target flux")
            axs.legend()

            fig, axs = plt.subplots()

            axs.plot(rounded, tar_imping_list - tar_recycle_list, '-o', color=color_list[0],
                     label='target flux diff')

            axs.plot(rounded, rfcore_list, '-o', color=color_list[3],
                     label='core particle flux')

            # axs.set_yscale("log")
            axs.set_title("compare target flux 2")
            axs.legend()

    def particle_balance(self):

        withshift = self.DF.withshift
        withseries = self.DF.withseries

        if withshift == False and withseries == True:

            # series_flag = self.DefaultSettings['series_flag']

            labels = list(self.data['dircomp']['Attempt'].keys())
            # colors = pylab.cm.gnuplot2(np.linspace(
            #     0, 1, len(labels) + 1))
            # color_map = dict(zip(labels, colors))

            self.particle_balance_method(iterlist=labels)
