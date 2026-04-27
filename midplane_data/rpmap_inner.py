# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:39:19 2026

@author: ychuang
"""


from load_coordinate.load_coord_method import load_coordgeo_method
from fit_data.fitting_method import fit_method_collection
from midplane_data.SOLPSplotter_PRmap import RP_mapping
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np


class Inner_RPmapping:

    def __init__(self, DF, data, lcm: load_coordgeo_method, fmc: fit_method_collection):

        self.DF = DF
        self.data = data
        self.fmc = fmc
        self.lcm = lcm

    def calc_RRsep_method(self, geo, radgrid, vertgrid, gfile_data, psiNinterp_function, plotRR, midplane_loc):

        pol_range = int(geo['nx'] + 2)
        rad_range = int(geo['ny'] + 2)

        RadLoc = radgrid
        VertLoc = vertgrid
        psiNinterp_RBS = psiNinterp_function

        # mag_axis_z = gfile_data['zmaxis']
        # mag_axis_z = 0.15
        # print(mag_axis_z)

        "Calculate midplane R for Z at the magnetic axis"

        "Calculate weight"

        DEV = self.DF.DEV

        if DEV == 'mast':
            mag_axis_z = gfile_data['zmaxis']

            crup = RadLoc[:, 58]
            crlow = RadLoc[:, 60]
            czup = VertLoc[:, 58]
            czlow = VertLoc[:, 60]

            average_pair = (58, 60)

        elif DEV == 'mastu':

            if midplane_loc == 'maxis':

                mag_axis_z = gfile_data['zmaxis']

                crup = RadLoc[:, 72]
                crlow = RadLoc[:, 74]
                czup = VertLoc[:, 72]
                czlow = VertLoc[:, 74]

                average_pair = (72, 74)

            elif midplane_loc == 'TSmeasure':

                mag_axis_z = 0.015

                crup = RadLoc[:, 72]
                crlow = RadLoc[:, 74]
                czup = VertLoc[:, 72]
                czlow = VertLoc[:, 74]

                average_pair = (72, 74)

        else:
            print("please check the device input")

        weight_mid = np.zeros(rad_range)
        for x in range(rad_range):
            weight_mid[x] = (mag_axis_z - czlow[x]) / (czup[x] - czlow[x])

        mid_choice = np.zeros((rad_range, 4))
        mid_choice[:, 0] = czup
        mid_choice[:, 1] = czlow
        mid_choice[:, 2] = crup
        mid_choice[:, 3] = crlow

        mid_R = np.zeros(rad_range)
        for xa in range(rad_range):
            mid_R[xa] = weight_mid[xa]*crup[xa] + \
                (1 - weight_mid[xa])*crlow[xa]

        mid_Z = np.zeros(rad_range)
        for xb in range(rad_range):
            mid_Z[xb] = weight_mid[xb]*czup[xb] + \
                (1 - weight_mid[xb])*czlow[xb]

        psi_solps_mid = np.zeros(rad_range)

        # psi_solps_cp = psiNinterp_RBS(crloc, czloc)
        for i in range(rad_range):
            psi_mid = psiNinterp_RBS(mid_R[i], mid_Z[i])
            psi_solps_mid[i] = psi_mid

        sep_index_dic = self.calc_sep_index(psi=psi_solps_mid, sep_value=1)
        sep_index_high = int(sep_index_dic['index'][1])

        print('index high number is {}'.format(sep_index_high))

        weight_psi = (psi_solps_mid[sep_index_high] - 1)/(
            psi_solps_mid[sep_index_high] - psi_solps_mid[sep_index_high - 1])

        R_sep = (1 - weight_psi) * \
            mid_R[sep_index_high] + weight_psi*mid_R[sep_index_high - 1]

        print('rsep is {:.3f}'.format(R_sep))

        # psi_high = psi_solps_mid[sep_index_high]
        # psi_low = psi_solps_mid[sep_index_high -1]

        # print('psi_high is {:.3f}'.format(psi_high))
        # print('psi_low is {:.3f}'.format(psi_low))

        R_Rsep = mid_R - R_sep

        dsa_psi_func = interp1d(R_Rsep, psi_solps_mid,
                                kind='quadratic', fill_value='extrapolate')

        midR_psi_func = interp1d(mid_R, psi_solps_mid,
                                 kind='quadratic', fill_value='extrapolate')

        rsep_fit = midR_psi_func(R_sep)

        print('rsep_mid return psiN value: {:.3f}'.format(rsep_fit))

        midplane_dic = {'weight': weight_mid, 'mid_choice': mid_choice,
                        'mid_R': mid_R, 'mid_Z': mid_Z,
                        'psi_solps_mid': psi_solps_mid,
                        'weight_psi': weight_psi, 'R_Rsep': R_Rsep,
                        'dsa_psi_func': dsa_psi_func,
                        'midR_psi_func': midR_psi_func,
                        'average_pair': average_pair}

        psi_dsa_dic = self.fmc.dsa_psi_fit(dsa=R_Rsep, psi=psi_solps_mid)

        psi_dsa_ratio = psi_dsa_dic['dsa_psi_fitcoe'][0]

        if DEV == 'mast':
            pol_list = [57, 58, 59, 60, 61]

        elif DEV == 'mastu':
            pol_list = [70, 71, 72, 73, 74]

        if plotRR:
            plt.figure()
            for in_pol in pol_list:
                crloc = RadLoc[:, int(in_pol)]
                czloc = VertLoc[:, int(in_pol)]
                plt.plot(crloc, czloc, color='g', label='R&Zlocation')
            plt.plot(mid_R, mid_Z, color='r', label='R_Rsep')
            plt.xlabel('R: [m]')
            plt.ylabel('Z: [m]')

        return midplane_dic, psi_dsa_ratio, sep_index_high

    def calc_RRsep(self, plotRR, plot_psi_dsa_align, midplane_loc):

        withshift = self.DF.withshift
        withseries = self.DF.withseries
        DEV = self.DF.DEV

        if withshift == False and withseries == False:
            b2fgeo = self.data['b2fgeo']
            radloc = self.data['grid']['RadLoc']
            vertloc = self.data['grid']['VertLoc']
            gfile = self.data['gfile']['g']
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']

            midplane_dic, psi_dsa_ratio, sep_index_high = self.calc_RRsep_method(geo=b2fgeo,
                                                                                 radgrid=radloc, vertgrid=vertloc, gfile_data=gfile, psiNinterp_function=psiNinterp_RBS,
                                                                                 plotRR=plotRR, midplane_loc=midplane_loc)

            if DEV == 'mast':

                self.data['DefaultSettings']['psi_dsa'] = psi_dsa_ratio
                self.data['midplane_calc'] = midplane_dic
                self.data['DefaultSettings']['sep_index_RRsep'] = sep_index_high

            elif DEV == 'mastu':

                if midplane_loc == 'maxis':

                    self.data['DefaultSettings']['psi_dsa_maxis'] = psi_dsa_ratio
                    self.data['midplane_calc_maxis'] = midplane_dic
                    self.data['DefaultSettings']['sep_index_RRsep_maxis'] = sep_index_high

                elif midplane_loc == 'TSmeasure':

                    self.data['DefaultSettings']['psi_dsa'] = psi_dsa_ratio
                    self.data['midplane_calc'] = midplane_dic
                    self.data['DefaultSettings']['sep_index_RRsep'] = sep_index_high

        elif withshift == True and withseries == False:

            sep_index_dic = {}
            midplane_series_dic = {}
            psi_dsa_dic = {}
            for aa in self.data['dircomp']['multi_shift']:

                b2fgeo = self.data['b2fgeo'][aa]
                radloc = self.data['grid']['RadLoc'][aa]
                vertloc = self.data['grid']['VertLoc'][aa]
                gfile = self.data['gfile']['g']
                psiNinterp_RBS = self.data['gfile']['gcomp'][aa]['interp_dic']['RBS']

                midplane_dic, psi_dsa_ratio, sep_index_high = self.calc_RRsep_method(geo=b2fgeo,
                                                                                     radgrid=radloc, vertgrid=vertloc, gfile_data=gfile, psiNinterp_function=psiNinterp_RBS,
                                                                                     plotRR=plotRR, midplane_loc=midplane_loc)

                midplane_series_dic[aa] = midplane_dic
                psi_dsa_dic[aa] = psi_dsa_ratio
                sep_index_dic[aa] = sep_index_high

            if DEV == 'mast':

                self.data['DefaultSettings']['psi_dsa'] = psi_dsa_dic
                self.data['midplane_calc'] = midplane_series_dic
                self.data['DefaultSettings']['sep_index_RRsep'] = sep_index_dic

            elif DEV == 'mastu':

                if midplane_loc == 'maxis':

                    self.data['DefaultSettings']['psi_dsa_maxis'] = psi_dsa_dic
                    self.data['midplane_calc_maxis'] = midplane_series_dic
                    self.data['DefaultSettings']['sep_index_RRsep_maxis'] = sep_index_dic

                elif midplane_loc == 'TSmeasure':

                    self.data['DefaultSettings']['psi_dsa'] = psi_dsa_dic
                    self.data['midplane_calc'] = midplane_series_dic
                    self.data['DefaultSettings']['sep_index_RRsep'] = sep_index_dic

            if plot_psi_dsa_align:
                color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                             'dot7': 'blue', 'one': 'purple'}
                plt.figure()
                for in_sh in self.data['dircomp']['multi_shift']:
                    # plt.plot(midplane_dic[in_sh]['psi_solps_mid'], midplane_dic[in_sh]['R_Rsep'], 'o-',
                    #   label= 'R-Rsep:[m] with modify {} m'.format(self.data['dircomp']['shift_dic'][in_sh]))
                    plt.plot(midplane_dic[in_sh]['psi_solps_mid'], midplane_dic[in_sh]['R_Rsep'], 'o-',
                             color=color_dic[in_sh])
                plt.xlabel('psiN')
                plt.title('psiN R-Rsep plot')
                # plt.legend()

        elif withshift == False and withseries == True:
            den_rep = list(self.data['dircomp']['Attempt'].keys())[0]

            b2fgeo = self.data['b2fgeo']
            radloc = self.data['grid']['RadLoc']
            vertloc = self.data['grid']['VertLoc']
            gfile = self.data['gfile']['g']
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']

            midplane_dic, psi_dsa_ratio, sep_index_high = self.calc_RRsep_method(geo=b2fgeo,
                                                                                 radgrid=radloc, vertgrid=vertloc, gfile_data=gfile, psiNinterp_function=psiNinterp_RBS,
                                                                                 plotRR=plotRR, midplane_loc=midplane_loc)

            if DEV == 'mast':

                self.data['DefaultSettings']['psi_dsa'] = psi_dsa_ratio
                self.data['midplane_calc'] = midplane_dic
                self.data['DefaultSettings']['sep_index_RRsep'] = sep_index_high

            elif DEV == 'mastu':

                if midplane_loc == 'maxis':

                    self.data['DefaultSettings']['psi_dsa_maxis'] = psi_dsa_ratio
                    self.data['midplane_calc_maxis'] = midplane_dic
                    self.data['DefaultSettings']['sep_index_RRsep_maxis'] = sep_index_high

                elif midplane_loc == 'TSmeasure':

                    self.data['DefaultSettings']['psi_dsa'] = psi_dsa_ratio
                    self.data['midplane_calc'] = midplane_dic
                    self.data['DefaultSettings']['sep_index_RRsep'] = sep_index_high

        elif self.withshift == True and self.withseries == True:
            print('calc_RRsep function is not there yet!')

        else:
            print('calc_RRsep function has a bug')
