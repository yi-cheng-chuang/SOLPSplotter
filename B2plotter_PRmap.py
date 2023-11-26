# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 17:27:01 2023

@author: user
"""

from B2plotter_load_data import load_data
import opacity_plot_method as opm
import matplotlib.pyplot as plt
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import numpy as np
import xarray as xr
import math


class RP_mapping(load_data):
    def __init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters):
        load_data.__init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters)
    
        
    def calc_RRsep(self, plotRR):
        
        if self.withshift == False and self.withseries == False:
        
            pol_range = int(self.data['b2fgeo']['nx'] + 2)
            rad_range = int(self.data['b2fgeo']['ny'] + 2)
            mag_axis_z = self.data['gfile']['g']['zmaxis']
            # print(mag_axis_z)
            
            Attempt = self.data['dircomp']['Attempt']
            DRT = self.data['dirdata']['outputdir']['Output']
            DRT2 = self.data['dirdata']['outputdir']['Output2']
            XDIM = self.data['b2fgeo']['nx'] + 2
            YDIM = self.data['b2fgeo']['ny'] + 2
            
            
            RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                        usecols = (3)).reshape((YDIM, XDIM))
            VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                          usecols = (3)).reshape((YDIM,XDIM))
            
            
            "Calculate midplane R for Z=0"
            
            "Calculate weight"
            
            crup = RadLoc[:, 58]
            crlow = RadLoc[:, 60]
            czup = VertLoc[:, 58]
            czlow = VertLoc[:, 60]
            
            pol_list = [52, 53, 54, 55, 56, 57]
            

            weight_mid = np.zeros(rad_range)
            for x in range(rad_range):
                weight_mid[x] = (mag_axis_z - czlow[x])/ (czup[x] - czlow[x])
            
            mid_choice = np.zeros((rad_range, 4))
            mid_choice[:, 0] = czup
            mid_choice[:, 1] = czlow
            mid_choice[:, 2] = crup
            mid_choice[:, 3] = crlow
            
            
            
            mid_R = np.zeros(rad_range)
            for xa in range(rad_range):
                mid_R[xa] = weight_mid[xa]*crup[xa] + (1 - weight_mid[xa])*crlow[xa]
            
            mid_Z = np.zeros(rad_range)
            for xb in range(rad_range):
                mid_Z[xb] = weight_mid[xb]*czup[xb] + (1 - weight_mid[xb])*czlow[xb]
                
            pol_list = [57, 58, 59, 60]
            
            
            
            psi_solps_mid = np.zeros(rad_range)
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            # psi_solps_cp = psiNinterp_RBS(crloc, czloc)
            for i in range(rad_range):
                psi_mid = psiNinterp_RBS(mid_R[i], mid_Z[i])
                psi_solps_mid[i] = psi_mid
            
            sep = []
            for in_sep in psi_solps_mid:
                if in_sep >= 1:
                    sep.append(list(psi_solps_mid).index(in_sep))
            
                    
            sep_index = sep[0]
            # print(sep)
            
            weight_psi = (psi_solps_mid[sep_index] - 1)/(psi_solps_mid[sep_index] -psi_solps_mid[sep_index -1])
            
            R_sep = (1 - weight_psi)*mid_R[sep_index] + weight_psi*mid_R[sep_index -1]
            R_Rsep = mid_R - R_sep
            
            midplane_dic = {'weight': weight_mid, 'mid_choice': mid_choice, 
                            'mid_R': mid_R, 'mid_Z': mid_Z, 
                            'psi_solps_mid': psi_solps_mid, 
                            'weight_psi': weight_psi, 'R_Rsep': R_Rsep}
            
            psi_dsa_dic = fm.dsa_psi_fit(dsa= R_Rsep, psi= psi_solps_mid)
            
            self.data['DefaultSettings']['psi_dsa'] = psi_dsa_dic['dsa_psi_fitcoe'][0]
                        
            self.data['midplane_calc'] = midplane_dic
            
            if plotRR:
                plt.figure(figsize=(7,7))
                for in_pol in pol_list:
                    crloc = RadLoc[:, int(in_pol)]
                    czloc = VertLoc[:, int(in_pol)]
                    plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
                plt.plot(mid_R, mid_Z, color = 'r', label= 'R_Rsep')
                plt.xlabel('R: [m]')
                plt.ylabel('Z: [m]')
            
            
               
        
        elif self.withshift == True and self.withseries == False:        
            pol_range_dic = {}
            rad_range_dic = {}
            SEP_dic = {}
            midplane_dic = {}
            psi_dsa_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                pol_range = int(self.data['b2fgeo'][aa]['nx'] + 2)
                # print('xdim is {}'.format(str(pol_range)))
                rad_range = int(self.data['b2fgeo'][aa]['ny'] + 2)
                # print('ydim is {}'.format(str(rad_range)))
                pol_range_dic[aa] = pol_range
                rad_range_dic[aa] = rad_range
                
                mag_axis_z = self.data['gfile']['g']['zmaxis']
                # print(mag_axis_z)
                
                Attempt = self.data['dircomp']['Attempt'][aa]
                DRT = self.data['dirdata']['infolderdir'][aa]['outputdir']['Output']
                DRT2 = self.data['dirdata']['infolderdir'][aa]['outputdir']['Output2']
                
                
                RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                            usecols = (3)).reshape((rad_range, pol_range))
                VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                              usecols = (3)).reshape((rad_range, pol_range))
                
                
                "Calculate midplane R for Z=0"
                
                "Calculate weight"
                
                if aa == 'one':
                    crup = RadLoc[:, 57]
                    crlow = RadLoc[:, 59]
                    czup = VertLoc[:, 57]
                    czlow = VertLoc[:, 59]
                else:
                    crup = RadLoc[:, 57]
                    crlow = RadLoc[:, 59]
                    czup = VertLoc[:, 57]
                    czlow = VertLoc[:, 59]
                
                
                weight_mid = np.zeros(rad_range)
                for x in range(rad_range):
                    weight_mid[x] = (czup[x] - mag_axis_z)/ (czup[x] - czlow[x])
                
                mid_choice = np.zeros((rad_range, 4))
                mid_choice[:, 0] = czup
                mid_choice[:, 1] = czlow
                mid_choice[:, 2] = crup
                mid_choice[:, 3] = crlow
                
                
                
                mid_R = np.zeros(rad_range)
                for xa in range(rad_range):
                    mid_R[xa] = (1 - weight_mid[xa])*crup[xa] + weight_mid[xa]*crlow[xa]
                
                mid_Z = np.zeros(rad_range)
                for xb in range(rad_range):
                    mid_Z[xb] = (1 - weight_mid[xb])*czup[xb] + weight_mid[xb]*czlow[xb]
                    
                
                pol_list = [57, 58, 59, 60]
                
                
                psi_solps_mid = np.zeros(rad_range)
                psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic'][aa]
                # psi_solps_cp = psiNinterp_RBS(crloc, czloc)
                for i in range(rad_range):
                    psi_mid = psiNinterp_RBS(mid_R[i], mid_Z[i])
                    psi_solps_mid[i] = psi_mid
                
                sep = []
                for in_sep in psi_solps_mid:
                    if in_sep >= 1:
                        sep.append(list(psi_solps_mid).index(in_sep))
                
                        
                sep_index = sep[0]
                SEP_dic[aa] = sep_index
                
                weight_psi = (psi_solps_mid[sep_index] - 1)/(psi_solps_mid[sep_index] -psi_solps_mid[sep_index -1])
                
                R_sep = (1 - weight_psi)*mid_R[sep_index] + weight_psi*mid_R[sep_index -1]
                R_Rsep = mid_R - R_sep
                
                midplane_dic[aa] = {'weight': weight_mid, 'mid_choice': mid_choice, 
                                'mid_R': mid_R, 'mid_Z': mid_Z, 
                                'psi_solps_mid': psi_solps_mid, 
                                'weight_psi': weight_psi, 'R_Rsep': R_Rsep}
                
                pd_dic = fm.dsa_psi_fit(dsa= R_Rsep, psi= psi_solps_mid)
                
                psi_dsa_dic[aa] = pd_dic['dsa_psi_fitcoe'][0]
                
                
                
                if plotRR:
                    plt.figure(figsize=(7,7))
                    for in_pol in pol_list:
                        crloc = RadLoc[:, int(in_pol)]
                        czloc = VertLoc[:, int(in_pol)]
                        plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
                    plt.plot(mid_R, mid_Z, color = 'r', label= 'R_Rsep')
                    plt.xlabel('R: [m]')
                    plt.ylabel('Z: [m]')
                    plt.title('{} R-Z plot'.format(aa))
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            plt.figure(figsize=(7,7))
            for in_sh in self.data['dircomp']['multi_shift']:
                # plt.plot(midplane_dic[in_sh]['psi_solps_mid'], midplane_dic[in_sh]['R_Rsep'], 'o-', 
                #   label= 'R-Rsep:[m] with modify {} m'.format(self.data['dircomp']['shift_dic'][in_sh]))
                plt.plot(midplane_dic[in_sh]['psi_solps_mid'], midplane_dic[in_sh]['R_Rsep'], 'o-', 
                  color= color_dic[in_sh])
            plt.xlabel('psiN')
            plt.title('psiN R-Rsep plot')
            # plt.legend()
                    
            self.data['DefaultSettings']['psi_dsa'] = psi_dsa_dic
            self.data['midplane_calc'] = midplane_dic
            self.data['DefaultSettings']['XDIM'] = pol_range_dic
            self.data['DefaultSettings']['YDIM'] = rad_range_dic
            self.data['DefaultSettings']['SEP'] = SEP_dic
            
        elif self.withshift == False and self.withseries == True:
        
            pol_range = self.data['b2fgeo']['nx'] + 2
            rad_range = self.data['b2fgeo']['ny'] + 2
            
            
            solps_dsa_dic = {}
            for aa in self.data['dircomp']['Attempt'].keys():
                solps_dsa_dic[aa] = lcm.read_dsa(self.data['dirdata']['simudir'][aa] + '/dsa')
            
            den_list = list(self.data['dircomp']['Attempt'].keys())
            aa = den_list[0]
            
            mag_axis_z = self.data['gfile']['g']['zmaxis']
            print(mag_axis_z)
            
            Attempt = self.data['dircomp']['Attempt'][aa]
            DRT = self.data['dirdata']['outputdir'][aa]['Output']
            # DRT2 = self.data['dirdata']['outputdir']['Output2']
            # XDIM = self.data['b2fgeo']['nx'] + 2
            # YDIM = self.data['b2fgeo']['ny'] + 2
            
            
            RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                        usecols = (3)).reshape((rad_range, pol_range))
            VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                          usecols = (3)).reshape((rad_range, pol_range))
            
            
            "Calculate midplane R for Z=0"
            
            "Calculate weight"
            crup = RadLoc[:, 58]
            crlow = RadLoc[:, 60]
            czup = VertLoc[:, 58]
            czlow = VertLoc[:, 60]
            
            weight_mid = np.zeros(rad_range)
            for x in range(rad_range):
                weight_mid[x] = (czup[x] - mag_axis_z)/ (czup[x] - czlow[x])
            
            mid_choice = np.zeros((rad_range, 4))
            mid_choice[:, 0] = czup
            mid_choice[:, 1] = czlow
            mid_choice[:, 2] = crup
            mid_choice[:, 3] = crlow
            
            
            
            mid_R = np.zeros(rad_range)
            for xa in range(rad_range):
                mid_R[xa] = (1 - weight_mid[xa])*crup[xa] + weight_mid[xa]*crlow[xa]
            
            mid_Z = np.zeros(rad_range)
            for xb in range(rad_range):
                mid_Z[xb] = (1 - weight_mid[xb])*czup[xb] + weight_mid[xb]*czlow[xb]
                
            
            pol_list = [57, 58, 59, 60]
            
            
            psi_solps_mid = np.zeros(rad_range)
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp']
            # psi_solps_cp = psiNinterp_RBS(crloc, czloc)
            for i in range(rad_range):
                psi_mid = psiNinterp_RBS(mid_R[i], mid_Z[i])
                psi_solps_mid[i] = psi_mid
            
            sep = []
            for in_sep in psi_solps_mid:
                if in_sep >= 1:
                    sep.append(list(psi_solps_mid).index(in_sep))
            
                    
            sep_index = sep[0]
            # print(sep)
            
            weight_psi = (psi_solps_mid[sep_index] - 1)/(psi_solps_mid[sep_index] -psi_solps_mid[sep_index -1])
            
            R_sep = (1 - weight_psi)*mid_R[sep_index] + weight_psi*mid_R[sep_index -1]
            R_Rsep = mid_R - R_sep
            
            psi_dsa_dic = fm.dsa_psi_fit(dsa= R_Rsep, psi= psi_solps_mid)           
            self.data['DefaultSettings']['psi_dsa'] = psi_dsa_dic['dsa_psi_fitcoe'][0]
            
            midplane_dic = {'weight': weight_mid, 'mid_choice': mid_choice, 
                            'mid_R': mid_R, 'mid_Z': mid_Z, 
                            'psi_solps_mid': psi_solps_mid, 
                            'weight_psi': weight_psi, 'R_Rsep': R_Rsep}
                
            self.data['midplane_calc'] = midplane_dic
            
            if plotRR:
                plt.figure(figsize=(7,7))
                for in_pol in pol_list:
                    crloc = RadLoc[:, int(in_pol)]
                    czloc = VertLoc[:, int(in_pol)]
                    plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
                plt.plot(mid_R, mid_Z, color = 'r', label= 'R_Rsep')
                plt.xlabel('R: [m]')
                plt.ylabel('Z: [m]')
                plt.title('{} R-Z plot'.format(aa))
            
    
    
    
    
    def calc_pol_angle(self, pol_list):
        
        if self.withshift == False and self.withseries == False:
            
            Attempt = self.data['dircomp']['Attempt']
            DRT = self.data['dirdata']['outputdir']['Output']
            XDIM = self.data['b2fgeo']['nx'] + 2
            YDIM = self.data['b2fgeo']['ny'] + 2
            
            
            RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                        usecols = (3)).reshape((YDIM, XDIM))
            VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                          usecols = (3)).reshape((YDIM,XDIM))
            
            mag_axis_z = self.data['gfile']['g']['zmaxis']
            print(mag_axis_z)
            mag_axis_r = self.data['gfile']['g']['rmaxis'] + self.data['dircomp']['shift_value']
            print(mag_axis_r)
            
            angle_list = []
            
            for pol_index in pol_list:
                rloc = RadLoc[:, int(pol_index)][0]
                zloc = VertLoc[:, int(pol_index)][0]
                
                x = rloc - mag_axis_r
                y = zloc - mag_axis_z
                
                
                angle = np.arctan2(y, x)*180 / np.pi
                if angle <= -90:
                    angle = angle + 360
                angle_list.append(angle)
                
            ln = len(pol_list)
            pol_loc = np.zeros(ln)
            i = 0
            
            for ii in pol_list:
                pol_loc[i] = int(ii)
                i = i + 1
            
            
            plt.figure(figsize=(7,7))
            plt.plot(pol_loc, angle_list, 'o-', color = 'r', label= 'poloidal angle')
            plt.xlabel('poloidal index')
            plt.title('poloidal angle verses poloidal index')
            plt.legend()
            
            self.data['angle'] = angle_list
            
        elif self.withshift == True and self.withseries == False:        
            angle_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                pol_range = int(self.data['b2fgeo'][aa]['nx'] + 2)
                # print('xdim is {}'.format(str(pol_range)))
                rad_range = int(self.data['b2fgeo'][aa]['ny'] + 2)
                # print('ydim is {}'.format(str(rad_range)))
                
                mag_axis_z = self.data['gfile']['g']['zmaxis']
                mag_axis_r = self.data['gfile']['g']['rmaxis'] + self.data['dircomp']['shift_dic'][aa]
                # print(mag_axis_z)
                
                Attempt = self.data['dircomp']['Attempt'][aa]
                DRT = self.data['dirdata']['infolderdir'][aa]['outputdir']['Output']
                
                          
                RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                            usecols = (3)).reshape((rad_range, pol_range))
                VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                              usecols = (3)).reshape((rad_range, pol_range))
                
                angle_list = []
                
                for pol_index in pol_list:
                    rloc = RadLoc[:, int(pol_index)][0]
                    zloc = VertLoc[:, int(pol_index)][0]
                    
                    x = rloc - mag_axis_r
                    y = zloc - mag_axis_z
                    
                    
                    angle = np.arctan2(y, x)*180 / np.pi
                    if angle <= -90:
                        angle = angle + 360
                    angle_list.append(angle)
                    
                    
                                  
                    angle_dic[aa] = angle_list
                    
            ln = len(pol_list)
            pol_loc = np.zeros(ln)
            i = 0
            
            for ii in pol_list:
                pol_loc[i] = int(ii)
                i = i + 1
            
            plt.figure(figsize=(7,7))
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            for ab in self.data['dircomp']['multi_shift']: 
                plt.plot(pol_loc, angle_dic[ab], 'o-', color = color_dic[ab])
                plt.xlabel('poloidal index')
                plt.title('poloidal angle verses poloidal index'.format(self.data['dircomp']['shift_dic'][ab]))
                
                
                
            
                
            self.data['angle'] = angle_dic
            
           
        elif self.withshift == False and self.withseries == True:
        
            pol_range = self.data['b2fgeo']['nx'] + 2
            rad_range = self.data['b2fgeo']['ny'] + 2
            
            
            solps_dsa_dic = {}
            for aa in self.data['dircomp']['Attempt'].keys():
                solps_dsa_dic[aa] = lcm.read_dsa(self.data['dirdata']['simudir'][aa] + '/dsa')
            
            den_list = list(self.data['dircomp']['Attempt'].keys())
            aa = den_list[0]
            
            mag_axis_z = self.data['gfile']['g']['zmaxis']
            print(mag_axis_z)
            
            Attempt = self.data['dircomp']['Attempt'][aa]
            DRT = self.data['dirdata']['outputdir'][aa]['Output']
            
            
            RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                        usecols = (3)).reshape((rad_range, pol_range))
            VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                          usecols = (3)).reshape((rad_range, pol_range))
            
            mag_axis_z = self.data['gfile']['g']['zmaxis']
            print(mag_axis_z)
            mag_axis_r = self.data['gfile']['g']['rmaxis']
            print(mag_axis_r)
            
            angle_list = []
            
            for pol_index in pol_list:
                rloc = RadLoc[:, int(pol_index)][0]
                zloc = VertLoc[:, int(pol_index)][0]
                
                x = rloc - mag_axis_r
                y = zloc - mag_axis_z
                
                
                angle = np.arctan2(y, x)*180 / np.pi
                if angle <= -90:
                    angle = angle + 360
                angle_list.append(angle)
                
            ln = len(pol_list)
            pol_loc = np.zeros(ln)
            i = 0
            
            for ii in pol_list:
                pol_loc[i] = int(ii)
                i = i + 1
            
            
            # plt.figure(figsize=(7,7))
            # plt.plot(pol_loc, angle_list, color = 'r', label= 'angle')
            # plt.xlabel('poloidal index')
            # plt.ylabel('angle')
            # plt.title('angle verses poloidal index')
            
            self.data['angle'] = angle_list
                
            
    
    
    
    
    """
    from stackoverflow
    How to interpolate a 2D curve in Python
    https://stackoverflow.com/questions/52014197/how-to-interpolate-a-2d-curve-in-python


    """

    def RR_sep_calculator(self, cr, cz):
        # Define some points:
        p = len(cr)
        points = np.zeros((p, 2))
        points[:, 0] = cr
        points[:, 1] = cz
        # points = np.array([[0, 1, 8, 2, 2],
        #                    [1, 0, 6, 7, 2]]).T  # a (nbre_points x nbre_dim) array
        
        "Linear length along the line:"
        
        # crzdiff = np.diff(points, axis=0)
        # crzsum = np.sum( np.diff(points, axis=0)**2, axis=1 )
        distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )))
        distance = np.insert(distance, 0, 0)/distance[-1]
        
        arclength = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )))
        arclength = np.insert(arclength, 0, 0)
        
        
        "Interpolation for different methods:"
        
        interpolations_methods = ['slinear', 'quadratic', 'cubic']
        alpha = np.linspace(0, 1, 38)
        
        interpolated_points = {}
        for method in interpolations_methods:
            interpolator =  interp1d(distance, points, kind=method, axis=0)
            interpolated_points[method] = interpolator(alpha)
        
        "Graph:"
        # plt.figure(figsize=(7,7))
        # for method_name, curve in interpolated_points.items():
        #     plt.plot(*curve.T, '-', label=method_name);
        
        # plt.plot(*points.T, 'ok', label='original points');
        # plt.axis('equal'); plt.legend(); plt.xlabel('r'); plt.ylabel('z');

        return arclength, interpolated_points       
    
    
    def calc_sep_index(self, psi, rad_range):
        index_low = []
        index_high = []
        index = []
        # print(psi)
        for y in range(len(psi)):
            if psi[y] <= 1:
                index_low.append(y)
            elif psi[y] > 1:
                index_high.append(y)
            else:
                pass

        
        # print(type(index_high))
        index.append(index_low[-1])
        index.append(index_high[0])
        # print('the following is index_low')
        # print(index_low)
        
        index_dic = {'index_low': index_low, 'index_high': index_high, 
                     'index': index}
        return index_dic
    
    
    def calc_dsa(self, pol_loc):
        
        if self.withshift == False and self.withseries == False:
            # pol_range = int(self.data['b2fgeo']['nx'] + 2)
            rad_range = int(self.data['b2fgeo']['ny'] + 2)
            
            
            Attempt = self.data['dircomp']['Attempt']
            DRT = self.data['dirdata']['outputdir']['Output']
            XDIM = self.data['b2fgeo']['nx'] + 2
            YDIM = self.data['b2fgeo']['ny'] + 2
            dsa = lcm.read_dsa(self.data['dirdata']['simudir'] + '/dsa')
            
            RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                        usecols = (3)).reshape((YDIM, XDIM))
            VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                          usecols = (3)).reshape((YDIM,XDIM))
            
            pol_index = int(pol_loc)
            crloc = RadLoc[:, pol_index]
            czloc = VertLoc[:, pol_index]
            # print(crloc)
            
            arclength, interpfunc_dic = self.RR_sep_calculator(cr = crloc, cz = czloc)
            
            SEP = self.data['DefaultSettings']['SEP']
            
            sep_dist = np.mean([arclength[int(SEP)+ 1], arclength[int(SEP)]])
            # print(sep_dist)
            
            
            RR_sep = arclength - sep_dist
            
            dsa_dic = {'arclength': arclength, 'interpfunc': interpfunc_dic,
                        'dsa_{}_val'.format(pol_index): RR_sep}
            
            self.data['dsa']['dsa_{}'.format(pol_index)] = dsa_dic
            
            
            "Debug tool"
            # comp_dsa = np.zeros((int(rad_range), 2))
            # comp_dsa[:, 0] = dsa
            # comp_dsa[:, 1] = RR_sep
            
            # dsa_dic = {'arclength': arclength, 'interpfunc': interpfunc_dic,
            #            'RadLoc': RadLoc, 'VertLoc': VertLoc,
            #            'dsa_{}_val'.format(pol_index): RR_sep, 'comp_dsa': comp_dsa}
            
            
            "Test for left only result"
            # LLdsa = crLowerLeft - LLsep
            # LRdsa = crLowerRight - LRsep
            # ULdsa = crUpperLeft - ULsep
            # URdsa = crUpperRight - URsep
            # avag_rad = np.zeros(rad_range)
            # for j in range(rad_range):
            #     avag_rad[j] = np.mean([LLdsa[j], ULdsa[j], URdsa[j], LRdsa[j]])
          
            # crzsum, points = fm.RR_sep_calculator(cr = crLowerLeft, cz = czLowerLeft)
            
            # Rad0Cor = np.loadtxt('{}/Rad0Cor{}'.format(DRT2, str(Attempt)),
            #             usecols = (3)).reshape((YDIM, XDIM))
            # Vert0Cor = np.loadtxt('{}/Vert0Cor{}'.format(DRT2, str(Attempt)), 
            #               usecols = (3)).reshape((YDIM,XDIM))
            
            # comp = np.zeros((int(rad_range), 2))
            # comp[:, 0] = Rad0Cor[:, pol_index]
            # comp[:, 1] = crLowerLeft
            
            
        
        elif self.withshift == True and self.withseries == False:
            dsa_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                pol_range = int(self.data['b2fgeo'][aa]['nx'] + 2)
                # print('xdim is {}'.format(str(pol_range)))
                rad_range = int(self.data['b2fgeo'][aa]['ny'] + 2)
                # print('ydim is {}'.format(str(rad_range)))
                
                Attempt = self.data['dircomp']['Attempt'][aa]
                DRT = self.data['dirdata']['infolderdir'][aa]['outputdir']['Output']
            
                RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                            usecols = (3)).reshape((rad_range, pol_range))
                VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                              usecols = (3)).reshape((rad_range, pol_range))
            
                pol_index = int(pol_loc)
                crloc = RadLoc[:, pol_index]
                czloc = VertLoc[:, pol_index]
                # print(crloc)
            
                arclength, interpfunc_dic = self.RR_sep_calculator(cr = crloc, cz = czloc)
                
                SEP = self.data['DefaultSettings']['SEP'][aa]
                
                sep_dist = np.mean([arclength[int(SEP)+ 1], arclength[int(SEP)]])
                # print(sep_dist)
            
            
                RR_sep = arclength - sep_dist
            
                dsa_dic[aa] = {'arclength': arclength, 'interpfunc': interpfunc_dic,
                            'dsa_{}_val'.format(pol_index): RR_sep}
            
            self.data['dsa']['dsa_{}'.format(pol_index)] = dsa_dic
    
    
    
        elif self.withshift == False and self.withseries == True:
            
            den_list = list(self.data['dircomp']['Attempt'].keys())
            aa = den_list[0] 
            
            # pol_range = int(self.data['b2fgeo']['nx'] + 2)
            rad_range = int(self.data['b2fgeo']['ny'] + 2)
            
            
                       
            Attempt = self.data['dircomp']['Attempt'][aa]
            DRT = self.data['dirdata']['outputdir'][aa]['Output']
            XDIM = self.data['b2fgeo']['nx'] + 2
            YDIM = self.data['b2fgeo']['ny'] + 2
            # dsa = lcm.read_dsa(self.data['dirdata']['simudir'] + '/dsa')
            
            RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                        usecols = (3)).reshape((YDIM, XDIM))
            VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                          usecols = (3)).reshape((YDIM,XDIM))
            
            pol_index = int(pol_loc)
            crloc = RadLoc[:, pol_index]
            czloc = VertLoc[:, pol_index]
            # print(crloc)
            
            arclength, interpfunc_dic = self.RR_sep_calculator(cr = crloc, cz = czloc)
            
            SEP = self.data['DefaultSettings']['SEP']
            
            sep_dist = np.mean([arclength[int(SEP)+ 1], arclength[int(SEP)]])
            # print(sep_dist)
            
            
            RR_sep = arclength - sep_dist
            
            dsa_dic = {'arclength': arclength, 'interpfunc': interpfunc_dic,
                        'dsa_{}_val'.format(pol_index): RR_sep}
            
            self.data['dsa']['dsa_{}'.format(pol_index)] = dsa_dic
