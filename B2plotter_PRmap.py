# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 17:27:01 2023

@author: user
"""

from B2plotter_class import load_data
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


class mapping_RRsep_to_psi(load_data):
    def __init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters):
        load_data.__init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters)
    
        
    def calc_RRsep(self):
        
        if self.withshift == False and self.withseries == False:
        
            pol_range = int(self.data['b2fgeo']['nx'] + 2)
            rad_range = int(self.data['b2fgeo']['ny'] + 2)
            
            
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
            
            crup = RadLoc[:, 53]
            crlow = RadLoc[:, 55]
            czup = VertLoc[:, 53]
            czlow = VertLoc[:, 55]
            
            pol_list = [52, 53, 54, 55, 56, 57]
            

            weight_mid = np.zeros(rad_range)
            for x in range(rad_range):
                weight_mid[x] = czup[x]/ (czup[x] - czlow[x])
            
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
                
            pol_list = [52, 53, 54, 55, 56, 57]
            
            # plt.figure(figsize=(7,7))
            # for in_pol in pol_list:
            #     crloc = RadLoc[:, int(in_pol)]
            #     czloc = VertLoc[:, int(in_pol)]
            #     plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
            # plt.plot(mid_R, mid_Z, color = 'r', label= 'R_Rsep')
            # plt.xlabel('R: [m]')
            # plt.ylabel('Z: [m]')
            
                       
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
            print(sep)
            
            weight_psi = (psi_solps_mid[sep_index] - 1)/(psi_solps_mid[sep_index] -psi_solps_mid[sep_index -1])
            
            R_sep = (1 - weight_psi)*mid_R[sep_index] + weight_psi*mid_R[sep_index -1]
            R_Rsep = mid_R - R_sep
            
            midplane_dic = {'weight': weight_mid, 'mid_choice': mid_choice, 
                            'mid_R': mid_R, 'mid_Z': mid_Z, 
                            'psi_solps_mid': psi_solps_mid, 
                            'weight_psi': weight_psi, 'R_Rsep': R_Rsep}
                
            self.data['midplane_calc'] = midplane_dic
            
               
        
        elif self.withshift == True and self.withseries == False:        
            pol_range_dic = {}
            rad_range_dic = {}
            SEP_dic = {}
            midplane_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                pol_range = int(self.data['b2fgeo'][aa]['nx'] + 2)
                # print('xdim is {}'.format(str(pol_range)))
                rad_range = int(self.data['b2fgeo'][aa]['ny'] + 2)
                # print('ydim is {}'.format(str(rad_range)))
                pol_range_dic[aa] = pol_range
                rad_range_dic[aa] = rad_range
                
                
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
                    crup = RadLoc[:, 53]
                    crlow = RadLoc[:, 55]
                    czup = VertLoc[:, 53]
                    czlow = VertLoc[:, 55]
                else:
                    crup = RadLoc[:, 53]
                    crlow = RadLoc[:, 55]
                    czup = VertLoc[:, 53]
                    czlow = VertLoc[:, 55]
                
                
                weight_mid = np.zeros(rad_range)
                for x in range(rad_range):
                    weight_mid[x] = czup[x]/ (czup[x] - czlow[x])
                
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
                    
                
                pol_list = [52, 53, 54, 55, 56, 57]
                
                # plt.figure(figsize=(7,7))
                # for in_pol in pol_list:
                #     crloc = RadLoc[:, int(in_pol)]
                #     czloc = VertLoc[:, int(in_pol)]
                #     plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
                # plt.plot(mid_R, mid_Z, color = 'r', label= 'R_Rsep')
                # plt.xlabel('R: [m]')
                # plt.ylabel('Z: [m]')
                # plt.title('{} R-Z plot'.format(aa))
                
                
                
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
                    
            
            self.data['midplane_calc'] = midplane_dic
            self.data['DefaultSettings']['XDIM'] = pol_range_dic
            self.data['DefaultSettings']['YDIM'] = rad_range_dic
            self.data['DefaultSettings']['SEP'] = SEP_dic
            
        if self.withshift == False and self.withseries == True:
        
            pol_range = self.data['b2fgeo']['nx'] + 2
            rad_range = self.data['b2fgeo']['ny'] + 2
            
            
            solps_dsa_dic = {}
            for aa in self.data['dircomp']['Attempt'].keys():
                solps_dsa_dic[aa] = lcm.read_dsa(self.data['dirdata']['simudir'][aa] + '/dsa')
            
            den_list = list(self.data['dircomp']['Attempt'].keys())
            aa = den_list[0]
            
            
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
            crup = RadLoc[:, 53]
            crlow = RadLoc[:, 55]
            czup = VertLoc[:, 53]
            czlow = VertLoc[:, 55]
            
            weight_mid = np.zeros(rad_range)
            for x in range(rad_range):
                weight_mid[x] = czup[x]/ (czup[x] - czlow[x])
            
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
                
            
            pol_list = [52, 53, 54, 55, 56, 57]
            
            # plt.figure(figsize=(7,7))
            # for in_pol in pol_list:
            #     crloc = RadLoc[:, int(in_pol)]
            #     czloc = VertLoc[:, int(in_pol)]
            #     plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
            # plt.plot(mid_R, mid_Z, color = 'r', label= 'R_Rsep')
            # plt.xlabel('R: [m]')
            # plt.ylabel('Z: [m]')
            # plt.title('{} R-Z plot'.format(aa))
            
            
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
            
            midplane_dic = {'weight': weight_mid, 'mid_choice': mid_choice, 
                            'mid_R': mid_R, 'mid_Z': mid_Z, 
                            'psi_solps_mid': psi_solps_mid, 
                            'weight_psi': weight_psi, 'R_Rsep': R_Rsep}
                
            self.data['midplane_calc'] = midplane_dic
            
    
    
    """
    from stackoverflow
    How to interpolate a 2D curve in Python
    https://stackoverflow.com/questions/52014197/how-to-interpolate-a-2d-curve-in-python


    """

    def RR_sep_calculator(cr, cz):
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
    
    
    
    
    def calc_dsa(self, pol_list):
        
        if self.withshift == False and self.withseries == False:
        
            dsa = lcm.read_dsa(self.data['dirdata']['simudir'] + '/dsa')
            
            
            LLdsa = crLowerLeft - LLsep
            LRdsa = crLowerRight - LRsep
            ULdsa = crUpperLeft - ULsep
            URdsa = crUpperRight - URsep
            avag_rad = np.zeros(rad_range)
            for j in range(rad_range):
                avag_rad[j] = np.mean([LLdsa[j], ULdsa[j], URdsa[j], LRdsa[j]])
          
            # crzsum, points = fm.RR_sep_calculator(cr = crLowerLeft, cz = czLowerLeft)
            
            # Rad0Cor = np.loadtxt('{}/Rad0Cor{}'.format(DRT2, str(Attempt)),
            #             usecols = (3)).reshape((YDIM, XDIM))
            # Vert0Cor = np.loadtxt('{}/Vert0Cor{}'.format(DRT2, str(Attempt)), 
            #               usecols = (3)).reshape((YDIM,XDIM))
            
            # comp = np.zeros((int(rad_range), 2))
            # comp[:, 0] = Rad0Cor[:, pol_index]
            # comp[:, 1] = crLowerLeft
            
            RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                        usecols = (3)).reshape((YDIM, XDIM))
            VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                          usecols = (3)).reshape((YDIM,XDIM))
            
            crloc = RadLoc[:, pol_index]
            czloc = VertLoc[:, pol_index]
            
            
            arclength, interpfunc_dic = fm.RR_sep_calculator(cr = crloc, cz = czloc)
            
            
            sep_dist = np.mean([arclength[int(SEP)-1], arclength[int(SEP)-2]])
            # print(sep_dist)
            
            
            RR_sep = arclength - sep_dist
            
            
            comp_dsa = np.zeros((int(rad_range), 2))
            comp_dsa[:, 0] = dsa
            comp_dsa[:, 1] = RR_sep
            
            dsa_dic = {'arclength': arclength, 'interpfunc': interpfunc_dic,
                       'RadLoc': RadLoc, 'VertLoc': VertLoc,
                       'dsa_{}_val'.format(pol_loc): RR_sep, 'comp_dsa': comp_dsa}
            self.data['dsa']['dsa_{}'.format(pol_loc)] = dsa_dic
        
        elif self.withshift == True and self.withseries == False:
        
            crloc = np.mean([crLowerLeft, crLowerRight, 
                               crUpperLeft, crUpperRight], axis= 0)
            
            czloc = np.mean([czLowerLeft, czLowerRight, 
                               czUpperLeft, czUpperRight], axis= 0)
            
            # comp_cr = np.zeros((int(rad_range), 5))
            # comp_cr[:, 0] = crloc
            # comp_cr[:, 1] = crLowerLeft
            # comp_cr[:, 2] = crLowerRight
            # comp_cr[:, 3] = crUpperLeft
            # comp_cr[:, 4] = crUpperRight
            
            # comp_cz = np.zeros((int(rad_range), 5))
            # comp_cz[:, 0] = czloc
            # comp_cz[:, 1] = czLowerLeft
            # comp_cz[:, 2] = czLowerRight
            # comp_cz[:, 3] = czUpperLeft
            # comp_cz[:, 4] = czUpperRight
            
            
            arclength, interpfunc_dic = fm.RR_sep_calculator(cr = crloc, cz = czloc)
               
            sep_dist = np.mean([arclength[int(SEP)-1], arclength[int(SEP)-2]])
            # print(sep_dist)
            
            
            RR_sep = arclength - sep_dist
            
            comp = np.zeros((int(rad_range), 2))
            comp[:, 0] = solps_dsa_dic[aa]
            comp[:, 1] = RR_sep
            
            
            dsa_dic[aa] = {'arclength': arclength, 'interpfunc': interpfunc_dic,
                           'comp': comp,
                           'dsa_{}_val'.format(pol_loc): RR_sep}

    
    
    
        elif self.withshift == False and self.withseries == True:
            crloc = np.mean([crLowerLeft, crLowerRight, 
                               crUpperLeft, crUpperRight], axis= 0)
            
            czloc = np.mean([czLowerLeft, czLowerRight, 
                               czUpperLeft, czUpperRight], axis= 0)
            
            # comp_cr = np.zeros((int(rad_range), 5))
            # comp_cr[:, 0] = crloc
            # comp_cr[:, 1] = crLowerLeft
            # comp_cr[:, 2] = crLowerRight
            # comp_cr[:, 3] = crUpperLeft
            # comp_cr[:, 4] = crUpperRight
            
            # comp_cz = np.zeros((int(rad_range), 5))
            # comp_cz[:, 0] = czloc
            # comp_cz[:, 1] = czLowerLeft
            # comp_cz[:, 2] = czLowerRight
            # comp_cz[:, 3] = czUpperLeft
            # comp_cz[:, 4] = czUpperRight
            
            arclength, interpfunc_dic = fm.RR_sep_calculator(cr = crloc, cz = czloc)
            
            
            sep_dist = np.mean([arclength[int(SEP)-1], arclength[int(SEP)-2]])
            # print(sep_dist)
            
            
            RR_sep = arclength - sep_dist
            
            
            # lk = len(self.data['dircomp']['Attempt'].keys())
            # comp_dsa = np.zeros((int(rad_range), lk + 1))
            # ii = 0
            # for i in self.data['dircomp']['Attempt'].keys():
            #     comp_dsa[:, ii] = solps_dsa_dic[i]
            #     ii = ii + 1
            # comp_dsa[:, ii] = RR_sep
            
            # dsa_dic = {'arclength': arclength, 'interpfunc': interpfunc_dic,
            #         'comp_cr': comp_cr, 'comp_cz': comp_cz, 'comp_dsa': comp_dsa,
            #             'dsa_{}_val'.format(pol_loc): RR_sep}
            
            dsa_dic = {'arclength': arclength, 'interpfunc': interpfunc_dic,
                        'dsa_{}_val'.format(pol_loc): RR_sep}


     
            
    def Opacity_study_poloidal_plot(self, pol_list, x_choice):
        self.data['poloidal_index'] = pol_list