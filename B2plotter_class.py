# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 11:54:55 2023

@author: user
"""
import os
import sys
import numpy as np
import xarray as xr
import glob
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors, cm
from scipy.io import loadmat
from D3DPreProcess import PsiNtoR
import equilibrium as eq
from D3DPreProcess import RhotoPsiN
import B2plotter_set as b2s
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.stats import binned_statistic
import matplotlib.pyplot as plt
import B2TransportParser as b2tp


class B2plotter:
    def __init__(self, DEV, withshift, withseries, DefaultSettings):
        
        self.DEV = DEV
        self.withshift = withshift
        self.withseries = withseries
                 
            
        "DefaultSettings"    
        if isinstance(DefaultSettings, dict):
            self.DefaultSettings = DefaultSettings
        else:
            print('parameter has to be a dictionary')
        
        if DefaultSettings is None:
            print('There is no input defaultsettings')
        else:
            self.DefaultSettings = DefaultSettings
         
        keylist = []
        for key, value in self.DefaultSettings.items():
            keylist.append(key)
        
        "Useful data"
        self.data = {'defaultkey':keylist,'dircomp':{}, 'DefaultSettings': {},
                     'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'dsa':{},
                     'gfile':{}, 'gridsettings': {}, 'psi':{}, 
                     'outputdata':{}}

        
    def _reset_object(self):
        # self.Shot=None
        # self.Attempts=None
        # self.Parameter=[]
        # self.PARAM={}
        # self.ExpDict={}
        # self.RadCoords={}
        self.data = {}
        
#-------------load-device-simulation-directory---------------------
    
    def load_mast_dir(self):
        if self.DEV == 'mast':
            if self.withshift == False and self.withseries == False:
                self.data['dircomp'] = b2s.mast_comp_dic()
                mast_basedir, Attempt_dic, shift_value = lmem.mast_base_dir()
                self.data['dirdata'] = mast_basedir
                self.data['dircomp']['Attempt'] = Attempt_dic
                self.data['dircomp']['shift_value'] = shift_value
                
            elif self.withshift == True and self.withseries == False:
                self.data['dircomp'] = b2s.mast_comp_dic_withshift()
                shift_dir, att_dic = lmem.mast_withshift_dir()
                self.data['dirdata'] = shift_dir
                self.data['dircomp']['Attempt'] = att_dic
            
            elif self.withshift == False and self.withseries == True:
                series_flag = self.DefaultSettings['series_flag']
                if series_flag == 'change_den':
                    self.data['dircomp'] = b2s.mast_comp_dir_series()
                    series_dir, att_dic = lmem.mast_series_dir(series_flag= series_flag)
                    self.data['dirdata'] = series_dir
                    self.data['dircomp']['Attempt'] = att_dic
                elif series_flag == 'eireneN':
                    self.data['dircomp'] = b2s.mast_comp_dir_eireneN()
                    series_dir, att_dic = lmem.mast_series_dir(series_flag= series_flag)
                    self.data['dirdata'] = series_dir
                    self.data['dircomp']['Attempt'] = att_dic
            elif self.withshift == True and self.withseries == True:
                print('load_mast_dir is not there yet, to be continue...')      
            else:
                print('There is a bug')

        else:
            print('DEV setting is not mast')
    
 
#-------------------------------------------------------------------       
 
    def loadb2mn_method(self, itername):
        
        'loadb2mn method'
        try:
            if self.withshift == False and self.withseries == False:
                b2mn = lcm.scrape_b2mn(self.data['dirdata']['simudir']
                                      + '/b2mn.dat')
            elif self.withshift == True and self.withseries == False:
                b2mn = lcm.scrape_b2mn(self.data['dirdata']['infolderdir'][itername]['simudir']
                                          + '/b2mn.dat')
            elif self.withshift == False and self.withseries == True:
                b2mn = lcm.scrape_b2mn(self.data['dirdata']['simudir'][itername]
                                      + '/b2mn.dat')
        except:
            print('can not generate b2mn')
        
        return b2mn
        


    def loadgeo_method(self, itername):
        
        
        'loadb2fgmtry method'
        try:
            if self.withshift == False and self.withseries == False:
                geo = lcm.read_b2fgmtry(self.data['dirdata']['simutop'] 
                                       + '/baserun/b2fgmtry')
            elif self.withshift == True and self.withseries == False:
                geo = lcm.read_b2fgmtry(self.data['dirdata']['infolderdir'][itername]['simutop']  
                                           + '/baserun/b2fgmtry')
            elif self.withshift == False and self.withseries == True:
                geo = lcm.read_b2fgmtry(self.data['dirdata']['simutop'] 
                                           + '/baserun/b2fgmtry')
            # print(type(geo))
        except:
            print('can not generate geo')


        g = lcm.loadg(self.data['dirdata']['gdir'][0])
        psiN = (g['psirz'] - g['simag']) / (g['sibry'] - g['simag'])

        dR = g['rdim'] / (g['nw'] - 1)
        dZ = g['zdim'] / (g['nh'] - 1)
        
        gZ = np.zeros(g['nh'])
        for i in range(g['nh']):
            gZ[i] = g['zmid'] - 0.5 * g['zdim'] + i * dZ
                
        if self.withshift == False and self.withseries == False:
            shift = self.data['dircomp']['shift_value']
            # print(shift)
        elif self.withshift == True and self.withseries == False:
            shift = self.data['dircomp']['shift_dic'][itername]
            # print(shift)
        elif self.withshift == False and self.withseries == True:
            shift = self.data['dircomp']['shift_value']
            # print(shift)
            
        
        gR = np.zeros(g['nw'])
        for i in range(g['nw']):
            gR[i] = g['rleft'] + i * dR + float(shift)
            

        psiNinterp_RBS = interpolate.RectBivariateSpline(gR, gZ, np.transpose(psiN))
        psiNinterp_2d = interpolate.interp2d(gR, gZ, psiN, kind = 'cubic')
        psiNinterp_RGI = interpolate.RegularGridInterpolator((gR, gZ), np.transpose(psiN))
        
        interp_dic = {}
        interp_dic['RBS'] = psiNinterp_RBS
        interp_dic['2d'] = psiNinterp_2d
        interp_dic['RGI'] = psiNinterp_RGI
        
        
        gfilesum = {'psiN': psiN, 'dR': dR, 'dZ': dZ, 'gR': gR, 'gZ': gZ,
   'check': 'yeah! shift is {} and series is {}'.format(self.withshift, self.withseries),
                    'interp_dic': interp_dic}
            
        return geo, gfilesum
            
        
    def load_solpsgeo(self):
        
        g = lcm.loadg(self.data['dirdata']['gdir'][0])
        self.data['gfile']['g'] = g
        
        
        if self.withshift == False and self.withseries == False:
            b2mn = self.loadb2mn_method(itername= None)
            self.data['b2mn'] = b2mn
            geo, gfilesum = self.loadgeo_method(itername= None)
            self.data['b2fgeo'] = geo
            self.data['gfile']['gcomp'] = gfilesum
            
        elif self.withshift == True and self.withseries == False:
            b2mn_dic = {}
            geo_dic = {}
            gfilesum_dic = {}
            for shiftname in self.data['dircomp']['multi_shift']:
                b2mn = self.loadb2mn_method(itername= shiftname)
                b2mn_dic[shiftname] = b2mn
                geo, gfilesum = self.loadgeo_method(itername= shiftname)
                geo_dic[shiftname] = geo
                gfilesum_dic[shiftname] = gfilesum
            
            self.data['b2mn'] = b2mn_dic
            self.data['b2fgeo'] = geo_dic      
            self.data['gfile']['gcomp'] = gfilesum_dic
            
            # b2mn_dic = {}
            # geo_dic = {}
            # for aa in self.data['dircomp']['multi_shift']:
            #     try: 
            #         geo = lcm.read_b2fgmtry(self.data['dirdata']['infolderdir'][aa]['simutop']  
            #                                + '/baserun/b2fgmtry')
            #         # print(type(geo))
            #     except:
            #         print('can not generate geo')
                
            #     geo_dic[aa] = geo
                
            #     try:
            #         b2mn = lcm.scrape_b2mn(self.data['dirdata']['infolderdir'][aa]['simudir']
            #                               + '/b2mn.dat')
            #     except:
            #         print('can not generate b2mn')
                
            #     b2mn_dic[aa] = b2mn
                
                
            # self.data['b2mn'] = b2mn_dic
            # self.data['b2fgeo'] = geo_dic
            
            # g = lcm.loadg(self.data['dirdata']['gdir'][0])
            # psiN = (g['psirz'] - g['simag']) / (g['sibry'] - g['simag'])
    
            # dR = g['rdim'] / (g['nw'] - 1)
            # dZ = g['zdim'] / (g['nh'] - 1)
            
            
            # gZ = np.zeros(g['nh'])
            # for i in range(g['nh']):
            #     gZ[i] = g['zmid'] - 0.5 * g['zdim'] + i * dZ
            
            # gR_dic = {}
            # interp_dic = {}
            # for ab in self.data['dircomp']['multi_shift']:
            #     shift = self.data['dircomp']['shift_dic'][ab]
            #     # print(shift)
        
            #     gR = np.zeros(g['nw'])
            #     for i in range(g['nw']):
            #         gR[i] = g['rleft'] + i * dR + float(shift)
                
            #     gR_dic[ab] = gR
                
                
            #     psiNinterp_RBS = interpolate.RectBivariateSpline(gR, gZ, np.transpose(psiN))
            #     # psiNinterp_2d = interpolate.interp2d(gR, gZ, psiN, kind = 'cubic')
            #     # psiNinterp_RGI = interpolate.RegularGridInterpolator((gR, gZ), np.transpose(psiN))
                
            #     interp_dic[ab] = psiNinterp_RBS
            
            
            # gfiledic = {'psiN': psiN, 'dR': dR, 'dZ': dZ, 'gR': gR_dic, 'gZ': gZ,
            #             'check': 'oh_yeah', 'interp_dic': interp_dic}
            
            # self.data['gfile']['g'] = g
            # self.data['gfile']['gcomp'] = gfiledic
        
            # return interp_dic
        
        elif self.withshift == False and self.withseries == True:
            b2mn_dic = {}
            for seriesname in self.data['dircomp']['Attempt'].keys():
                b2mn = self.loadb2mn_method(itername= seriesname)
                b2mn_dic[seriesname] = b2mn
            
            
            self.data['b2mn'] = b2mn
                
            geo, gfilesum = self.loadgeo_method(itername= None) 
            self.data['b2fgeo'] = geo
            self.data['gfile']['gcomp'] = gfilesum
                            
        elif self.withshift == True and self.withseries == True:
            print('load_solpsgeo is not there yet, to be continue...')
        
        else:
            print('There is a bug')


    def calcpsi_method(self, itername):
        
        if self.withshift == False and self.withseries == False:
            pol_range = int(self.data['b2fgeo']['nx'] + 2)
            # print('xdim is {}'.format(pol_range))
            rad_range = int(self.data['b2fgeo']['ny'] + 2)
            # print('ydim is {}'.format(rad_range))
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            # psiNinterp_RGI = self.data['gfile']['gcomp']['interp_dic']['RGI'] 
            # psiNinterp_2d = self.data['gfile']['gcomp']['interp_dic']['2d']
            Attempt = self.data['dircomp']['Attempt']
            DRT = self.data['dirdata']['outputdir']['Output']
            
        elif self.withshift == True and self.withseries == False:
            pol_range = int(self.data['b2fgeo'][itername]['nx'] + 2)
            # print('xdim is {}'.format(str(pol_range)))
            rad_range = int(self.data['b2fgeo'][itername]['ny'] + 2)
            # print('ydim is {}'.format(str(rad_range)))
            psiNinterp_RBS = self.data['gfile']['gcomp'][itername]['interp_dic']['RBS']
            Attempt = self.data['dircomp']['Attempt'][itername]
            DRT = self.data['dirdata']['infolderdir'][itername]['outputdir']['Output']
            
        elif self.withshift == False and self.withseries == True:
            pol_range = int(self.data['b2fgeo']['nx'] + 2)
            # print('xdim is {}'.format(pol_range))
            rad_range = int(self.data['b2fgeo']['ny'] + 2)
            # print('ydim is {}'.format(rad_range))
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            # psiNinterp_RGI = self.data['gfile']['gcomp']['interp_dic']['RGI'] 
            # psiNinterp_2d = self.data['gfile']['gcomp']['interp_dic']['2d']
            Attempt = self.data['dircomp']['Attempt'][itername]
            DRT = self.data['dirdata']['outputdir'][itername]['Output']
            

        psival = np.zeros((rad_range, pol_range))
       
        RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                    usecols = (3)).reshape((rad_range, pol_range))
        VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                      usecols = (3)).reshape((rad_range, pol_range))
        
        coord_dic = {'RadLoc': RadLoc, 'VertLoc': VertLoc}
            
        for pol_loc in range(pol_range):
            for i in range(rad_range):
                # print(i)
                psival[i, pol_loc] = psiNinterp_RBS(RadLoc[i, pol_loc], 
                                                      VertLoc[i, pol_loc])
        return coord_dic, psival
                
    def calcpsi(self):
            
        if self.withshift == False and self.withseries == False:
            
            coord_dic, psival = self.calcpsi_method(itername= None)
            self.data['grid'] = coord_dic
            self.data['psi']['psival'] = psival
        
        elif self.withshift == True and self.withseries == False:
            pol_range_dic = {}
            rad_range_dic = {}
            SEP_dic = {}
            dsa_dic = {}
            psival_dic = {}
            solps_dsa_dic = {}
            coord_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                
                geo = self.data['b2fgeo'][aa]
                pol_range = int(self.data['b2fgeo'][aa]['nx'] + 2)
                # print('xdim is {}'.format(str(pol_range)))
                rad_range = int(self.data['b2fgeo'][aa]['ny'] + 2)
                # print('ydim is {}'.format(str(rad_range)))
                pol_range_dic[aa] = pol_range
                rad_range_dic[aa] = rad_range
            
                solps_dsa_dic[aa] = lcm.read_dsa(self.data['dirdata']['infolderdir'][aa]['simudir'] + '/dsa')
    
                psiNinterp_RBS = self.data['gfile']['gcomp'][aa]['interp_dic']['RBS']

                psival = np.zeros((rad_range, pol_range))
                
                Attempt = self.data['dircomp']['Attempt'][aa]
                DRT = self.data['dirdata']['infolderdir'][aa]['outputdir']['Output']
                
                
                RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                            usecols = (3)).reshape((rad_range, pol_range))
                VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                              usecols = (3)).reshape((rad_range, pol_range))
                
                coord_dic = {'RadLoc': RadLoc, 'VertLoc': VertLoc}
                
                self.data['grid'][aa] = coord_dic
                
                for pol_loc in range(pol_range):
                    for i in range(rad_range):
                        # print(i)
                        psival[i, pol_loc] = psiNinterp_RBS(RadLoc[i, pol_loc], 
                                                              VertLoc[i, pol_loc])
                      
                self.data['psi'][aa] = psival
            
        elif self.withshift == False and self.withseries == True:
            # print('we are working on calcpsi for series case')
            series_rep = list(self.data['dircomp']['Attempt'].keys())[0]
            self.calcpsi_method(itername= series_rep)
    
    def plot_seperatrix(self):
        psi_1d = self.data['psi']['psival'][0, :]
        # self.data['psi']['psi1d'] = psi_1d
        
        pol_range = int(self.data['b2fgeo']['nx'] + 2)
        rad_range = int(self.data['b2fgeo']['ny'] + 2)
        
        
        index_low = []
        index_high = []
        index = np.zeros(2)
        for y in range(rad_range):
            if psi_1d[y] <= 1:
                index_low.append(y)
            if psi_1d[y] >= 1:
                index_high.append(y)
    
        
        index[0] = index_low[-1]
        index[1] = index_high[0]
        
        index_dic = {'index_low': index_low, 'index_high': index_high, 
                     'index': index}
        self.data['index'] = index_dic
        
        
        
    
    def calcpsi_1D(self, pol_loc):
        
        if self.withshift == False and self.withseries == False:
            geo = self.data['b2fgeo']
            pol_range = int(self.data['b2fgeo']['nx'] + 2)
            # print('xdim is {}'.format(str(pol_range)))
            rad_range = int(self.data['b2fgeo']['ny'] + 2)
            # print('ydim is {}'.format(str(rad_range)))
            
            self.data['DefaultSettings']['XDIM'] = pol_range
            self.data['DefaultSettings']['YDIM'] = rad_range
            
    
            psiNinterp_RGI = self.data['gfile']['gcomp']['interp_dic']['RGI'] 
            psiNinterp_2d = self.data['gfile']['gcomp']['interp_dic']['2d']
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            # print(type(psiNinterp_RBS))
            psival = np.zeros((pol_range, rad_range))
            
            pol_index = int(pol_loc)
            
            crLowerLeft = geo['crx'][pol_index,:,0]
            crLowerRight = geo['crx'][pol_index,:,1]
            crUpperLeft = geo['crx'][pol_index,:,2]
            crUpperRight = geo['crx'][pol_index,:,3]
            czLowerLeft = geo['cry'][pol_index,:,0]
            czLowerRight = geo['cry'][pol_index,:,1]
            czUpperLeft = geo['cry'][pol_index,:,2]
            czUpperRight = geo['cry'][pol_index,:,3]
                
            # LLsep = np.mean([crLowerLeft[int(SEP)-1], crLowerLeft[int(SEP)-2]])
            # LRsep = np.mean([crLowerRight[int(SEP)-1], crLowerRight[int(SEP)-2]])
            # ULsep = np.mean([crUpperLeft[int(SEP)-1], crUpperLeft[int(SEP)-2]])
            # URsep = np.mean([crUpperRight[int(SEP)-1], crUpperRight[int(SEP)-2]])
            
            # weight = dsa[int(SEP)-1]/ (dsa[int(SEP)-1] - dsa[int(SEP)-2])
            # print(weight)
            
            
            "dsa check"
            # kk = np.zeros((int(rad_range),3))
            # kk[:, 0] = dsa
            # kk[:, 1] = avag_rad
            
            # del_dsa = np.zeros(rad_range)
            # for ia in range(rad_range):
            #     del_dsa[ia] = avag_rad[ia] - dsa[ia]
            
            # kk[:, 2] = del_dsa
            # self.data['psi']['dsa_{}_val'.format(pol_loc)] = avag_rad
            
            
            crz_LL = np.stack([crLowerLeft.ravel(), czLowerLeft.ravel()], -1)  # shape (N, 2) in 2d
            crz_UL = np.stack([crUpperLeft.ravel(), czUpperLeft.ravel()], -1)
            crz_LR = np.stack([crLowerRight.ravel(), czLowerRight.ravel()], -1)  # shape (N, 2) in 2d
            crz_UR = np.stack([crUpperRight.ravel(), czUpperRight.ravel()], -1)
            psi_solps_LL = psiNinterp_RGI(crz_LL)
            psi_solps_UL = psiNinterp_RGI(crz_UL)
            psi_solps_LR = psiNinterp_RGI(crz_LR)
            psi_solps_UR = psiNinterp_RGI(crz_UR)
            
            
            psi_solps = np.zeros(rad_range)
            for i in range(rad_range):
                psi_solps[i] = np.mean([psi_solps_LL[i], psi_solps_UL[i], 
                                        psi_solps_LR[i], psi_solps_UR[i]])
            
            psi_solps_2d = np.zeros(rad_range)
            for i in range(rad_range):
                psi_LL_2d = psiNinterp_2d(crLowerLeft[i], czLowerLeft[i])
                psi_UL_2d = psiNinterp_2d(crUpperLeft[i], czUpperLeft[i])
                psi_LR_2d = psiNinterp_2d(crLowerRight[i], czLowerRight[i])
                psi_UR_2d = psiNinterp_2d(crUpperRight[i], czUpperRight[i])
                psi_solps_2d[i] = np.mean([psi_LL_2d, psi_UL_2d, 
                                           psi_LR_2d, psi_UR_2d])
             
            
            psi_solps_RBS = np.zeros(rad_range)
            for i in range(rad_range):
                psi_LL_RBS = psiNinterp_RBS(crLowerLeft[i], czLowerLeft[i])
                psi_UL_RBS = psiNinterp_RBS(crUpperLeft[i], czUpperLeft[i])
                psi_LR_RBS = psiNinterp_RBS(crLowerRight[i], czLowerRight[i])
                psi_UR_RBS = psiNinterp_RBS(crUpperRight[i], czUpperRight[i])
                psi_solps_RBS[i] = np.mean([psi_LL_RBS, psi_UL_RBS, 
                                            psi_LR_RBS, psi_UR_RBS])
            
            "Only work for original case"
            # GF = eq.equilibrium(gfile= self.data['dirdata']['gdir'][0])
            # print(type(GF.psiN._func))
            # psi_solps_GF = np.zeros(rad_range)
            # for i in range(rad_range):
            #     psi_LL_GF = GF.psiN(crLowerLeft[i], czLowerLeft[i])
            #     psi_UL_GF = GF.psiN(crUpperLeft[i], czUpperLeft[i])
            #     psi_LR_GF = GF.psiN(crLowerRight[i], czLowerRight[i])
            #     psi_UR_GF = GF.psiN(crUpperRight[i], czUpperRight[i])
            #     psi_solps_GF[i] = np.mean([psi_LL_GF, psi_UL_GF, 
            #                                 psi_LR_GF, psi_UR_GF])
            
            
            Attempt = self.data['dircomp']['Attempt']
            DRT = self.data['dirdata']['outputdir']['Output']
            XDIM = self.data['b2fgeo']['nx'] + 2
            YDIM = self.data['b2fgeo']['ny'] + 2
            
            RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                        usecols = (3)).reshape((YDIM, XDIM))
            VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                          usecols = (3)).reshape((YDIM,XDIM))
                 
            crloc = RadLoc[:, pol_index]
            czloc = VertLoc[:, pol_index]
            
            
            psi_solps_cp = np.zeros(rad_range)
            # psi_solps_cp = psiNinterp_RBS(crloc, czloc)
            for i in range(rad_range):
                psi_CP = psiNinterp_RBS(crloc[i], czloc[i])
                psi_solps_cp[i] = psi_CP
            
            psival = np.zeros((int(rad_range), 4))
            psival[:, 0] = psi_solps
            psival[:, 1] = psi_solps_2d
            psival[:, 2] = psi_solps_RBS
            psival[:, 3] = psi_solps_cp
            
            index_dic = self.calc_sep_index(psi = psi_solps_RBS, 
                                            rad_range = rad_range)
            
            self.data['DefaultSettings']['SEP'] = index_dic['index'][0]             
            self.data['psi']['psi_{}_val'.format(pol_loc)] = psival
        
        elif self.withshift == True and self.withseries == False:
            pol_range_dic = {}
            rad_range_dic = {}
            SEP_dic = {}
            dsa_dic = {}
            psival_dic = {}
            solps_dsa_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                geo = self.data['b2fgeo'][aa]
                pol_range = int(self.data['b2fgeo'][aa]['nx'] + 2)
                # print('xdim is {}'.format(str(pol_range)))
                rad_range = int(self.data['b2fgeo'][aa]['ny'] + 2)
                # print('ydim is {}'.format(str(rad_range)))
                pol_range_dic[aa] = pol_range
                rad_range_dic[aa] = rad_range
          
                # if rad_range % 2 == 0:
                #     SEP = rad_range/ 2 + 1
                # else:
                #     SEP = round(rad_range/ 2)
                # SEP_dic[aa] = SEP       
                solps_dsa_dic[aa] = lcm.read_dsa(self.data['dirdata']['infolderdir'][aa]['simudir'] + '/dsa')

                psiNinterp_RBS = self.data['gfile']['gcomp'][aa]['interp_dic']['RBS']
                # print(type(psiNinterp_RBS))
                
                pol_index = int(pol_loc)
                
                crLowerLeft = geo['crx'][pol_index,:,0]
                crLowerRight = geo['crx'][pol_index,:,1]
                crUpperLeft = geo['crx'][pol_index,:,2]
                crUpperRight = geo['crx'][pol_index,:,3]
                czLowerLeft = geo['cry'][pol_index,:,0]
                czLowerRight = geo['cry'][pol_index,:,1]
                czUpperLeft = geo['cry'][pol_index,:,2]
                czUpperRight = geo['cry'][pol_index,:,3]
                                  
                # weight = dsa[int(SEP)-1]/ (dsa[int(SEP)-1] - dsa[int(SEP)-2])
                # print(weight)
                
                
                psi_solps_RBS = np.zeros(rad_range)
                for i in range(rad_range):
                    psi_LL_RBS = psiNinterp_RBS(crLowerLeft[i], czLowerLeft[i])
                    psi_UL_RBS = psiNinterp_RBS(crUpperLeft[i], czUpperLeft[i])
                    psi_LR_RBS = psiNinterp_RBS(crLowerRight[i], czLowerRight[i])
                    psi_UR_RBS = psiNinterp_RBS(crUpperRight[i], czUpperRight[i])
                    psi_solps_RBS[i] = np.mean([psi_LL_RBS, psi_UL_RBS, 
                                                psi_LR_RBS, psi_UR_RBS])
                
                
                psival = np.zeros(int(rad_range))
                psival = psi_solps_RBS
                
                psival_dic[aa] = psival
                
                
                
                index_dic = self.calc_sep_index(psi = psival, 
                                                rad_range = rad_range)
                SEP_dic[aa] = index_dic['index'][0]
                
                
            self.data['DefaultSettings']['XDIM'] = pol_range_dic
            self.data['DefaultSettings']['YDIM'] = rad_range_dic
            self.data['DefaultSettings']['SEP'] = SEP_dic
            self.data['dsa']['dsa_{}'.format(pol_loc)] = dsa_dic
            self.data['psi']['psi_{}_val'.format(pol_loc)] = psival_dic 
            
        elif self.withshift == False and self.withseries == True:
            geo = self.data['b2fgeo']
            pol_range = int(self.data['b2fgeo']['nx'] + 2)
            # print('xdim is {}'.format(str(pol_range)))
            rad_range = int(self.data['b2fgeo']['ny'] + 2)
            # print('ydim is {}'.format(str(rad_range)))
            # for aa in self.data['dircomp']['Attemptkey']:
            
            solps_dsa_dic = {}
            for aa in self.data['dircomp']['Attempt'].keys():
                solps_dsa_dic[aa] = lcm.read_dsa(self.data['dirdata']['simudir'][aa] + '/dsa')

            psiNinterp_RBS = self.data['gfile']['gcomp']['interp']
            # print(type(psiNinterp_RBS))
            
            pol_index = int(pol_loc)
            
            crLowerLeft = geo['crx'][pol_index,:,0]
            crLowerRight = geo['crx'][pol_index,:,1]
            crUpperLeft = geo['crx'][pol_index,:,2]
            crUpperRight = geo['crx'][pol_index,:,3]
            czLowerLeft = geo['cry'][pol_index,:,0]
            czLowerRight = geo['cry'][pol_index,:,1]
            czUpperLeft = geo['cry'][pol_index,:,2]
            czUpperRight = geo['cry'][pol_index,:,3]
                
            
            psi_solps_RBS = np.zeros(rad_range)
            for i in range(rad_range):
                psi_LL_RBS = psiNinterp_RBS(crLowerLeft[i], czLowerLeft[i])
                psi_UL_RBS = psiNinterp_RBS(crUpperLeft[i], czUpperLeft[i])
                psi_LR_RBS = psiNinterp_RBS(crLowerRight[i], czLowerRight[i])
                psi_UR_RBS = psiNinterp_RBS(crUpperRight[i], czUpperRight[i])
                psi_solps_RBS[i] = np.mean([psi_LL_RBS, psi_UL_RBS, 
                                            psi_LR_RBS, psi_UR_RBS])
            
            
            # psival = np.zeros(int(rad_range))
            psival = psi_solps_RBS
        
        
            index_dic = self.calc_sep_index(psi = psival, 
                                            rad_range = rad_range)
            SEP = index_dic['index'][0]
        
            
            self.data['DefaultSettings']['XDIM'] = pol_range
            self.data['DefaultSettings']['YDIM'] = rad_range
            self.data['DefaultSettings']['SEP'] = SEP
            self.data['psi']['psi_{}_val'.format(pol_loc)] = psival
        
        elif self.withshift == True and self.withseries == True:
            print('calcpsi_1D is not there yet, to be continue...')
            
        else:
            print('There is a bug')


    def RR_cal_dsa(self, pol_loc):
        geo = self.data['b2fgeo']
        pol_range = int(self.data['b2fgeo']['nx'] + 2)
        # print('xdim is {}'.format(str(pol_range)))
        rad_range = int(self.data['b2fgeo']['ny'] + 2)
        
        pol_index = int(pol_loc)
        
        crLowerLeft = geo['crx'][pol_index,:,0]
        
        Attempt = self.data['dircomp']['Attempt']
        DRT = self.data['dirdata']['infolderdir']['outputdir']['Output']
        DRT2 = self.data['dirdata']['infolderdir']['outputdir']['Output2']
        XDIM = self.data['b2fgeo']['nx'] + 2
        YDIM = self.data['b2fgeo']['ny'] + 2
        
        
        
        Rad0Cor = np.loadtxt('{}/Rad0Cor{}'.format(DRT2, str(Attempt)),
                    usecols = (3)).reshape((YDIM, XDIM))
        Vert0Cor = np.loadtxt('{}/Vert0Cor{}'.format(DRT2, str(Attempt)), 
                      usecols = (3)).reshape((YDIM,XDIM))
        
        comp = np.zeros((int(rad_range), 2))
        comp[:, 0] = Rad0Cor[:, pol_index]
        comp[:, 1] = crLowerLeft
        
        RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                    usecols = (3)).reshape((YDIM, XDIM))
        VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                      usecols = (3)).reshape((YDIM,XDIM))
        
        crloc = RadLoc[:, pol_index]
        czloc = VertLoc[:, pol_index]

                
    def calc_flux_expansion(self, pol_loc, ped_index, iter_index):
        
        if self.withshift == False and self.withseries == False:
            
            arcR = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]
            RR_sep = self.data['midplane_calc']['R_Rsep']
            
            arcR_inv = list(reversed(arcR))
            RRsep_inv = list(reversed(RR_sep))
                       
            arcR_cut = []
            RRsep_cut = []
            
            for p_in in ped_index:
                arcR_cut.append(arcR_inv[p_in])
                RRsep_cut.append(RRsep_inv[p_in])
                
                       
            flux_fit_dic = fm.flux_expand_fit(RRsep = RRsep_cut, arclength = arcR_cut)
            
            flux_expand = flux_fit_dic['flux_fitcoe'][0]
            
            return flux_expand
        
        elif self.withshift == True and self.withseries == False:
            
            arcR = self.data['dsa']['dsa_{}'.format(pol_loc)][iter_index]['dsa_{}_val'.format(pol_loc)]
            RR_sep = self.data['midplane_calc'][iter_index]['R_Rsep']
            
            arcR_inv = list(reversed(arcR))
            RRsep_inv = list(reversed(RR_sep))
            
            
            arcR_cut = []
            RRsep_cut = []
            
            for p_in in ped_index:
                arcR_cut.append(arcR_inv[p_in])
                RRsep_cut.append(RRsep_inv[p_in])
                
                       
            flux_fit_dic = fm.flux_expand_fit(RRsep = RRsep_cut, arclength = arcR_cut)
            
            flux_expand = flux_fit_dic['flux_fitcoe'][0]
            
            return flux_expand
        
        elif self.withshift == False and self.withseries == True:
            
            arcR = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]
            RR_sep = self.data['midplane_calc']['R_Rsep']
            
            arcR_inv = list(reversed(arcR))
            RRsep_inv = list(reversed(RR_sep))
            
            arcR_cut = []
            RRsep_cut = []
            
            for p_in in ped_index:
                arcR_cut.append(arcR_inv[p_in])
                RRsep_cut.append(RRsep_inv[p_in])
                
                       
            flux_fit_dic = fm.flux_expand_fit(RRsep = RRsep_cut, arclength = arcR_cut)
            
            flux_expand = flux_fit_dic['flux_fitcoe'][0]
            
            return flux_expand
        
        elif self.withshift == True and self.withseries == True:
            print('calc_flux_expansion is not there yet, to be continue...')
            
        else:
            print('There is a bug')
    
    

'Way to generate align transport coefficient'
    
"""
m = len(yd)

one_trans = b2tp.InputfileParser(one_list[0], plot= False)
ond = one_trans['1'].T
onki = one_trans['3'].T
onke = one_trans['4'].T
onx= ond[:,0]  #the coordinate here is R-R_sep
fd = ond[:,1]
fki = onki[:,1]
fke = onke[:,1]

d_func = interpolate.interp1d(x, yd, fill_value = 'extrapolate')
ond[:,1] = d_func(onx)
ki_func = interpolate.interp1d(x, yki, fill_value = 'extrapolate')
onki[:,1] = ki_func(onx)
ke_func = interpolate.interp1d(x, yke, fill_value = 'extrapolate')
onke[:,1] = ke_func(onx)


b = b2tp.Generate(cod, CoeffID=1, SpeciesID=2, M=[1])
c = b2tp.WriteInputfile(file='b2.transport.inputfile_align_{}_{}'.format(shift_a, n), points= one_trans ,M_1 = True, M=[1])
"""
            
            


        
            
            
        
    
