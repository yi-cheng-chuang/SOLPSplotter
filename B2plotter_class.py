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
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from scipy.io import loadmat
from scipy import interpolate
from D3DPreProcess import PsiNtoR
from D3DPreProcess import RhotoPsiN
import B2plotter_set as b2s
import mastload_data_method as mdm
import coord_sutils as cs 
from scipy import interpolate


class B2plotter:
    def __init__(self, DEV, DefaultSettings):
        
        self.DEV = DEV   

                 
            
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
                     'grid':{}, 'dirdata':{}, 'ExpDict': {}, 'RadCoords':{},
                     'gfile':{}, 'gridsettings': {}, 'psi':{}, 'outputdata':{}}

        
    def _reset_object(self):
        self.Shot=None
        self.Attempts=None
        self.Parameter=[]
        self.PARAM={}
        self.ExpDict={}
        self.RadCoords={}
        
#-------------load-device-simulation-directory---------------------
    
    def load_mast_dir(self):
        if self.DEV == 'mast':
            mast_basedir, Attempt_dic = mdm.mast_base_dir()
            self.data['dirdata'] = mast_basedir
        else:
            print('DEV setting is not mast')
        
        self.data['dircomp'] = b2s.mast_comp_dic()
        self.data['dircomp']['Attempt'] = Attempt_dic
 
#-------------------------------------------------------------------       
 
    def load_solpsgeo(self):
        try: 
            geo = cs.read_b2fgmtry(self.data['dirdata']['tbase'] 
                                   + '/baserun/b2fgmtry')
            # print(type(geo))
        except:
            print('can not generate geo')
        
        try:
            b2mn = cs.scrape_b2mn(self.data['dirdata']['adir']['simudir']
                                  + '/b2mn.dat')
        except:
            print('can not generate b2mn')
        
        
        self.data['b2mn'] = b2mn
        self.data['b2fgeo'] = geo
        
        
        g = cs.loadg(self.data['dirdata']['gdir'][0])
        psiN = (g['psirz'] - g['simag']) / (g['sibry'] - g['simag'])

        dR = g['rdim'] / (g['nw'] - 1)
        dZ = g['zdim'] / (g['nh'] - 1)
        
        
        ast = self.data['dircomp']['a_shift']
        shift = self.data['dircomp']['shiftdic'][ast]
        # print(shift)

        gR = []
        for i in range(g['nw']):
            gR.append(g['rleft'] + i * dR + float(shift))

        gZ = []
        for i in range(g['nh']):
            gZ.append(g['zmid'] - 0.5 * g['zdim'] + i * dZ)
        
        psiNinterp = interpolate.RectBivariateSpline(gR, gZ, psiN)
            
        gfiledic = {'psiN': psiN, 'dR': dR, 'dZ': dZ, 'gR': gR, 'gZ': gZ,
                    'check': 'yeah'}
        
        self.data['gfile']['g'] = g
        self.data['gfile']['gcomp'] = gfiledic
        
        return psiNinterp
    
    
    def calcpsi(self):
        # for j in range(len(Y)):
        #     for i in range(len(X)):
        #         PsinLoc.values[j,i,n] = GF.psiN(RadLoc.loc[Y[j],X[i],Attempt].values,VertLoc.loc[Y[j],X[i],Attempt].values,)
        
        geo = self.data['b2fgeo']
        pol_range = int(self.data['b2fgeo']['nx'] + 2)
        print('xdim is {}'.format(pol_range))
        rad_range = int(self.data['b2fgeo']['ny'] + 2)
        print('ydim is {}'.format(rad_range))
        
        self.data['DefaultSettings']['XDIM'] = pol_range
        self.data['DefaultSettings']['YDIM'] = rad_range
        # gR = self.data['gfile']['gcomp']['gR']
        # gZ = self.data['gfile']['gcomp']['gZ']
        # psiN = self.data['gfile']['gcomp']['psiN']
        psiNinterp = self.load_solpsgeo()
        psival = np.zeros((pol_range, rad_range))
        for pol_loc in range(pol_range):
            crLowerLeft = geo['crx'][pol_loc,:,0]
            crUpperLeft = geo['crx'][pol_loc,:,2]
            czLowerLeft = geo['cry'][pol_loc,:,0]
            czUpperLeft = geo['cry'][pol_loc,:,2]
            
            # self.data['psi']['crLowerLeft'] = crLowerLeft
            # self.data['psi']['czLowerLeft'] = czLowerLeft
            # self.data['psi']['crUpperLeft'] = crUpperLeft
            # self.data['psi']['czUpperLeft'] = czUpperLeft
            
            psi_solps = np.zeros(rad_range)
            for i in range(rad_range):
                psi_solps_LL = psiNinterp(crLowerLeft[i], czLowerLeft[i])
                psi_solps_UL = psiNinterp(crUpperLeft[i], czUpperLeft[i])
                psi_solps[i] = np.mean([psi_solps_LL,psi_solps_UL])
            
            psival[pol_loc, :] = psi_solps
            
        self.data['psi']['psival'] = psival
    
    def calcpsi_sp(self, pol_loc):
        
        geo = self.data['b2fgeo']
        pol_range = int(self.data['b2fgeo']['nx'] + 2)
        print('xdim is {}'.format(pol_range))
        rad_range = int(self.data['b2fgeo']['ny'] + 2)
        print('ydim is {}'.format(rad_range))
        
        self.data['DefaultSettings']['XDIM'] = pol_range
        self.data['DefaultSettings']['YDIM'] = rad_range
        
        psiNinterp = self.load_solpsgeo()
        psival = np.zeros((pol_range, rad_range))
        
        crLowerLeft = geo['crx'][pol_loc,:,0]
        crUpperLeft = geo['crx'][pol_loc,:,2]
        czLowerLeft = geo['cry'][pol_loc,:,0]
        czUpperLeft = geo['cry'][pol_loc,:,2]
            
            # self.data['psi']['crLowerLeft'] = crLowerLeft
            # self.data['psi']['czLowerLeft'] = czLowerLeft
            # self.data['psi']['crUpperLeft'] = crUpperLeft
            # self.data['psi']['czUpperLeft'] = czUpperLeft
            
            psi_solps = np.zeros(rad_range)
            for i in range(rad_range):
                psi_solps_LL = psiNinterp(crLowerLeft[i], czLowerLeft[i])
                psi_solps_UL = psiNinterp(crUpperLeft[i], czUpperLeft[i])
                psi_solps[i] = np.mean([psi_solps_LL,psi_solps_UL])
            
            psival[pol_loc, :] = psi_solps
            
        self.data['psi']['psival'] = psival
    
    
    
            
    
    def load_vessel(self):
        # try:
        #     WallFile = np.loadtxt('{}/mesh.extra'.format(self.data['dirdata']['tbase']))
        # except:
        #     print('mesh.extra file not found! Using vvfile.ogr instead')
        #     WallFile=None
            
        VVFILE = np.loadtxt('{}/vvfile.ogr'.format(self.data['dirdata']['tbase']))
        
        self.data['vessel'] = VVFILE
        
        
    def creat_grid(self):
        N = len([self.data['dircomp']['a_shift']])
        print(N)
        ast = self.data['dircomp']['a_shift']
        shift = self.data['dircomp']['shiftdic'][ast]
        Attempts = [str(shift)]
        XGrid = self.data['b2fgeo']['nx']
        print(XGrid)
        XMin= 1
        XMax= XGrid
        X_Core = int(self.data['b2fgeo']['rightcut'] - self.data['b2fgeo']['leftcut'])
        
        CoreBound = []
        CoreBound.append(int(self.data['b2fgeo']['leftcut']))
        CoreBound.append(int(self.data['b2fgeo']['rightcut'] -1))
        
        YSurf = int(self.data['b2fgeo']['ny'])
    
        # Create X and Y mesh grid arrays
        
        X = np.linspace(XMin,XMax,XGrid)
        Y = np.linspace(1,YSurf,YSurf)
        Xx, Yy = np.meshgrid(X,Y)
        
        self.data['gridsettings']['CoreBound'] = CoreBound
        self.data['gridsettings']['XMin'] = XMin
        self.data['gridsettings']['XMax'] = XMax
        self.data['gridsettings']['X'] = X
        self.data['gridsettings']['Y'] = Y
        self.data['gridsettings']['Xx'] = Xx
        self.data['gridsettings']['Yy'] = Yy
        
        Ya = YSurf + 1
        Xa = X_Core + 1
        
        
        RadLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), 
                              coords=[Y,X,Attempts], 
                        dims=['Radial_Location','Poloidal_Location','Attempt'], 
                              name = r'Radial Coordinate $m$')
        VertLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), 
                               coords=[Y,X,Attempts], 
                        dims=['Radial_Location','Poloidal_Location','Attempt'], 
                               name = r'Vertical Coordinate $m$')
        
        Core_Corners = xr.DataArray(np.zeros((Ya, Xa, N,2)), 
                                    coords=[np.concatenate(([0],Y)), np.linspace(CoreBound[0],CoreBound[1]+1,X_Core+1),Attempts,['X','Y']], 
                                    dims=['Radial Index','Poloidal Index','Attempt','Point Coordinates'], 
                                    name = r'Core Corner Coordinates $m$')
        Div1_Corners = xr.DataArray(np.zeros((YSurf+1,CoreBound[0]+1,N,2)), 
                                    coords=[np.concatenate(([0],Y)), 
                np.linspace(0,CoreBound[0],CoreBound[0]+1),Attempts,['X','Y']], 
          dims=['Radial Index','Poloidal Index','Attempt','Point Coordinates'], 
                              name = r'Inner Divertor Corner Coordinates $m$')
        Div2_Corners = xr.DataArray(np.zeros((YSurf+1,XMax-CoreBound[1],N,2)), 
                                    coords=[np.concatenate(([0],Y)), 
        np.linspace(CoreBound[1]+1,XMax,XMax-CoreBound[1]),Attempts,['X','Y']], 
        dims=['Radial Index','Poloidal Index','Attempt','Point Coordinates'], 
                              name = r'Outer Divertor Corner Coordinates $m$')

        YYLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], 
                        dims=['Radial_Location','Poloidal_Location','Attempt'], 
                             name = r'Radial Grid Point $N$')
        
        PsinLoc = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], 
                        dims=['Radial_Location','Poloidal_Location','Attempt'], 
                               name = r'Normalized Psi $\psi_N$')
        
        PolLbl = ['XXLoc', 'Theta', 'dXP','dXP_norm']
        
        PolVec = xr.DataArray(np.zeros((YSurf,XGrid,N,4)), 
                              coords=[Y,X,Attempts,PolLbl], 
        dims=['Radial_Location','Poloidal_Location','Attempt','Poloidal Metric'], 
                              name = 'Poloidal Coordinate Data')
        
        grid_dic = {'RadLoc': RadLoc, 'VertLoc': VertLoc, 
                    'Core_Corners': Core_Corners, 'Div1_Corners': Div1_Corners,
                    'Div2_Corners': Div2_Corners, 'YYLoc': YYLoc,
                    'PsinLoc': PsinLoc, 'PolLbl': PolLbl, 'PolVec': PolVec}
        
        # self.data['grid'] = grid_dic
        return grid_dic
        
        
    def load_output_geo(self, grid_dic):
        n = 0
        DRT = self.data['dirdata']['adir']['outputdir']['Output']
        DRT2 = self.data['dirdata']['adir']['outputdir']['Output2']
        XDIM = self.data['b2fgeo']['nx'] + 2
        YDIM = self.data['b2fgeo']['ny'] + 2
        Attempt = self.data['dircomp']['Attempt']
        
        self.data['DefaultSettings']['XDIM'] = XDIM
        self.data['DefaultSettings']['YDIM'] = YDIM
        XMin = self.data['gridsettings']['XMin']
        XMax = self.data['gridsettings']['XMax']
        CoreBound = self.data['gridsettings']['CoreBound']
        Yy = self.data['gridsettings']['Yy']
        Core_Corners = grid_dic['Core_Corners']
        YYLoc = grid_dic['YYLoc']
        RadLoc = grid_dic['RadLoc']
        VertLoc = grid_dic['VertLoc']
        Div1_Corners = grid_dic['Div1_Corners']
        Div2_Corners = grid_dic['Div2_Corners']
        
                    
        YYLoc.values[:,:, n] = Yy
        RadLoc.values[:,:,n] = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                    usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
        VertLoc.values[:,:,n] = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                      usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
        
        Rad0Cor = np.loadtxt('{}/Rad0Cor{}'.format(DRT2, str(Attempt)), 
                    usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
        Vert0Cor = np.loadtxt('{}/Vert0Cor{}'.format(DRT2, str(Attempt)), 
                    usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
        
        Rad1Cor = np.loadtxt('{}/Rad1Cor{}'.format(DRT2, str(Attempt)), 
              usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
        Vert1Cor = np.loadtxt('{}/Vert1Cor{}'.format(DRT2, str(Attempt)), 
              usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]

        Rad2Cor = np.loadtxt('{}/Rad2Cor{}'.format(DRT2, str(Attempt)), 
              usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
        Vert2Cor = np.loadtxt('{}/Vert2Cor{}'.format(DRT2, str(Attempt)), 
              usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]

        Rad3Cor = np.loadtxt('{}/Rad3Cor{}'.format(DRT2, str(Attempt)), 
              usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
        Vert3Cor = np.loadtxt('{}/Vert3Cor{}'.format(DRT2, str(Attempt)), 
              usecols = (3)).reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
        
        Core_Corners.values[:-1,:-1,n,0] = Rad0Cor[:,CoreBound[0]:CoreBound[1]+1]
        Core_Corners.values[:-1,-1,n,0] = Rad1Cor[:,CoreBound[1]]
        Core_Corners.values[-1,:-1,n,0] = Rad2Cor[-1,CoreBound[0]:CoreBound[1]+1]
        Core_Corners.values[-1,-1,n,0] = Rad3Cor[-1,CoreBound[1]]
        
        Core_Corners.values[:-1,:-1,n,1] = Vert0Cor[:,CoreBound[0]:CoreBound[1]+1]
        Core_Corners.values[:-1,-1,n,1] = Vert1Cor[:,CoreBound[1]]
        Core_Corners.values[-1,:-1,n,1] = Vert2Cor[-1,CoreBound[0]:CoreBound[1]+1]
        Core_Corners.values[-1,-1,n,1] = Vert3Cor[-1,CoreBound[1]]
        
        Div1_Corners.values[:-1,:-1,n,0] = Rad0Cor[:,XMin-1:CoreBound[0]]
        Div1_Corners.values[:-1,-1,n,0] = Rad1Cor[:,CoreBound[0]-1]
        Div1_Corners.values[-1,:-1,n,0] = Rad2Cor[-1,XMin-1:CoreBound[0]]
        Div1_Corners.values[-1,-1,n,0] = Rad3Cor[-1,CoreBound[0]-1]
        
        Div1_Corners.values[:-1,:-1,n,1] = Vert0Cor[:,XMin-1:CoreBound[0]]
        Div1_Corners.values[:-1,-1,n,1] = Vert1Cor[:,CoreBound[0]-1]
        Div1_Corners.values[-1,:-1,n,1] = Vert2Cor[-1,XMin-1:CoreBound[0]]
        Div1_Corners.values[-1,-1,n,1] = Vert3Cor[-1,CoreBound[0]-1]
        
        Div2_Corners.values[:-1,:-1,n,0] = Rad0Cor[:,CoreBound[1]+1:]
        Div2_Corners.values[:-1,-1,n,0] = Rad1Cor[:,-1]
        Div2_Corners.values[-1,:-1,n,0] = Rad2Cor[-1,CoreBound[1]+1:]
        Div2_Corners.values[-1,-1,n,0] = Rad3Cor[-1,-1]
        
        Div2_Corners.values[:-1,:-1,n,1] = Vert0Cor[:,CoreBound[1]+1:]
        Div2_Corners.values[:-1,-1,n,1] = Vert1Cor[:,-1]
        Div2_Corners.values[-1,:-1,n,1] = Vert2Cor[-1,CoreBound[1]+1:]
        Div2_Corners.values[-1,-1,n,1] = Vert3Cor[-1,-1]
        
        # self.data['grid']['Core_Corners'] = Core_Corners.values[:, :, 0, :]
        # self.data['grid']['YYLoc'] = YYLoc.values[:, :, 0]
        # self.data['grid']['RadLoc'] = RadLoc.values[:, :, 0]
        # self.data['grid']['VertLoc'] = VertLoc.values[:, :, 0]
        # self.data['grid']['Div1_Corners'] = Div1_Corners.values[:, :, 0, :] 
        # self.data['grid']['Div2_Corners'] = Div2_Corners[:, :, 0, :]
        
        self.data['grid']['Core_Corners'] = Core_Corners
        self.data['grid']['YYLoc'] = YYLoc
        self.data['grid']['RadLoc'] = RadLoc
        self.data['grid']['VertLoc'] = VertLoc
        self.data['grid']['Div1_Corners'] = Div1_Corners 
        self.data['grid']['Div2_Corners'] = Div2_Corners
        
        
    
    def pop_pol_coord(self, grid_dic):
        #Populate Poloidal Coordinate Data array
        n = 0
        Attempt = len([self.data['dircomp']['a_shift']])
        # print(Attempt)
        Xx = self.data['gridsettings']['Xx']
        PolVec = grid_dic['PolVec']
        X = self.data['gridsettings']['X']
        RadLoc = self.data['grid']['RadLoc'] 
        VertLoc = self.data['grid']['VertLoc'] 
        PolVec.loc[:,:, 0,'XXLoc'] = Xx
        JXA = self.data['b2mn']['jxa']
        CoreBound = self.data['DefaultSettings']['CoreBound']
        YDIM = int(self.data['DefaultSettings']['YDIM'])
        
        if YDIM % 2 == 0:
            SEP = YDIM/ 2 + 1
        else:
            SEP = round(YDIM/ 2)
        
        YVector=np.zeros((len(X),2))
        YVector[:,0] = RadLoc.values[1,:,n] - RadLoc.values[0,:,n]
        YVector[:,1] = VertLoc.values[1,:,n] - VertLoc.values[0,:,n]
                    
    
        for i in range(len(X)):
            PolVec.loc[:,X[i],Attempt,'Theta'] = np.degrees(np.math.atan2(np.linalg.det([YVector[JXA-1,:],YVector[i,:]]),np.dot(YVector[JXA-1,:],YVector[i,:])))
            if PolVec.loc[:,X[i],Attempt,'Theta'].values[0] < 0 and X[i] < JXA:
                PolVec.loc[:,X[i],Attempt,'Theta'] = PolVec.loc[:,X[i],Attempt,'Theta'] + 360          
        
        XP_range=np.array([CoreBound[0]-1,CoreBound[0],CoreBound[1],CoreBound[1]+1])
        X_xp=np.mean(RadLoc.loc[SEP,XP_range,Attempt].values)
        Y_xp=np.mean(VertLoc.loc[SEP,XP_range,Attempt].values)
        Xpoint=np.array([X_xp,Y_xp])
        
        for index in np.arange(CoreBound[0],CoreBound[1]+1):
            if index == CoreBound[0]:
                PolVec.loc[:,index,Attempt,'dXP'] = round(np.sqrt((RadLoc.loc[SEP,index,Attempt].values-X_xp)**2 + (VertLoc.loc[SEP,index,Attempt].values-Y_xp)**2),5)
            else:
                NL = np.sqrt((RadLoc.loc[SEP,index,Attempt].values-RadLoc.loc[SEP,index-1,Attempt].values)**2 + (VertLoc.loc[SEP,index,Attempt].values-VertLoc.loc[SEP,index-1,Attempt].values)**2)
                PolVec.loc[:,index,Attempt,'dXP']= PolVec.loc[:,index-1,Attempt,'dXP']+NL
     
        PolVec.loc[:,:,Attempt,'dXP_norm'] = PolVec.loc[:,:,Attempt,'dXP'].values/np.max(PolVec.loc[:,:,Attempt,'dXP'].values)
        
        self.data['grid']['PolVec'] = PolVec



    def data_dense(self):
        
        if AddNew is not None:
            if AddNew not in self.Parameter:
                self.Parameter.append(AddNew)
                P = len(self.Parameter)
        
        for p in self.Parameter:
            if p not in self.PARAM.keys() or 'AVG' in Attempts: 
                self.PARAM[p] = xr.DataArray(np.zeros((YSurf,XGrid,N)), coords=[Y,X,Attempts], dims=['Radial_Location','Poloidal_Location','Attempt'], name = self.PARAMDICT[p])
                for n in range(N):
                    Attempt = Attempts[n]
                    if Attempt == 'AVG':
                        self.PARAM[p].values[:,:,n] = self.PARAM[p].values[:,:,:-1].mean(2)
                    else:    
                        DRT = '{}/Attempt{}'.format(BASEDRT, str(Attempt))   #Generate path
                        try:
                            RawData = np.loadtxt('{}/Output/{}{}'.format(DRT, p, str(Attempt)),usecols = (3))
                        except Exception as err:
                            print(err)
                            try:
                                 RawData = np.loadtxt('{}/Output2/{}{}'.format(DRT, p, str(Attempt)),usecols = (3))
                            except Exception as err:
                                print(err)
                                print('Parameter {} not found for Attempt {}. Creating NAN Array'.format(p, str(Attempt)))
                                self.PARAM[p].values[:,:,n] = np.nan
                                RawData=[]
                        
                        if len(RawData) > 0:        
                            if RawData.size == XDIM*YDIM:
                                self.PARAM[p].values[:,:,n] = RawData.reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
                            elif RawData.size == XDIM*YDIM*2:
                                self.PARAM[p].values[:,:,n] = RawData.reshape((2*YDIM,XDIM))[1+YDIM:2*YDIM-1,XMin:XMax+1]
                                

            

    def set_plot(self):
        Publish = 'b2plottersetting'
        if Publish == 'b2plottersetting':
            plt.rcParams.update({'font.weight': 'bold'})
            plt.rcParams.update({'font.size': 16})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
        
        self.data['DefaultSettings']['Publish'] = Publish
        # if Publish != []:
        #     plt.rc('font',size=30)
        #     plt.rc('lines',linewidth=5,markersize=15)
        # else:
        #     plt.rc('font',size=14)
    
    def set_rad_pol(self):
        if RadSlc == 'all':
            RadSlc = self.PARAM.coords['Radial_Location'].values
        if RadSlc == None:
            RadSlc = SEP
            
        if PolSlc == 'all':
            PolSlc = self.PARAM.coords['Poloidal_Location'].values
        if PolSlc == None:
            PolSlc = JXA
            
            

class load_expdata(B2plotter):
    
    def __init__(self, DEV, DefaultSettings, loadDS):
        B2plotter.__init__(self, DEV, DefaultSettings)
        # Employee.__init__(self, first, last, pay)
        self.loadDS = loadDS
        
    def loadmastdata(self):
        if self.loadDS['EXP']:
            mastloc = '{}/{}'.format(self.data['dirdata']['gbase'],
                                     self.loadDS['expfilename'])
            expdic = mdm.read_mastfile(mastloc)
            self.data['ExpDict'] = expdic
            self.data['dirdata']['mastloc'] = mastloc
        
        if self.loadDS['fit']:
            fitloc = '{}/{}/{}'.format(self.data['dirdata']['basedrt'], 
                                    self.DEV, self.loadDS['fitfname'])
            fitdic = mdm.read_fitfile(fitloc)
            self.data['fitprofile'] = fitdic
            self.data['dirdata']['fitloc'] = fitloc
        
        
    # N = len(Attempts)
    # print('{} Attempt(s) Entered'.format(N))
    # P = len(self.Parameter)
            
            
class load_data(load_expdata):
    
    def __init__(self, DEV, DefaultSettings, loadDS, Parameters):
        load_expdata.__init__(self, DEV, DefaultSettings, loadDS)
        # Employee.__init__(self, first, last, pay)
        
        "Parameters"
        if isinstance(Parameters, dict):
            self.Parameters = Parameters
        else:
            print('parameter has to be a dictionary')
            
        if Parameters is None:
            print('There is no parameters input')
        else:
            self.Parameters = Parameters
            
        Plist = []
        for pkey, pvalue in self.Parameters.items():
            Plist.append(pkey)
        
        self.data['paramkey'] = Plist
        
    "Add and remove elements from parameter or defaultsettings"
    def add_dic(self, new_set, assign='default'):
        if assign == 'param':
            self.Parameters = new_set | self.Parameters
            Plist = []
            for pkey, pvalue in self.Parameters.items():
                Plist.append(pkey)
        
        # elif assign == 'default':
        #     self.DefaultSettings = new_set | self.DefaultSettings
        #     keylist = []
        #     for key, value in self.DefaultSettings.items():
        #         keylist.append(key)
        else:
            print('assign parameter is incorrect')
            
        
    def remove_dic(self, new_set, assign='param'):
        if assign == 'param':
            if new_set.keys() in self.data['defaultkey']:
                del self.Parameters[new_set.keys()]
            Plist = []
            for pkey, pvalue in self.Parameters.items():
                Plist.append(pkey)
                
        # elif assign == 'default':
        #     if new_set.keys() in self.data['defaultkey']:
        #         del self.DefaultSettings[new_set.keys()]
        #     keylist = []
        #     for key, value in self.DefaultSettings.items():
        #         keylist.append(key)
        else:
            print('assign parameter incorrect')

    
    def load_output_data(self, param):
        BASEDRT = self.data['dirdata']['adir']['outputdir']['Output']
        Attempt = self.data['dircomp']['Attempt']
        Attempts = len([self.data['dircomp']['a_shift']])
        N = len(self.data['dircomp']['a_shift'])
        XGrid = int(self.data['b2fgeo']['nx'])
        X = self.data['gridsettings']['X']
        XMin= 1
        XMax= XGrid
        # print(XGrid)
        XDIM = int(self.data['DefaultSettings']['XDIM'])
        YSurf = int(self.data['b2fgeo']['ny'])
        Y = self.data['gridsettings']['Y']
        YDIM = int(self.data['DefaultSettings']['YDIM'])
        n = 0
        
        
        # DRT = '{}/Attempt{}'.format(BASEDRT, str(Attempt))   #Generate path
        test = param in self.Parameters.keys()
        self.data['outputdata'][param] = np.zeros([YDIM, XDIM], dtype= np.float32)
 #        self.data['outputdata'][param] = xr.DataArray(np.zeros((YSurf,XGrid,N)), 
 #                                     coords=[Y,X,Attempts], 
 # dims=['Radial_Location','Poloidal_Location','Attempt'], name = param)
        if test:
            print('yes, {} is in parameter'.format(param))
            RawData = np.loadtxt('{}/{}{}'.format(BASEDRT, param, str(Attempt)),usecols = (3))
        elif test == False:
            print('no, {} is not in parameter'.format(param))
        else:
            print('there might be a bug')
        
        if len(RawData) > 0:        
            if RawData.size == XDIM*YDIM:
                # self.data['outputdata'][param].values[:,:,n] = RawData.reshape((YDIM,XDIM))[1:YDIM-1,XMin:XMax+1]
                self.data['outputdata'][param] = RawData.reshape((YDIM,XDIM))
            elif RawData.size != XDIM*YDIM:
                print('rawdata size is not equal to {}'.format(str(XDIM*YDIM)))
            # elif RawData.size == XDIM*YDIM*2:
            #     self.data['outputdata'][param].values[:,:,n] = RawData.reshape((2*YDIM,XDIM))[1+YDIM:2*YDIM-1,XMin:XMax+1]
        else:
            print('we have a problem loading rawdata')
            
            
class single_plot(load_data):
    def __init__(self, DEV, DefaultSettings, loadDS, Parameters):
        load_data.__init__(self, DEV, DefaultSettings, loadDS, Parameters)
    
    def simple_plot(self):
        
        
    
    
    

            
            
            
    
        
        
        
            
            
        
    
