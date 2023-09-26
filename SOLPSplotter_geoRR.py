# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:21:40 2023

@author: user
"""
from B2plotter_class import B2plotter
import xarray as xr
import matplotlib.pyplot as plt
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate
from scipy.optimize import curve_fit
from matplotlib.path import Path
from scipy.stats import binned_statistic
import numpy as np


class geo_RR(B2plotter):
    def __init__(self, DEV, withshift, withseries, DefaultSettings):
        B2plotter.__init__(self, DEV, withshift, withseries, DefaultSettings)


    def creat_grid(self):
        N = len([self.data['dircomp']['a_shift']])
        print(N)
        ast = self.data['dircomp']['a_shift']
        shift = self.data['dircomp']['shift_dic'][ast]
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
        DRT = self.data['dirdata']['outputdir']['Output']
        DRT2 = self.data['dirdata']['outputdir']['Output2']
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

    
    def dsa_1D(self, pol_loc):
        
        if self.withshift == False and self.withseries == False:
            geo = self.data['b2fgeo']
            pol_range = int(self.data['b2fgeo']['nx'] + 2)
            # print('xdim is {}'.format(str(pol_range)))
            rad_range = int(self.data['b2fgeo']['ny'] + 2)
            # print('ydim is {}'.format(str(rad_range)))
            
            self.data['DefaultSettings']['XDIM'] = pol_range
            self.data['DefaultSettings']['YDIM'] = rad_range
            
            if rad_range % 2 == 0:
                SEP = rad_range/ 2 + 1
            else:
                SEP = round(rad_range/ 2)
            
            self.data['DefaultSettings']['SEP'] = SEP
    
            dsa = lcm.read_dsa(self.data['dirdata']['simudir'] + '/dsa')
            
            
            psiNinterp_RGI, psiNinterp_2d, psiNinterp_RBS = self.load_solpsgeo()
            # print(type(psiNinterp_RBS))
            psival = np.zeros((pol_range, rad_range))
            
            pol_index = int(pol_loc) + 1
            
            crLowerLeft = geo['crx'][pol_index,:,0]
            crLowerRight = geo['crx'][pol_index,:,1]
            crUpperLeft = geo['crx'][pol_index,:,2]
            crUpperRight = geo['crx'][pol_index,:,3]
            czLowerLeft = geo['cry'][pol_index,:,0]
            czLowerRight = geo['cry'][pol_index,:,1]
            czUpperLeft = geo['cry'][pol_index,:,2]
            czUpperRight = geo['cry'][pol_index,:,3]
                
             
            LLsep = np.mean([crLowerLeft[int(SEP)-1], crLowerLeft[int(SEP)-2]])
            LRsep = np.mean([crLowerRight[int(SEP)-1], crLowerRight[int(SEP)-2]])
            ULsep = np.mean([crUpperLeft[int(SEP)-1], crUpperLeft[int(SEP)-2]])
            URsep = np.mean([crUpperRight[int(SEP)-1], crUpperRight[int(SEP)-2]])
            
            weight = dsa[int(SEP)-1]/ (dsa[int(SEP)-1] - dsa[int(SEP)-2])
            # print(weight)
            
            LLdsa = crLowerLeft - LLsep
            LRdsa = crLowerRight - LRsep
            ULdsa = crUpperLeft - ULsep
            URdsa = crUpperRight - URsep
            avag_rad = np.zeros(rad_range)
            for j in range(rad_range):
                avag_rad[j] = np.mean([LLdsa[j], ULdsa[j], URdsa[j], LRdsa[j]])
          
            
            RadLoc = self.data['grid']['RadLoc']
            VertLoc = self.data['grid']['VertLoc']
            
          
            "P0 = Grab X,Y Point of Inner Most (Core) Cell at PolPos"
          
            # P0 = np.zeros(2)
            
            # P0[0] = np.mean(crLowerLeft[1], crLowerRight[1], 
            #                 crUpperLeft[1], crUpperRight[1])
            # P0[1] = np.mean(czLowerLeft[1], czLowerRight[1], 
            #                 czUpperLeft[1], czUpperRight[1])
            pol_lc = int(pol_loc) + 1
            # print(pol_lc)
            
            # RR_SEP = RadLoc.loc[1, pol_lc, '0'].values
            
            P0 = np.array([RadLoc.loc[1, pol_lc, '0'].values,VertLoc.loc[1, pol_lc, '0'].values])
            
            # "P1 = Grab X,Y Point of Separatrix at PolPos"
            # P1 = np.zeros(2)
            
            # P1[0] = np.mean(LLsep, LRsep, ULsep, URsep)
            
            # LLzsep = np.mean([czLowerLeft[int(SEP)-1], czLowerLeft[int(SEP)-2]])
            # LRzsep = np.mean([czLowerRight[int(SEP)-1], czLowerRight[int(SEP)-2]])
            # ULzsep = np.mean([czUpperLeft[int(SEP)-1], czUpperLeft[int(SEP)-2]])
            # URzsep = np.mean([czUpperRight[int(SEP)-1], czUpperRight[int(SEP)-2]])
            
            # P1[0] = np.mean(LLzsep, ULzsep, URzsep, LRzsep)
            
            P1 = np.array([RadLoc.loc[int(SEP), pol_lc, '0'].values,VertLoc.loc[int(SEP), pol_lc, '0'].values])
            
            R_SEP=np.sqrt((P1[0]-P0[0])**2 + (P1[1]-P0[1])**2)
            
            
            
            PR=RadLoc.loc[:,:, '0'].values.flatten()
            PZ=VertLoc.loc[:,:, '0'].values.flatten()
            
            Mag=0.5
            Thresh=0.01
            PP=P1-P0
            Theta=np.arctan2(PP[1],PP[0])
            displace=Thresh*np.array([-np.sin(Theta),np.cos(Theta)])
            P2=np.array([(Mag*np.cos(Theta)+P0[0]),(Mag*np.sin(Theta)+P0[1])])
            P0A=P0+displace
            P0B=P0-displace
            P2A=P2+displace
            P2B=P2-displace
            
            Bounds = Path([P0A,P2A,P2B,P0B])
            print(type(Bounds))
            Mask = Bounds.contains_points(np.vstack((PR,PZ)).T)
            print(type(Mask))
            ChordX = np.ma.array(PR,mask=~Mask).compressed()
            print(type(ChordX))
            ChordY = np.ma.array(PZ,mask=~Mask).compressed()
            print(type(ChordY))
            RR_SEP = np.sqrt((ChordX-P1[0])**2 + (ChordY-P1[1])**2)*np.sign(np.sqrt((ChordX-P0[0])**2 + (ChordY-P0[1])**2)-R_SEP)
            
            
            RR = binned_statistic(RR_SEP, RR_SEP,statistic='mean',bins=38)[0]
            RR_SEP_avg = RR[~np.isnan(RR)]
            
            # rr_sep = np.zeros((int(rad_range), 2))
            # rr_sep[:, 0] = dsa
            # rr_sep[1:23, 1] = RR_SEP_avg
            
            P3 = np.zeros((36, 2))
            P3[:, 0] = RadLoc.loc[:, pol_lc, '0'].values.T
            P3[:, 1] = VertLoc.loc[:, pol_lc, '0'].values.T
            
            plt.figure(1)
            plt.plot(P3[:, 0], P3[:, 1],'o-', color= 'r', label= 'RR_x_coord')
            plt.plot(crLowerLeft, czLowerLeft, 'o-', color= 'g', label= 'xport_x_coord')
            # plt.xlabel('shift: [m]')
            # plt.ylabel('neutral penetration length: [m]')
            # plt.title('neutral penetration length verses shift distance')
            plt.legend()
            
            dsa_dic = {'dsa_{}_val'.format(pol_loc): P3}
            
            
            # self.data['psi']['dsa_{}_val'.format(pol_loc)] = rr_sep
            self.data['dsa'] = dsa_dic
            # self.data['psi']['RR_sep'] = RR_SEP
            
        else:
            print('more work need to be done')



    

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
                                
                                
    def set_rad_pol(self):
        if RadSlc == 'all':
            RadSlc = self.PARAM.coords['Radial_Location'].values
        if RadSlc == None:
            RadSlc = SEP
            
        if PolSlc == 'all':
            PolSlc = self.PARAM.coords['Poloidal_Location'].values
        if PolSlc == None:
            PolSlc = JXA