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
from scipy.interpolate import interp1d, griddata
from D3DPreProcess import PsiNtoR
from D3DPreProcess import RhotoPsiN
import B2plotter_tool as b2t
import mastload_data_method as mdm
import coord_sutils as cs
import coord_plot as cp
from scipy import interpolate


class B2plotter:
    def __init__(self, DEV, Shot, shift, series, Attempts, Parameters, 
                 DefaultSettings):
        
        self.DEV = DEV   
        self.Shot = Shot
        self.shift = shift
        self.series = series
        
        "Attempts"
        if isinstance(Attempts, list):
            self.Attempts = Attempts
        else:    
            print('Attempt has to be a list')
            
        
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
        self.data = {'paramkey': Plist, 'defaultkey':keylist, 
                     'ExpDict': {}, 'RadCoords':{}, 'solpsdata': {}}

                    
    "Add and remove elements from parameter or defaultsettings"
    def add_dic(self, new_set, assign='default'):
        if assign == 'param':
            self.Parameters = new_set | self.Parameters
            Plist = []
            for pkey, pvalue in self.Parameters.items():
                Plist.append(pkey)
        
        elif assign == 'default':
            self.DefaultSettings = new_set | self.DefaultSettings
            keylist = []
            for key, value in self.DefaultSettings.items():
                keylist.append(key)
        else:
            print('assign parameter is incorrect')
        
    
    def remove_dic(self, new_set, assign='param'):
        if assign == 'param':
            if new_set.keys() in self.data['defaultkey']:
                del self.Parameters[new_set.keys()]
            Plist = []
            for pkey, pvalue in self.Parameters.items():
                Plist.append(pkey)
                
        elif assign == 'default':
            if new_set.keys() in self.data['defaultkey']:
                del self.DefaultSettings[new_set.keys()]
            keylist = []
            for key, value in self.DefaultSettings.items():
                keylist.append(key)
        else:
            print('assign parameter incorrect')
        
        
    def _reset_object(self):
        self.Shot=None
        self.Attempts=None
        self.Parameter=[]
        self.PARAM={}
        self.ExpDict={}
        self.RadCoords={}
    
        
class mastdata(B2plotter):
    
    def __init__(self, DEV, Shot, shift, series, 
                 Attempts, Parameters, loadDS, DefaultSettings):
        B2plotter.__init__(self, DEV, Shot, shift, series, Attempts, Parameters, DefaultSettings)
        # Employee.__init__(self, first, last, pay)
        self.loadDS = loadDS
        
    def loadmastdata(self, a_shift):
        ROOTSHOT = self.loadDS['ROOTSHOT']
        
        if self.DEV == 'mast':
            mastdic = mdm.mast_b2_dir(a_shift)
            self.data['dirdic'] = mastdic
        
        if self.loadDS['EXP']:
            mastloc = '{}/{}'.format(self.data['dirdic']['gbase'],
                                     self.loadDS['expfilename'])
            expdic = mdm.read_mastfile(mastloc)
            self.data['ExpDict'] = expdic
            self.data['dirdic']['mastloc'] = mastloc
        
        if self.loadDS['fit']:
            fitloc = '{}/{}/{}'.format(self.data['dirdic']['basedrt'], 
                                    self.DEV, self.loadDS['fitfname'])
            fitdic = mdm.read_fitfile(fitloc)
            self.data['fitprofile'] = fitdic
            self.data['dirdic']['fitloc'] = fitloc
    def mastcalcpsi(self, a_shift='org', geo=None, b2mn=None, 
                   dsa=None, shift= 0):
            if dsa is None:
                dsa = cs.read_dsa(self.data['dirdic']['simudir'] + '/dsa')
            if geo is None:
                geo = cs.read_b2fgmtry(self.data['dirdic']['tbase'] + '/baserun/b2fgmtry')
            if b2mn is None:
                b2mn = cs.scrape_b2mn(self.data['dirdic']['simudir'] + '/b2mn.dat')
            
            crLowerLeft = geo['crx'][b2mn['jxa']+1,:,0]
            crUpperLeft = geo['crx'][b2mn['jxa']+1,:,2]
            czLowerLeft = geo['cry'][b2mn['jxa']+1,:,0]
            czUpperLeft = geo['cry'][b2mn['jxa']+1,:,2]
            
            self.data['geo'] = geo
            self.data['b2mn'] = b2mn   
            
            ncells = len(czLowerLeft)

            g = cs.loadg(self.data['dirdic']['gdir'][0])
            psiN = (g['psirz'] - g['simag']) / (g['sibry'] - g['simag'])

            dR = g['rdim'] / (g['nw'] - 1)
            dZ = g['zdim'] / (g['nh'] - 1)

            gR = []
            for i in range(g['nw']):
                gR.append(g['rleft'] + i * dR + float(shift))

            gZ = []
            for i in range(g['nh']):
                gZ.append(g['zmid'] - 0.5 * g['zdim'] + i * dZ)

            gR = np.array(gR)
            gZ = np.array(gZ)

            R_solps_top = 0.5 * (np.array(crLowerLeft) + np.array(crUpperLeft))
            Z_solps_top = 0.5 * (np.array(czLowerLeft) + np.array(czUpperLeft))

            psiNinterp = interpolate.interp2d(gR, gZ, psiN, kind = 'cubic')

            psi_solps = np.zeros(ncells)
            for i in range(ncells):
                psi_solps_LL = psiNinterp(crLowerLeft[i], czLowerLeft[i])
                psi_solps_UL = psiNinterp(crUpperLeft[i], czUpperLeft[i])
                psi_solps[i] = np.mean([psi_solps_LL,psi_solps_UL])
                
            psi_list = psi_solps.tolist()
            RLL_list = crLowerLeft.tolist()
            ZLL_list = czLowerLeft.tolist()
            RUL_list = crUpperLeft.tolist()
            ZUL_list = czUpperLeft.tolist()
            
            print(type(ZUL_list))
            print(type(dsa))
            print(type(psi_list))
            

            self.data['solpsdata']['crLowerLeft'] = RLL_list
            self.data['solpsdata']['czLowerLeft'] = ZLL_list
            self.data['solpsdata']['crUpperLeft'] = RUL_list
            self.data['solpsdata']['czUpperLeft'] = ZUL_list
            self.data['solpsdata']['dsa'] = dsa
            self.data['solpsdata']['psiSOLPS'] = psi_list
            
            datakey = ['crLowerLeft','czLowerLeft', 'crUpperLeft','czUpperLeft', 
                       'dsa', 'psiSOLPS']
            cn = len(datakey)
            print(len(RLL_list))
            print(len(ZLL_list))
            print(len(RUL_list))
            print(len(ZUL_list))
            print(len(dsa))
            print(len(psi_list))
            
            dataindex = [RLL_list, ZLL_list, RUL_list, ZUL_list, dsa, psi_list]
            
            
            with open('C:/Users/user/Documents/GitHub/efold_test/coord.txt', 'w') as file:
                colcount = 0    # Track the column number

                for x in datakey:
                    # First thing we do is add one to the column count when
                    # starting the loop. Since we're doing some math on it below
                    # we want to make sure we don't divide by zero.
                    colcount += 1
                
                    # After each entry, add a tab character ("\t")
                    file.write(x + "\t\t\t")
                
                    # Now, check the column count against the MAX_COLUMNS. We
                    # use a modulus operator (%) to get the remainder after dividing;
                    # any number divisible by 3 in our example will return '0'
                    # via modulus.
                    if colcount == cn:
                        # Now write out a new-line ("\n") to move to the next line.
                        file.write("\n")

                for i in range(len(dsa)):
                    for p in dataindex:
                        file.write(str(p[i]) + "\t\t\t")
                    file.write("\n")
        
        
        
            
            
        
    
