# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 18:24:05 2023

@author: user
"""






    def load_vessel(self):
        # try:
        #     WallFile = np.loadtxt('{}/mesh.extra'.format(self.data['dirdata']['tbase']))
        # except:
        #     print('mesh.extra file not found! Using vvfile.ogr instead')
        #     WallFile=None
            
        VVFILE = np.loadtxt('{}/vvfile.ogr'.format(self.data['dirdata']['simutop']))
        
        self.data['vessel'] = VVFILE
        
        # if plot:
        #     plt.plot