# -*- coding: utf-8 -*-
"""
Created on Sun Jul  6 16:57:11 2025

@author: ychuang
"""



import numpy as np
from load_coordinate.load_coord_method import load_coordgeo_method
from load_coordinate.SOLPSplotter_geo import load_geometry
from scipy import interpolate


class CrossMachine_load_geo:
    
    def __init__(self, DF, data, lcm: load_coordgeo_method, lg: load_geometry):
        
        self.DF = DF
        self.data = data
        self.lcm = lcm
        self.lg = lg
    
    
    def CM_load_solpsgeo(self):
    
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        terminal = self.DF.terminal
        
        if self.DF.DEV == 'cross_machine' and self.DF.Dnames == 'mastu_mast':
            
            device_list = ['mast', 'mastu']
            
            for kk in device_list:
                
                if terminal == True:
                    
                    print('check the g file dir')
                    print(self.data['{}_dirdata']['gdir'].format(kk))
                    g_loc = self.data['{}_dirdata'.format(kk)]['gdir'][0]
                    
                elif terminal == False:
                    g_loc = self.data['{}_dirdata'.format(kk)]['gdir']
                
                print('at CM_load_solpsgeo, this is {}'.format(kk))
                print(g_loc)
                gfile_data = self.lcm.loadg(g_loc)
                g_dic = {'g': gfile_data}
                
                
                if withshift == False and withseries == False:
                    simudir = self.data['{}_dirdata'.format(kk)]['simudir']
                    simutop = self.data['{}_dirdata'.format(kk)]['simutop']
                    shift = self.data['{}_dircomp'.format(kk)]['shift']
                    
                    if kk == 'mast':
                        
                        b2mn, geo, gfilesum = self.lg.loadgeo_method(attempt_loc = simudir, 
                                simufile_loc = simutop, g_data = gfile_data, shift_value = shift)
                    
                    elif kk == 'mastu':
                        
                        b2mn, geo, gfilesum = self.lg.loadsymgeo_method(attempt_loc = simudir, 
                                simufile_loc = simutop, g_data = gfile_data, shift_value = shift)
                    
                    
                    self.data['{}_b2mn'.format(kk)] = b2mn
                    self.data['{}_b2fgeo'.format(kk)] = geo
                    
                    g_dic['gcomp'] = gfilesum
                    
                    self.data['{}_gfile'.format(kk)] = g_dic
    
    
    
    
    def CM_calcpsi_avcr(self):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        if self.DF.DEV == 'cross_machine' and self.DF.Dnames == 'mastu_mast':
            
            device_list = ['mast', 'mastu']
            
            for kk in device_list:
                
                
                if withshift == False and withseries == False:
                    
                    b2fgeo = self.data['{}_b2fgeo'.format(kk)]
                    psiNinterp_RBS = self.data['{}_gfile'.format(kk)]['gcomp']['interp_dic']['RBS']
                    
                    
                    RadLoc, VertLoc, psival, pol_range, rad_range = self.lg.calcpsi_method_avcr(geo = b2fgeo, 
                                                                            psi_RBS = psiNinterp_RBS)
                    
                    
                    coord_dic = {'RadLoc': RadLoc, 'VertLoc': VertLoc}
                    self.data['{}_grid'.format(kk)] = coord_dic
                    psi_dic = {'psival': psival}
                    self.data['{}_psi'.format(kk)] = psi_dic
                    range_dic = {'XDIM': pol_range, 'YDIM': rad_range}
                    self.data['{}_DefaultSettings'.format(kk)] = range_dic
        
        
        
        
        
        

                        
                        
                        
                        
                    
                    
                    

                    

                
            
            
        
        
    
""" 
    
    
    def load_solpsgeo(self):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        terminal = self.DF.terminal
        
        
        if terminal == True:
            
            print('check the g file dir')
            print(self.data['dirdata']['gdir'])
            g_loc = self.data['dirdata']['gdir'][0]
            
            
        elif terminal == False:
            g_loc = self.data['dirdata']['gbase'] + '/MAST__RMP_results/g027205.00275_efitpp'
            

        gfile_data = self.lcm.loadg(g_loc)
        self.data['gfile']['g'] = gfile_data
        
        
        if withshift == False and withseries == False:
            simudir = self.data['dirdata']['simudir']
            simutop = self.data['dirdata']['simutop']
            shift = self.data['dircomp']['shift_value']
            
            b2mn, geo, gfilesum = self.loadgeo_method(attempt_loc = simudir, 
                    simufile_loc = simutop, g_data = gfile_data, shift_value = shift)
            
            self.data['b2mn'] = b2mn
            self.data['b2fgeo'] = geo
            self.data['gfile']['gcomp'] = gfilesum


"""








