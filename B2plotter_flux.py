# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 21:26:02 2023

@author: user
"""
from B2plotter_plot import Opacity_study



class flux_adjustment(Opacity_study):
    def __init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters, Publish):
        Opacity_study.__init__(self, DEV, withshift, withseries, 
                               DefaultSettings, loadDS, Parameters, Publish)
    
    def flux_comparison_plot(self):
        
        pol_range = self.data['b2fgeo']['nx'] + 2
        print('xdim is {}'.format(str(pol_range)))
        rad_range = self.data['b2fgeo']['ny'] + 2
        print('ydim is {}'.format(str(rad_range)))
        
        self.data['DefaultSettings']['XDIM'] = pol_range
        self.data['DefaultSettings']['YDIM'] = rad_range
               
        self.load_output_data(param= 'IonFlx')
        self.load_output_data(param= 'IonPol')
        
        d0sol_radflux = self.data['outputdata']['IonFlx']['D_0'][37, 1:97]
        print(d0sol_radflux.size)
        d1sol_radflux = self.data['outputdata']['IonFlx']['D_1'][37, 1:97]
        print(d1sol_radflux.size)
        radflux_dic = {'D_0': d0sol_radflux, 'D_1': d1sol_radflux}
        self.data['radflux'] = radflux_dic
        
        
        d0tar_polflux = self.data['outputdata']['IonPol']['D_0'][1:37, 1]
        print(d0tar_polflux.size)
        d1tar_polflux = self.data['outputdata']['IonPol']['D_1'][1:37, 96]
        print(d1tar_polflux.size)
        polflux_dic = {'D_0': d0tar_polflux, 'D_1': d1tar_polflux}
        self.data['polflux'] = polflux_dic
        
        
        "for D_0"
        
        
        
        
        # ln = len(pol_list)
        # pol_loc = np.zeros(ln)
        # i = 0
        
        # for ii in self.data['poloidal_index']:
        #     pol_loc[i] = int(ii)
        #     i = i + 1
    
    