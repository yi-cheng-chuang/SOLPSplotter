# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 23:22:49 2025

@author: ychuang
"""

from SOLPS_input.Amal_dir_input import Amal_directory_input
from load_directory.SOLPS_load_directory import load_directory
from load_directory.load_dirdata_method import load_dir_method
from load_directory.grab_attempt_number import grab_aptn_method
from load_coordinate.load_coord_method import load_coordgeo_method
from load_coordinate.SOLPSplotter_geo import load_geometry
from set_plot.plot_format import plot_setting
from fit_data.fitting_method import fit_method_collection
from midplane_data.SOLPSplotter_PRmap import RP_mapping
from load_experimental_data.load_expdata_method import read_expdata_method
from load_experimental_data.SOLPSplotter_load_expdata import load_expdata
from load_experimental_data.SOLPSplotter_mastuTS import preprocess_mastuTS
from load_experimental_data.Amal_mastu_TS import Amal_Steven_mastuTS
from load_simulation_data.read_B2simulation_data import load_B2simu_data
from load_simulation_data.read_ft_files import load_ftfiles_data
from fit_data.SOLPSplotter_fit import profile_fit
from midplane_data.measure_minor_radius_plot import minor_radius_measurement
from load_experimental_data.fitTS_usinglmfit import lmfit_TS




class mastu_load_prepare_module:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
      

    def mastu_prepare_for_plots(self):
        
        xdi = Amal_directory_input(DF = self.DF, data = self.data)
        xga = grab_aptn_method(DF = self.DF, data = self.data)
        xldm = load_dir_method(DF = self.DF, data = self.data, di= xdi, gam = xga)
        xld = load_directory(DF = self.DF, data = self.data, di= xdi, ldm= xldm)
        xlc = load_coordgeo_method(DF = self.DF, data = self.data)
        xlg = load_geometry(DF = self.DF, data = self.data, lcm = xlc)
        xps = plot_setting(DF = self.DF, data = self.data)
        xfm = fit_method_collection(DF = self.DF, data = self.data)
        xrp = RP_mapping(DF = self.DF, data = self.data, lcm = xlc, fmc= xfm)
        xre = read_expdata_method(DF = self.DF, data = self.data)
        xle = load_expdata(DF = self.DF, data = self.data, fmc= xfm, rem = xre, lg= xlg)
        xpm = preprocess_mastuTS(DF = self.DF, data = self.data, fmc = xfm)
        xlb = load_B2simu_data(DF = self.DF, data = self.data, ldm= xldm)
        xlf = load_ftfiles_data(DF = self.DF, data = self.data, ldm= xldm)
        xpf = profile_fit(DF = self.DF, data = self.data, fmc = xfm, rp= xrp)
        xmrm = minor_radius_measurement(DF = self.DF, data = self.data, rp= xrp)
        xlt = lmfit_TS(DF = self.DF, data = self.data, fmc = xfm, pm = xpm)
        xasmt = Amal_Steven_mastuTS(DF = self.DF, data = self.data, fmc = xfm)

        
        
        
        
        
        if self.DF.DEV == 'mastu':
            
            xld.load_mastu_dir()
            xlg.load_mastusolpsgeo()
            # xlg.load_vessel()
            xlg.calcpsi_avcr()
            print('working on RRsep TSmeasure now...')
            xrp.calc_RRsep(plotRR= True, plot_psi_dsa_align= False, midplane_loc = 'TSmeasure')
            print('working on RRsep maxis now...')
            xrp.calc_RRsep(plotRR= True, plot_psi_dsa_align= False, midplane_loc = 'maxis')
            xmrm.calc_minor_radius(plotRR= True)          
            xasmt.Amal_load_mastu_TS(plot_OD = False, plot_P = True, writefile = True)
            
            
            
            
            
            
            