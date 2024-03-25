# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 23:11:05 2024

@author: ychuang
"""

from SOLPSplotter_PRmap import RP_mapping
import opacity_plot_method as opm
import matplotlib.pyplot as plt
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate
from scipy.optimize import curve_fit
import numpy as np


class profile_fit(RP_mapping):
    
    def __init__(self, DefaultSettings, loadDS):
        RP_mapping.__init__(self, DefaultSettings, loadDS)
    
       
    
    def opacity_data_fit_method(self, b2fstate, Neuden, psiN, 
                psi_dsa_ratio, pol_list, itername): 
        # i = 0
        ln = len(pol_list)
        efold = np.zeros(ln)
        efold_l = np.zeros(ln)
        delta = np.zeros(ln)
        delta_l = np.zeros(ln)
        opq = np.zeros(ln)
        # pol_loc = np.zeros(ln)
        neu_den = np.zeros(ln)
        ne_ped = np.zeros(ln)
        tdelta = np.zeros(ln)
        fluxexp = np.zeros(ln)
        
        
        # self.load_output_data(param= 'NeuDen')
        # self.load_output_data(param= 'Ne')
        # self.load_output_data(param= 'Te')
        
        Ne_data = b2fstate['ne'].transpose()
        Te_J = b2fstate['te'].transpose()
        ev = 1.6021766339999999 * pow(10, -19)
        Te_data = Te_J / ev
        
        # self.data['simu_data_check'] = Te_data
        
        
        for k in pol_list:
            
            
            # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
            # SEP = int(self.data['DefaultSettings']['SEP'])
            
            pol_in = int(k)
            i = pol_list.index(k)
            psi = psiN[:, pol_in]
            
            
            Nd = Neuden[:, pol_in]
            Ne = Ne_data[:, pol_in]
            Te = Te_data[:, pol_in]
            

            rd = fm.Opacity_calculator(x_coord= psi, ne = Ne, te = Te, 
                                   neuden = Nd)
            ped_index = rd['sep_index']
            
            
            fit_dic = {'tanh_ne_fit': rd['tanh_ne_fit'], 'tanh_te_fit': rd['tanh_te_fit'],
                       'exp_fit': rd['exp_fit']}
            
            flux_expand = self.calc_flux_expansion(pol_loc = k, ped_index = ped_index, 
                                                   iter_index = itername)

            
            efold[i] = rd['efold_length']
            delta[i] = rd['pedestal_width']
            opq[i] = rd['dimensionless_opaqueness']
            neu_den[i] = rd['n_sep_fit']
            ne_ped[i] = rd['electron_pedestal_density']
            tdelta[i] = rd['temperature_pedestal_width']
            fluxexp[i] = flux_expand
            efold_l[i] = rd['efold_length']*psi_dsa_ratio*flux_expand
            delta_l[i] = rd['pedestal_width']*psi_dsa_ratio*flux_expand
            
            
            # pol_loc[i] = int(k)
            # i = i + 1
        
        result = {'efold_length_psiN': efold, 'pedestal_width_psiN': delta,
                  'dimensionless_opaqueness': opq, 
                  'neutral_density': neu_den, 
                  'electron_pedestal_density': ne_ped,
                  'temperature_pedestal_width': tdelta,
                  'flux_expansion': fluxexp,
                  'efold_length': efold_l, 'pedestal_width': delta_l,                            
                  }
        
        self.data['poloidal_itemname'] = list(result.keys())
        
        return result
    
    
    def opacity_data_fit(self, pol_list):
        
        # self.load_ft44()
        self.load_output_data(param= 'NeuDen')
        
        if self.withshift == False and self.withseries == False:
            
            for p in pol_list:
                self.calc_dsa(pol_loc= p)
            # data = self.data['ft44']['dab2']
            # Neuden_data = np.transpose(data[:, :, 0])
            Neuden_data = self.data['outputdata']['NeuDen']
            fstate = self.data['b2fstate']
            psiN_map = self.data['psi']['psival']
            pd = self.data['DefaultSettings']['psi_dsa']
            
            fitresult = self.opacity_data_fit_method(b2fstate = fstate, Neuden = Neuden_data, 
                       psiN = psiN_map, psi_dsa_ratio = pd, pol_list = pol_list, 
                                    itername = None)
            
            self.data['opacity_poloidal'] = fitresult
            self.data['poloidal_itemname'] = list(fitresult.keys())
        
        elif self.withshift == True and self.withseries == False:
            
            fitresult_dic = {}
            
            for p in pol_list:
                self.calc_dsa(pol_loc= p)
            
            for aa in self.data['dircomp']['multi_shift']:
                
                Neuden_data = self.data['outputdata']['NeuDen'][aa]
                # data = self.data['ft44'][aa]['dab2']
                # Neuden_data = np.transpose(data[:, :, 0])
                fstate = self.data['b2fstate'][aa]
                psiN_map = self.data['psi']['psival'][aa]
                pd = self.data['DefaultSettings']['psi_dsa'][aa]
                
                fitresult = self.opacity_data_fit_method(b2fstate = fstate, Neuden = Neuden_data,
                        psiN = psiN_map, psi_dsa_ratio = pd, pol_list = pol_list,
                                        itername = aa)
                
                fitresult_dic[aa] = fitresult
            
            self.data['opacity_poloidal'] = fitresult_dic
            self.data['poloidal_itemname'] = list(fitresult_dic['org'].keys())
        
        elif self.withshift == False and self.withseries == True:
            
            fitresult_dic = {}
            
            for p in pol_list:
                self.calc_dsa(pol_loc= p)
            
            for aa in list(self.data['dircomp']['Attempt'].keys()):
                
                Neuden_data = self.data['outputdata']['NeuDen'][aa]
                # data = self.data['ft44'][aa]['dab2']
                # Neuden_data = np.transpose(data[:, :, 0])
                fstate = self.data['b2fstate'][aa]
                psiN_map = self.data['psi']['psival']
                pd = self.data['DefaultSettings']['psi_dsa']
                
                fitresult = self.opacity_data_fit_method(b2fstate = fstate, Neuden = Neuden_data,
                            psiN = psiN_map, psi_dsa_ratio = pd, pol_list = pol_list,
                                        itername = None)
                
                fitresult_dic[aa] = fitresult
            
            self.data['opacity_poloidal'] = fitresult_dic
            self.data['poloidal_itemname'] = list(fitresult.keys())
        
        elif self.withshift == True and self.withseries == True:
            print('opacity_data_fit is not there yet!')
        
        
        else:
            print('opacity_data_fit has a bug')
                
                
    def radial_data_fit_method(self, b2fstate, Neuden, psiN, pol_list):
        
        # self.load_output_data(param= 'NeuDen')
        # self.load_output_data(param= 'Ne')
        # self.load_output_data(param= 'Te')
        
        Ne_data = b2fstate['ne'].transpose()
        Te_J = b2fstate['te'].transpose()
        ev = 1.6021766339999999 * pow(10, -19)
        Te_data = Te_J / ev
        pol_in = int(pol_list[0])
        
        psi = psiN[:, pol_in]
        
        
        Nd = Neuden[:, pol_in]
        Ne = Ne_data[:, pol_in]
        Te = Te_data[:, pol_in]
        

        rd = fm.Opacity_calculator(x_coord= psi, ne = Ne, te = Te, 
                               neuden = Nd)
             
        fit_dic = {'tanh_ne_fit': rd['tanh_ne_fit'], 'tanh_te_fit': rd['tanh_te_fit'],
             'exp_fit': rd['exp_fit'], 'x_coord_cut': rd['x_coord_cut'],
                     'ne_symmetry_point': rd['ne_symmetry_point'], 
                      'te_symmetry_point': rd['te_symmetry_point'],
               'n_sep_fit': rd['n_sep_fit'], 'NeuDen': Nd, 'Ne': Ne, 'Te': Te}
        
        return fit_dic
    
    
    def radial_data_fit(self, pol_list):
        
        # self.load_output_data(param= 'NeuDen')
        
        self.load_ft44()
        
        if self.withshift == False and self.withseries == False:
            
            data = self.data['ft44']['dab2']
            Neuden_data = np.transpose(data[:, :, 0])
            # Neuden_data = self.data['outputdata']['NeuDen']
            fstate = self.data['b2fstate']
            psiN_map = self.data['psi']['psival']
            
            fitresult = self.radial_data_fit_method(b2fstate = fstate, 
                        Neuden = Neuden_data, psiN = psiN_map, pol_list = pol_list)
            
            self.data['radial_fit_data'] = fitresult
        
        elif self.withshift == True and self.withseries == False:
            
            fitresult_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                
                Neuden_data = self.data['outputdata']['NeuDen'][aa]
                fstate = self.data['b2fstate'][aa]
                psiN_map = self.data['psi']['psival'][aa]
                
                fitresult = self.radial_data_fit_method(b2fstate = fstate, 
                            Neuden = Neuden_data, psiN = psiN_map, pol_list = pol_list)
                
                fitresult_dic[aa] = fitresult
            
            self.data['radial_fit_data'] = fitresult_dic
        
        elif self.withshift == False and self.withseries == True:
            
            fitresult_dic = {}
            
            for aa in list(self.data['dircomp']['Attempt'].keys()):
                
                Neuden_data = self.data['outputdata']['NeuDen'][aa]
                fstate = self.data['b2fstate'][aa]
                psiN_map = self.data['psi']['psival']
                
                fitresult = self.radial_data_fit_method(b2fstate = fstate, 
                            Neuden = Neuden_data, psiN = psiN_map, pol_list = pol_list)
                
                fitresult_dic[aa] = fitresult
            
            self.data['radial_fit_data'] = fitresult_dic
        
        elif self.withshift == True and self.withseries == True:
            print('radial_data_fit is not there yet!')
        
        
        else:
            print('radial_data_fit has a bug')





            
# ----------------------------------------------------------------------------           
    
    """
    backup:
        

    
    def opacity_data_method_single(self, pol_list): 
        i = 0
        ln = len(pol_list)
        efold = np.zeros(ln)
        efold_l = np.zeros(ln)
        delta = np.zeros(ln)
        delta_l = np.zeros(ln)
        opq = np.zeros(ln)
        # pol_loc = np.zeros(ln)
        neu_den = np.zeros(ln)
        ne_ped = np.zeros(ln)
        tdelta = np.zeros(ln)
        fluxexp = np.zeros(ln)
        
        
        for k in pol_list:
            psi = self.data['psi']['psi_{}_val'.format(k)][:, 2]
            
            # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
            # SEP = int(self.data['DefaultSettings']['SEP'])
            
            pol_in = int(k)
            Nd = self.data['outputdata']['NeuDen'][:, pol_in]
            Ne = self.data['outputdata']['Ne'][:, pol_in]
            Te = self.data['outputdata']['Te'][:, pol_in]
            

            rd = fm.Opacity_calculator(x_coord= psi, ne = Ne, te = Te, 
                                   neuden = Nd)
            ped_index = rd['sep_index']
            
            fe = self.calc_flux_expansion(pol_loc= k, 
                            ped_index= ped_index, iter_index= None)
            pd = self.data['DefaultSettings']['psi_dsa']
            

            efold[i] = rd['efold_length']
            delta[i] = rd['pedestal_width']
            opq[i] = rd['dimensionless_opaqueness']
            neu_den[i] = rd['n_sep_fit']
            ne_ped[i] = rd['electron_pedestal_density']
            tdelta[i] = rd['temperature_pedestal_width']
            fluxexp[i] = fe
            efold_l[i] = rd['efold_length']*pd*fe
            delta_l[i] = rd['pedestal_width']*pd*fe
            
            
            # pol_loc[i] = int(k)
            i = i + 1
        
        result = {'efold_length_psiN': efold, 'pedestal_width_psiN': delta,
                  'dimensionless_opaqueness': opq, 
                  'neutral_density': neu_den, 
                  'electron_pedestal_density': ne_ped,
                  'temperature_pedestal_width': tdelta,
                  'flux_expansion': fluxexp,
                  'efold_length': efold_l, 'pedestal_width': delta_l,
                                  
                  }
        
        return result

    def opacity_data_method_multi(self, pol_list, iter_list):
        efold_dic = {}
        efold_leng_dic = {}
        delta_dic = {}
        delta_leng_dic = {}
        opq_dic = {}
        neu_den_dic = {}
        ne_ped_dic = {}
        tdelta_dic = {}
        flux_expand_dic = {}
        xcoord_cut_dic = {}
        
        for aa in iter_list:
            i = 0
            ln = len(pol_list)
            efold = np.zeros(ln)
            efold_l = np.zeros(ln)
            delta = np.zeros(ln)
            delta_l = np.zeros(ln)
            opq = np.zeros(ln)
            neu_den = np.zeros(ln)
            ne_ped = np.zeros(ln)
            tdelta = np.zeros(ln)
            flux_exp = np.zeros(ln)
            xcoord_cut = {}
        
            
            
            for k in pol_list:
                
                if self.withshift == True and self.withseries == False:
                    psi = self.data['psi']['psi_{}_val'.format(k)][aa][:, 2]
                elif self.withshift == False and self.withseries == True:
                    psi = self.data['psi']['psi_{}_val'.format(k)][:, 2]
                else:
                    print('out of expectation')
                
                
                # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
                # SEP = int(self.data['DefaultSettings']['SEP'])
                
                pol_in = int(k)
                Nd = self.data['outputdata']['NeuDen'][aa][:, pol_in]
                Ne = self.data['outputdata']['Ne'][aa][:, pol_in]
                Te = self.data['outputdata']['Te'][aa][:, pol_in]
                
                rd = fm.Opacity_calculator(x_coord = psi, ne = Ne, te = Te,
                                           neuden = Nd)
                
                xcoord_cut_index = {k: rd['x_coord_cut']}
                xcoord_cut = xcoord_cut | xcoord_cut_index
                
                ped_index = rd['sep_index']
                
                fe = self.calc_flux_expansion(pol_loc= k, 
                                ped_index= ped_index, iter_index= aa)
                
                if self.withshift == True and self.withseries == False:
                    pd = self.data['DefaultSettings']['psi_dsa'][aa]
                elif self.withshift == False and self.withseries == True:
                    pd = self.data['DefaultSettings']['psi_dsa']
                else:
                    print('there is a bug')

                
                
                efold[i] = rd['efold_length']
                efold_l[i] = rd['efold_length']*fe*pd
                delta[i] = rd['pedestal_width']
                delta_l[i] = rd['pedestal_width']*fe*pd
                opq[i] = rd['dimensionless_opaqueness']
                neu_den[i] = rd['n_sep_fit']
                ne_ped[i] = rd['electron_pedestal_density']
                tdelta[i] = rd['temperature_pedestal_width']
                flux_exp[i] = fe

                # pol_loc[i] = int(k)
                i = i + 1
            
            efold_dic[aa] = efold
            efold_leng_dic[aa] = efold_l
            delta_dic[aa] = delta
            delta_leng_dic[aa] = delta_l
            opq_dic[aa] = opq
            neu_den_dic[aa] = neu_den
            ne_ped_dic[aa] = ne_ped
            tdelta_dic[aa] = tdelta
            flux_expand_dic[aa] = flux_exp
            xcoord_cut_dic[aa] = xcoord_cut
            
            
        
        result = {'efold_length_psiN': efold_dic, 
                  'pedestal_width_psiN': delta_dic,
                  'efold_length': efold_leng_dic, 
                  'pedestal_width': delta_leng_dic,
                  'dimensionless_opaqueness': opq_dic, 
                  'neutral_density': neu_den_dic, 
                  'electron_pedestal_density': ne_ped_dic,
                  'temperature_pedestal_width': tdelta_dic,
                  'flux_expansion': flux_expand_dic               
                  }
        self.data['xcoord_cut'] = xcoord_cut_dic
        
        return result

    """