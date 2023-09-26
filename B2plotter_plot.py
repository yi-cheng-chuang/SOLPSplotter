# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 11:47:48 2023

@author: user
"""
from B2plotter_class import load_data
import opacity_plot_method as opm
import matplotlib.pyplot as plt
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate
from scipy.optimize import curve_fit
import numpy as np
import xarray as xr


class Opacity_study(load_data):
    def __init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters, Publish):
        load_data.__init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters)
        
        self.Publish = Publish
        self.data['DefaultSettings']['Publish'] = self.Publish
    
    
    def set_plot(self):
        if self.Publish == 'b2plottersetting':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth=5, markersize=9)
            plt.rcParams.update({'font.size': 20})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
            # font = {‘family’ : ‘monospace’,
            #         ‘weight’ : ‘italic’,
            #         ‘size’ : ‘medium’}
            # # pass in the font dict as kwargs
            # rc(‘font’, **font)
    
 
    
        else:
            print('Publish setting is incorrect or add another setting')
    

    def Opacity_study_poloidal_plot(self, pol_list, x_choice):
        self.data['poloidal_index'] = pol_list
        
        for j in pol_list:
            self.calcpsi_1D(pol_loc= j)
        
        self.load_output_data(param= 'NeuDen')
        self.load_output_data(param= 'Ne')
        self.load_output_data(param= 'Te')
        
        
        
        ln = len(pol_list)
        pol_loc = np.zeros(ln)
        i = 0
        
        for ii in self.data['poloidal_index']:
            pol_loc[i] = int(ii)
            i = i + 1
        
        if self.withshift == False and self.withseries == False:
            result = self.opacity_data_method_single(pol_list = pol_list, 
                                                     x_choice= x_choice)
            
            self.data['opacity_poloidal'] = result
            
            unit = {'efold_length': 'efold length: $\lambda_{n_D}$: [m]',
                    'pedestal_width': 'Pedestal width: $\Delta n$: [m]',
                      'dimensionless_opaqueness': 'dimensionless opaqueness', 
                      'neutral_density': 'neutral density ${n_D}$ (m$^{-3}$)', 
                      'electron_pedestal_density': 'electron pedestal density: $n_{ped}$ (m$^{-3}$)',
                      'temperature_pedestal_width': 'temperature pedestal width: $\Delta T$: [m]',
                      
                      'efold_length_method2': 'efold length: $\lambda_{n_D}$: [m]',
                      'dimensionless_opaqueness_method2': 'dimensionless opaqueness',
                      
                      'efold_length_method3': 'efold length: $\lambda_{n_D}$: [m]',
                      'dimensionless_opaqueness_method3': 'dimensionless opaqueness',
                      }
            
            char = {}
            char['withshift'] = self.withshift
            char['withseries'] = self.withseries
            
            log_flag = False
            
            # print(result.keys())
            # print(unit.keys())
            
            opm.opacity_plot(pol_loc = pol_loc, result_dic = result, unit_dic = unit,
                             log_flag = False, charactor= char,
                             iter_list = None, 
                             change_ver_dic = None)
            
        elif self.withshift == True and self.withseries == False:
            result = self.opacity_data_method_multi(pol_list = pol_list, 
                                iter_list = self.data['dircomp']['multi_shift'])
            
            self.data['opacity_poloidal'] = result
            
            ll = len(self.data['dircomp']['Attempt'].keys())
            mm = len(pol_list)
            
            series_list = self.data['dircomp']['Attempt'].keys()
            
            # data_collect_opq = xr.DataArray(np.zeros((ll, mm)), 
            #                       coords=[series_list, pol_list], 
            #                 dims=['different_density','Poloidal_Location'], 
            #                  name = r'dimensionless opaqueness $m$')
            
            
            
            data_collect_opq = np.zeros((mm, ll))
            i = 0
            for la in self.data['dircomp']['multi_shift']:
                lb = np.asarray(result['dimensionless_opaqueness_method2'][la])
                data_collect_opq[:, i] = lb
                i = i + 1
            
            self.data['data_collect'] = data_collect_opq
            
            shift_list = np.zeros(ll)
            ka = 0
            for k in self.data['dircomp']['multi_shift']:
                shift_list[ka] = float(self.data['dircomp']['shift_dic'][k])
                ka = ka + 1

            char = {}
            char['withshift'] = self.withshift
            char['withseries'] = self.withseries

            opm.data_reorder(iter_list = pol_list, change_var = shift_list,
                             data_collect = data_collect_opq, char = char)
            
            
            
            unit = {'efold_length': 'efold length: $\lambda_{n_D}$: [m]', 
                    'pedestal_width': 'Pedestal width: $\Delta n$: [m]',
                      'dimensionless_opaqueness': 'dimensionless opaqueness', 
                      'neutral_density': 'neutral density ${n_D}$ (m$^{-3}$)', 
                      'electron_pedestal_density': 'electron pedestal density: $n_{ped}$ (m$^{-3}$)',
                      'temperature_pedestal_width': 'temperature pedestal width: $\Delta T$: [m]',
                      
                      'efold_length_method2': 'efold length: $\lambda_{n_D}$: [m]',
                      'dimensionless_opaqueness_method2': 'dimensionless opaqueness',
                      
                      # 'efold_length_method3': 'efold length: $\lambda_{n_D}$: [m]',
                      # 'dimensionless_opaqueness_method3': 'dimensionless opaqueness',
                      }
            
            # char = {}
            # char['withshift'] = self.withshift
            # char['withseries'] = self.withseries
            
            # log_flag = False
            
            shift_dic = {}
            for k in self.data['dircomp']['multi_shift']:
                p = str(self.data['dircomp']['shift_dic'][k])
                shift_dic[k] = p
            
            opm.opacity_plot(pol_loc = pol_loc, result_dic = result, unit_dic = unit,
                             log_flag = False, charactor= char,
                             iter_list = self.data['dircomp']['multi_shift'], 
                             change_ver_dic = shift_dic)
        
        elif self.withshift == False and self.withseries == True:
            result = self.opacity_data_method_multi(pol_list = pol_list, 
                                iter_list = self.data['dircomp']['Attempt'].keys())
            
            self.data['opacity_poloidal'] = result
            
            ll = len(self.data['dircomp']['Attempt'].keys())
            mm = len(pol_list)
            
            series_list = self.data['dircomp']['Attempt'].keys()
            
            # data_collect_opq = xr.DataArray(np.zeros((ll, mm)), 
            #                       coords=[series_list, pol_list], 
            #                 dims=['different_density','Poloidal_Location'], 
            #                  name = r'dimensionless opaqueness $m$')
            
            
            
            data_collect_opq = np.zeros((mm, ll))
            i = 0
            for la in self.data['dircomp']['Attempt'].keys():
                lb = np.asarray(result['dimensionless_opaqueness_method2'][la])
                data_collect_opq[:, i] = lb
                i = i + 1
            
            self.data['data_collect'] = data_collect_opq
            
            density_list = np.zeros(ll)
            ka = 0
            for k in self.data['dircomp']['Attempt'].keys():
                density_list[ka] = float(k)
                ka = ka + 1

            char = {}
            char['withshift'] = self.withshift
            char['withseries'] = self.withseries

            opm.data_reorder(iter_list = pol_list, change_var = density_list,
                             data_collect = data_collect_opq, char = char)
            
            
            unit = {'efold_length': 'efold length: $\lambda_{n_D}$: [m]', 
                    'pedestal_width': 'Pedestal width: $\Delta n$: [m]',
                      'dimensionless_opaqueness': 'dimensionless opaqueness', 
                      'neutral_density': 'neutral density ${n_D}$ (m$^{-3}$)', 
                      'electron_pedestal_density': 'electron pedestal density: $n_{ped}$ (m$^{-3}$)',
                      'temperature_pedestal_width': 'temperature pedestal width: $\Delta T$: [m]',
                      
                      'efold_length_method2': 'efold length: $\lambda_{n_D}$: [m]',
                      'dimensionless_opaqueness_method2': 'dimensionless opaqueness',
                      
                      'efold_length_method3': 'efold length: $\lambda_{n_D}$: [m]',
                      'dimensionless_opaqueness_method3': 'dimensionless opaqueness',
                      }
            
            
            
            # log_flag = False
            
            density_dic = {}
            for k in self.data['dircomp']['Attempt'].keys():
                kk = float(k)*pow(10, 19)
                density_dic[k] = kk
            
            opm.opacity_plot(pol_loc = pol_loc, result_dic = result, unit_dic = unit,
                             log_flag = False, charactor= char,
                             iter_list = self.data['dircomp']['Attempt'].keys(), 
                             change_ver_dic = density_dic)
        
        elif self.withshift == True and self.withseries == True:
            print('Opacity_study_poloidal_plot is not there yet, to be continue...')
            
            
        else:
            print('more work need to be done')
            
            
    def opacity_data_method_single(self, pol_list, x_choice): 
        i = 0
        ln = len(pol_list)
        efold = np.zeros(ln)
        delta = np.zeros(ln)
        opq = np.zeros(ln)
        # pol_loc = np.zeros(ln)
        neu_den = np.zeros(ln)
        ne_ped = np.zeros(ln)
        tdelta = np.zeros(ln)
        
        efold_m2 = np.zeros(ln)
        opq_m2 = np.zeros(ln)
        std_m2 = np.zeros(ln)
        
        efold_m3 = np.zeros(ln)
        opq_m3 = np.zeros(ln)
        
        
        for k in pol_list:
            psi = self.data['psi']['psi_{}_val'.format(k)][:, 2]
            
            # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
            # SEP = int(self.data['DefaultSettings']['SEP'])
            
            pol_in = int(k) + 1
            Nd = self.data['outputdata']['NeuDen'][:, pol_in]
            Ne = self.data['outputdata']['Ne'][:, pol_in]
            Te = self.data['outputdata']['Te'][:, pol_in]
            
            if x_choice == 'psiN':
                rd = fm.Opacity_calculator(x_choice= 'psiN', x_coord= psi, ne = Ne, te = Te, 
                                       neuden = Nd)
                
                # Opacity_calculator(x_choice, x_coord, ne, te, neuden)
                
            elif x_choice == 'RRsep':
                print('RRsep version has not been built...')
                # rd = fm.Opacity_calculator(x_choice= 'RRsep', x_coord = dsa_pol_loc, ne = Ne, te = Te, 
                #                        neuden = Nd, psi = psi)
                
            
            
            
            efold[i] = rd['efold_length']
            delta[i] = rd['pedestal_width']
            opq[i] = rd['dimensionless_opaqueness']
            neu_den[i] = max(rd['exp_fit'])
            ne_ped[i] = rd['electron_pedestal_density']
            tdelta[i] = rd['temperature_pedestal_width']
            
            efold_m2[i] = rd['efold_length_method2']
            opq_m2[i] = rd['dimensionless_opaqueness_method2']
            std_m2[i] = rd['std_m2']
            
            efold_m3[i] = rd['efold_length_method3']
            opq_m3[i] = rd['dimensionless_opaqueness_method3']
            # pol_loc[i] = int(k)
            i = i + 1
        
        result = {'efold_length': efold, 'pedestal_width': delta,
                  'dimensionless_opaqueness': opq, 
                  'neutral_density': neu_den, 
                  'electron_pedestal_density': ne_ped,
                  'temperature_pedestal_width': tdelta,
                  
                  'efold_length_method2': efold_m2, 
                  'dimensionless_opaqueness_method2': opq_m2,
                  'std_m2': std_m2,
                  
                  # 'efold_length_method3': efold_m3, 
                  # 'dimensionless_opaqueness_method3': opq_m3,
                  
                  
                  
                  
                  }
        
        return result

    def opacity_data_method_multi(self, pol_list, iter_list):
        efold_dic = {}
        delta_dic = {}
        opq_dic = {}
        neu_den_dic = {}
        ne_ped_dic = {}
        tdelta_dic = {}
        
        efold_m2_dic = {}
        opq_m2_dic = {}
        std_m2_dic = {}
        
        efold_m3_dic = {}
        opq_m3_dic = {}
        
        for aa in iter_list:
            i = 0
            ln = len(pol_list)
            efold = np.zeros(ln)
            delta = np.zeros(ln)
            opq = np.zeros(ln)
            # pol_loc = np.zeros(ln)
            neu_den = np.zeros(ln)
            ne_ped = np.zeros(ln)
            tdelta = np.zeros(ln)
            
            efold_m2 = np.zeros(ln)
            opq_m2 = np.zeros(ln)
            std_m2 = np.zeros(ln)
            
            efold_m3 = np.zeros(ln)
            opq_m3 = np.zeros(ln)
            
            
            for k in pol_list:
                
                if self.withshift == True and self.withseries == False:
                    psi = self.data['psi']['psi_{}_val'.format(k)][aa]
                    # dsa_pol_loc = self.data['dsa']['dsa_{}'.format(k)][aa]['dsa_{}_val'.format(k)]
                elif self.withshift == False and self.withseries == True:
                    psi = self.data['psi']['psi_{}_val'.format(k)]
                    dsa_pol_loc = self.data['dsa']['dsa_{}'.format(k)]['dsa_{}_val'.format(k)]
                else:
                    print('out of expectation')
                
                
                # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
                # SEP = int(self.data['DefaultSettings']['SEP'])
                
                pol_in = int(k) + 1
                Nd = self.data['outputdata']['NeuDen'][aa][:, pol_in]
                Ne = self.data['outputdata']['Ne'][aa][:, pol_in]
                Te = self.data['outputdata']['Te'][aa][:, pol_in]
                
                rd = fm.Opacity_calculator(x_choice= 'psiN', x_coord = psi, ne = Ne, te = Te,
                                           neuden = Nd)
                
                # Opacity_calculator(x_choice, x_coord, ne, te, neuden)
                
                
                efold[i] = rd['efold_length']
                delta[i] = rd['pedestal_width']
                opq[i] = rd['dimensionless_opaqueness']
                neu_den[i] = max(rd['exp_fit'])
                ne_ped[i] = rd['electron_pedestal_density']
                tdelta[i] = rd['temperature_pedestal_width']
                efold_m2[i] = rd['efold_length_method2']
                opq_m2[i] = rd['dimensionless_opaqueness_method2']
                std_m2[i] = rd['std_m2']
                
                efold_m3[i] = rd['efold_length_method3']
                opq_m3[i] = rd['dimensionless_opaqueness_method3']
                # pol_loc[i] = int(k)
                i = i + 1
            
            efold_dic[aa] = efold
            delta_dic[aa] = delta
            opq_dic[aa] = opq
            neu_den_dic[aa] = neu_den
            ne_ped_dic[aa] = ne_ped
            tdelta_dic[aa] = tdelta
            
            efold_m2_dic[aa] = efold_m2
            opq_m2_dic[aa] = opq_m2
            std_m2_dic[aa] = std_m2
            
            efold_m3_dic[aa] = efold_m3
            opq_m3_dic[aa] = opq_m3
        
        result = {'efold_length': efold_dic, 'pedestal_width': delta_dic,
                  'dimensionless_opaqueness': opq_dic, 
                  'neutral_density': neu_den_dic, 
                  'electron_pedestal_density': ne_ped_dic,
                  'temperature_pedestal_width': tdelta_dic,
                  
                  'efold_length_method2': efold_m2_dic, 
                  'dimensionless_opaqueness_method2': opq_m2_dic,
                  'std_m2': std_m2_dic,
                  
                  'efold_length_method3': efold_m3_dic, 
                  'dimensionless_opaqueness_method3': opq_m3_dic,
                  }
        
        return result
          
            
    def Opacity_study_radial_plot(self, pol_loc, x_choice):
        self.load_output_data(param= 'NeuDen')
        self.load_output_data(param= 'Ne')
        self.load_output_data(param= 'Te')
        
        if self.withshift == False and self.withseries == False:
          
            psi = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 2]
            # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
            SEP = int(self.data['DefaultSettings']['SEP'])
            
            pol_index = int(pol_loc) + 1
            Nd = self.data['outputdata']['NeuDen'][:, pol_index]
            Ne = self.data['outputdata']['Ne'][:, pol_index]
            Te = self.data['outputdata']['Te'][:, pol_index]
            
            
            # result_dic = fm.Opacity_calculator(dsa_pol_loc, Ne, Te, Nd, psi)
            result_dic = fm.Opacity_calculator(x_choice = x_choice, x_coord = psi, 
                                           ne = Ne, te = Te, neuden = Nd)
            
            self.data['opacity_study'] = result_dic
            
            P = self.data['Parameter']
            opm.opacity_radial_method_single(result_dic, SEP, x_choice = x_choice, 
                                             x_coord = psi, Nd = Nd, Ne = Ne, Te = Te, 
                                                  P = P, log_flag = True)
        
        elif self.withshift == True and self.withseries == False:
            self.load_output_data(param= 'NeuDen')
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
            
            efold_dic = {}
            delta_dic = {}
            tdelta_dic = {}
            opq_dic = {}
            xcoord_cut_dic = {}
            # xcoord_dic = {}
            sym_pt_dic = {}
            te_sym_pt_dic = {}
            
            
            psi_dic = {}
            dsa_pol_loc_dic = {}
            
            # psi_psn_dic = {}    
            # dsa_psn_dic = {}
            
            
            Nd_dic = {}
            Ne_dic = {}
            Te_dic = {}
            exp_an_fit_dic = {}
            tanh_ne_fit_dic = {}
            tanh_te_fit_dic = {}
            
            exp_fit_m2_dic = {}
            x_m2_dic = {}
            efold_m2_dic = {}
            opq_m2_dic = {}
            
            exp_fit_m3_dic = {}
            x_m3_dic = {}
            efold_m3_dic = {}
            opq_m3_dic = {}
            
            
            for aa in self.data['dircomp']['multi_shift']:
                psi_dic[aa] = self.data['psi']['psi_{}_val'.format(pol_loc)][aa]
                # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
                
                # if x_choice == 'psiN':
                #     xcoord_dic[aa] = psi_dic[aa]
                
                pol_index = int(pol_loc) + 1
                Nd_dic[aa] = self.data['outputdata']['NeuDen'][aa][:, pol_index]
                Ne_dic[aa] = self.data['outputdata']['Ne'][aa][:, pol_index]
                Te_dic[aa] = self.data['outputdata']['Te'][aa][:, pol_index]
                
                
                result_dic = fm.Opacity_calculator(x_choice = x_choice, x_coord = psi_dic[aa], 
                            ne = Ne_dic[aa], te = Te_dic[aa], neuden = Nd_dic[aa])
                
                
                
                
                tanh_ne_fit_dic[aa] = result_dic['tanh_ne_fit']
                tanh_te_fit_dic[aa] = result_dic['tanh_te_fit']
                exp_an_fit_dic[aa] = result_dic['exp_fit']
                delta_dic[aa] = result_dic['pedestal_width']
                tdelta_dic[aa] = result_dic['temperature_pedestal_width']
                efold_dic[aa] = result_dic['efold_length']
                opq_dic[aa] = result_dic['dimensionless_opaqueness']
                xcoord_cut_dic[aa] = result_dic['x_coord_cut']
                sym_pt_dic[aa] = result_dic['ne_symmetry_point']
                te_sym_pt_dic[aa] = result_dic['te_symmetry_point']
                
                exp_fit_m2_dic[aa] = result_dic['exp_fit_m2']
                x_m2_dic[aa] = result_dic['x_m2']
                efold_m2_dic[aa] = result_dic['efold_length_method2']
                opq_m2_dic[aa] = result_dic['dimensionless_opaqueness_method2']
                
                exp_fit_m3_dic[aa] = result_dic['exp_fit_m3']
                x_m3_dic[aa] = result_dic['x_m3']
                efold_m3_dic[aa] = result_dic['efold_length_method3']
                opq_m3_dic[aa] = result_dic['dimensionless_opaqueness_method3']
                
                # dsa_psn_dic[aa] = result_dic['dsa_psn']
                # psi_psn_dic[aa] = result_dic['psi_psn']
                
            
            
            # result = {'efold_length': efold_dic, 'pedestal_width': delta_dic,
            #           'dimensionless_opaqueness': opq_dic, 
            #           'temperature_pedestal_width': tdelta_dic,
            #           'tanh_ne_fit': tanh_ne_fit_dic,
            #           'tanh_te_fit': tanh_te_fit_dic,
            #           'exp_fit': exp_an_fit_dic,
            #           'xcoord_cut': xcoord_cut_dic,
            #           'neutral_density': Nd_dic,
            #           'electron_density': Ne_dic,
            #           'electron_temperature_density': Te_dic,
            #           'psiN': psi_dic,'dsa': dsa_pol_loc_dic, 
            #           'psi_psn': psi_psn_dic, 'dsa_psn': dsa_psn_dic
            #           }
            
            
            result = {'efold_length': efold_dic, 'pedestal_width': delta_dic,
                      'dimensionless_opaqueness': opq_dic, 
                      'temperature_pedestal_width': tdelta_dic,
                      'tanh_ne_fit': tanh_ne_fit_dic,
                      'tanh_te_fit': tanh_te_fit_dic,
                      'exp_fit': exp_an_fit_dic,
                      'xcoord_cut': xcoord_cut_dic,
                      'neutral_density': Nd_dic,
                      'electron_density': Ne_dic,
                      'electron_temperature_density': Te_dic,
                      'psiN': psi_dic,'dsa': dsa_pol_loc_dic,
                      'te_symmetry_point': te_sym_pt_dic,
                      'ne_symmetry_point': sym_pt_dic,
                      
                      'exp_fit_m2': exp_fit_m2_dic,
                      'x_m2': x_m2_dic,
                      'efold_length_method2': efold_m2_dic,
                      'dimensionless_opaqueness_method2': opq_m2_dic,
                      
                      'exp_fit_m3': exp_fit_m3_dic,
                      'x_m3': x_m3_dic,
                      'efold_length_method3': efold_m3_dic,
                      'dimensionless_opaqueness_method3': opq_m3_dic
                      }
            
            
            
            rl = self.data['dircomp']['multi_shift']
            SEP = self.data['DefaultSettings']['SEP']
            
            self.data['opacity_study'] = result
            
            char = {}
            char['withshift'] = self.withshift
            char['withseries'] = self.withseries
            
            shift_dic = {}
            for k in self.data['dircomp']['multi_shift']:
                p = str(self.data['dircomp']['shift_dic'][k])
                shift_dic[k] = p
            
            log_flag = True
            ii = 1
            
            P = self.data['Parameter']
            opm.opacity_radial_method_multi(x_choice = x_choice, result_dic = result, SEP = SEP, 
                            iter_list = self.data['dircomp']['multi_shift'], 
                             change_var_dic = shift_dic, log_flag = True, 
                             char = char, P = P)
        
            
        elif self.withshift == False and self.withseries == True:
            self.load_output_data(param= 'NeuDen')
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
            
            efold_dic = {}
            delta_dic = {}
            tdelta_dic = {}
            opq_dic = {}
            psi_dic = {}
            dsa_pol_loc_dic = {}
            xcoord_cut_dic = {}
            
            sym_pt_dic = {}
            te_sym_pt_dic = {}
            
            # psi_psn_dic = {}  
            # dsa_psn_dic = {}
            
            Nd_dic = {}
            Ne_dic = {}
            Te_dic = {}
            exp_an_fit_dic = {}
            tanh_ne_fit_dic = {}
            tanh_te_fit_dic = {}
            
            exp_fit_m2_dic = {}
            x_m2_dic = {}
            efold_m2_dic = {}
            opq_m2_dic = {}
            
            exp_fit_m3_dic = {}
            x_m3_dic = {}
            efold_m3_dic = {}
            opq_m3_dic = {}
            
            
            for aa in self.data['dircomp']['Attempt'].keys():
                psi_dic[aa] = self.data['psi']['psi_{}_val'.format(pol_loc)]
                dsa_pol_loc_dic[aa] = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]
                # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
                
                
                pol_index = int(pol_loc) + 1
                Nd_dic[aa] = self.data['outputdata']['NeuDen'][aa][:, pol_index]
                Ne_dic[aa] = self.data['outputdata']['Ne'][aa][:, pol_index]
                Te_dic[aa] = self.data['outputdata']['Te'][aa][:, pol_index]
                
                
                result_dic = fm.Opacity_calculator(x_choice = x_choice, 
                    x_coord = psi_dic[aa], ne = Ne_dic[aa], te = Te_dic[aa], 
                    neuden = Nd_dic[aa])
                
                
                tanh_ne_fit_dic[aa] = result_dic['tanh_ne_fit']
                tanh_te_fit_dic[aa] = result_dic['tanh_te_fit']
                exp_an_fit_dic[aa] = result_dic['exp_fit']
                delta_dic[aa] = result_dic['pedestal_width']
                tdelta_dic[aa] = result_dic['temperature_pedestal_width']
                xcoord_cut_dic[aa] = result_dic['x_coord_cut']
                efold_dic[aa] = result_dic['efold_length']
                opq_dic[aa] = result_dic['dimensionless_opaqueness']
                sym_pt_dic[aa] = result_dic['ne_symmetry_point']
                te_sym_pt_dic[aa] = result_dic['te_symmetry_point']
                
                exp_fit_m2_dic[aa] = result_dic['exp_fit_m2']
                x_m2_dic[aa] = result_dic['x_m2']
                efold_m2_dic[aa] = result_dic['efold_length_method2']
                opq_m2_dic[aa] = result_dic['dimensionless_opaqueness_method2']
                
                exp_fit_m3_dic[aa] = result_dic['exp_fit_m3']
                x_m3_dic[aa] = result_dic['x_m3']
                efold_m3_dic[aa] = result_dic['efold_length_method3']
                opq_m3_dic[aa] = result_dic['dimensionless_opaqueness_method3']
                
                # dsa_cut_dic[aa] = result_dic['dsa_cut']
                # dsa_psn_dic[aa] = result_dic['dsa_psn']
                # psi_psn_dic[aa] = result_dic['psi_psn']
            
            
            # result = {'efold_length': efold_dic, 'pedestal_width': delta_dic,
            #           'dimensionless_opaqueness': opq_dic, 
            #           'temperature_pedestal_width': tdelta_dic,
            #           'tanh_ne_fit': tanh_ne_fit_dic,
            #           'tanh_te_fit': tanh_te_fit_dic,
            #           'exp_fit_in_width': exp_an_fit_dic,
            #           'dsa_cut': dsa_cut_dic,
            #           'neutral_density': Nd_dic,
            #           'electron_density': Ne_dic,
            #           'electron_temperature_density': Te_dic,
            #           'psiN': psi_dic, 'dsa': dsa_pol_loc_dic, 
            #           'psi_psn': psi_psn_dic, 'dsa_psn': dsa_psn_dic
            #           }
            
            result = {'efold_length': efold_dic, 'pedestal_width': delta_dic,
                      'dimensionless_opaqueness': opq_dic, 
                      'temperature_pedestal_width': tdelta_dic,
                      'tanh_ne_fit': tanh_ne_fit_dic,
                      'tanh_te_fit': tanh_te_fit_dic,
                      'exp_fit': exp_an_fit_dic,
                      'xcoord_cut': xcoord_cut_dic,
                      'neutral_density': Nd_dic,
                      'electron_density': Ne_dic,
                      'electron_temperature_density': Te_dic,
                      'psiN': psi_dic, 'dsa': dsa_pol_loc_dic,
                      'te_symmetry_point': te_sym_pt_dic,
                      'ne_symmetry_point': sym_pt_dic,
                      
                      'exp_fit_m2': exp_fit_m2_dic,
                      'x_m2': x_m2_dic,
                      'efold_length_method2': efold_m2_dic,
                      'dimensionless_opaqueness_method2': opq_m2_dic,
                      
                      'exp_fit_m3': exp_fit_m3_dic,
                      'x_m3': x_m3_dic,
                      'efold_length_method3': efold_m3_dic,
                      'dimensionless_opaqueness_method3': opq_m3_dic
                      }
            
            SEP = self.data['DefaultSettings']['SEP']
            
            self.data['opacity_study'] = result
            
            char = {}
            char['withshift'] = self.withshift
            char['withseries'] = self.withseries
            
            density_dic = {}
            for k in self.data['dircomp']['Attempt'].keys():
                kk = float(k)*pow(10, 19)
                density_dic[k] = kk
            
            
            P = self.data['Parameter']
            opm.opacity_radial_method_multi(x_choice = x_choice, result_dic = result, 
                SEP = SEP, iter_list = self.data['dircomp']['Attempt'].keys(), 
                change_var_dic = density_dic, log_flag = True, char= char, P = P)
            
        
            
            
        elif self.withshift == True and self.withseries == True:
            print('Opacity_study_radial_plot_psi is not there yet, to be continue...')    
            
        else:
            print('Opacity_study_radial_plot_psi has a bug')
    
            

        
    def Opacity_study_radial_plot_backup(self, pol_loc):
           
        if self.withshift == False and self.withseries == False:
            self.load_output_data(param= 'NeuDen')
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
            
            psi = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 2]
            dsa_pol_loc = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]
            # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
            SEP = int(self.data['DefaultSettings']['SEP'])
            
            pol_index = int(pol_loc) + 1
            Nd = self.data['outputdata']['NeuDen'][:, pol_index]
            Ne = self.data['outputdata']['Ne'][:, pol_index]
            Te = self.data['outputdata']['Te'][:, pol_index]
            
            self.loadmastdata()
            fitpsi = self.data['fitprofile']['psi_normal']
            fitNe = self.data['fitprofile']['electron_density(m^(-3))']
            
            psi_to_dsa_func = interpolate.interp1d(psi, dsa_pol_loc, fill_value = 'extrapolate')
                  
            cutpsi = []
            cutNe = []
            n = len(fitpsi)
            for i in range(n):
                if fitpsi[i] >= min(psi):
                    cutpsi.append(fitpsi[i])
                    cutNe.append(fitNe[i])
            
            cutpsi = np.asarray(cutpsi)
            cutNe = np.asarray(cutNe)
            
            fitdsa = psi_to_dsa_func(cutpsi)
            
            result = fm.Opacity_calculator(x_choice = 'psiN', x_coord = psi, 
                                           ne = Ne, te = Te, neuden = Nd)
            
            efold = result['efold_length']
            sym_pt = result['ne_symmetry_point']
            te_sym_pt = result['te_symmetry_point']
            exp_an_fit = result['exp_fit']
            xcoord_cut = result['x_coord_cut']
            dn = result['pedestal_width']
            tanh_ne_fit = result['tanh_ne_fit']
            tanh_te_fit = result['tanh_te_fit']
            dtn = result['temperature_pedestal_width']
            
            x = [-efold + max(xcoord_cut), max(xcoord_cut)]
            y = [max(exp_an_fit), max(exp_an_fit)]
            xd = [-dn + sym_pt, dn + sym_pt]
            yd = [tanh_ne_fit[SEP], tanh_ne_fit[SEP]]
            xt = [-dtn + te_sym_pt, dtn + te_sym_pt]
            yt = [tanh_te_fit[SEP], tanh_te_fit[SEP]]
            
            
            log_flag = True
            
            plt.figure(figsize=(7,7))
            if log_flag:
                plt.yscale('log')
            else:
                pass
            plt.plot(psi, Nd,'o-', color = 'green', label= 'solps neutral density')
            # plt.plot(psi_RGI, Nd,'o-', color = 'b', label= 'RGI_solps neutral density')
            plt.plot(xcoord_cut, exp_an_fit, color='r',lw= 5, label= 'exponential fit')
            plt.axvline(x= max(xcoord_cut), color='orange',lw=3)
            plt.plot(x,y, color='orange', lw=3, label= 'Neutral penetration length [m]: $\lambda_{n_D}$')
            plt.axvline(x=-efold + max(xcoord_cut), color='orange',lw=3)
            plt.axvline(x=dn + sym_pt, color='black',lw=3, ls='--')
            plt.axvline(x=-dn + sym_pt, color='black',lw=3, ls='--')
            # plt.text(0.1, 5*pow(10, 16), 'dimensionless opaqueness: {}'.format(opq))
            plt.xlabel('psiN')
            plt.ylabel(self.data['Parameter']['NeuDen'])
            plt.title('Neutral density with fits')
            # plt.title(plot_dic['an3da.last10'][0],fontdict={"family":"Calibri","size": 20})
            plt.legend()
                
            # plt.subplot(211, sharex= ax1)
            plt.figure(figsize=(7,7))
            # plt.plot(psi_xport, Ne,'o-', color = 'r', label= 'solps_electron density')
            plt.plot(psi, Ne,'o-', color = 'b', label= 'solps electron density')
            # plt.plot(fitdsa, cutNe,'o-', color = 'g', label= 'experiment electron density')
            plt.plot(psi, tanh_ne_fit, color='r',lw= 3, label= 'tanh fit')
            plt.plot(xd, yd, color='black', lw=3, label= 'Pedestal width [m]: $\Delta n_e$')
            plt.axvline(x=dn + sym_pt, color='black',lw=3)
            plt.axvline(x=-dn + sym_pt, color='black',lw=3)
            plt.axvline(x=max(xcoord_cut), color='orange',lw=3, ls='--')
            plt.axvline(x=-efold + max(xcoord_cut), color='orange',lw=3, ls='--')
            plt.xlabel('psiN')
            plt.ylabel(self.data['Parameter']['Ne'])
            plt.title('Electron density with fits')
            # plt.title(plot_dic['ne3da.last10'][0],fontdict={"family":"Calibri","size": 20})
            plt.legend()
            
            plt.figure(figsize=(7,7))
            # plt.plot(psi_xport, Ne,'o-', color = 'r', label= 'solps_electron density')
            plt.plot(psi, Te,'o-', color = 'b', label= 'solps electron tempurature')
            # plt.plot(fitdsa, cutNe,'o-', color = 'g', label= 'experiment electron density')
            plt.plot(psi, tanh_te_fit, color='r',lw= 3, label= 'tanh fit')
            plt.plot(xt, yt, color='black', lw=3, label= 'temperature pedestal width [m]: $\Delta n_e$')
            plt.axvline(x=dtn + te_sym_pt, color='black',lw=3)
            plt.axvline(x=-dtn + te_sym_pt, color='black',lw=3)
            plt.xlabel('psiN')
            plt.ylabel(self.data['Parameter']['Te'])
            plt.title('Electron temperature with fits')
            plt.legend()
            
            
            plt.show()
        
        
        elif self.withshift == True and self.withseries == False:
            self.load_output_data(param= 'NeuDen')
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
            
            efold_dic = {}
            delta_dic = {}
            opq_dic = {}
            
            dsa_pol_loc_dic = {}
            dsa_cut_dic = {}
            Nd_dic = {}
            Ne_dic = {}
            Te_dic = {}
            exp_an_fit_dic = {}
            tanh_ne_fit_dic = {}
            
            
            for aa in self.data['dircomp']['multi_shift']:
                psi = self.data['psi']['psi_{}_val'.format(pol_loc)][aa]
                dsa_pol_loc_dic[aa] = self.data['psi']['dsa_{}_val'.format(pol_loc)][aa]
                # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
                
                pol_index = int(pol_loc) + 1
                Nd_dic[aa] = self.data['outputdata']['NeuDen'][aa][:, pol_index]
                Ne_dic[aa] = self.data['outputdata']['Ne'][aa][:, pol_index]
                Te_dic[aa] = self.data['outputdata']['Te'][aa][:, pol_index]
                
                
                fit_tanh_dic = fm.tanh_dsa_fit(dsa= dsa_pol_loc_dic[aa], 
                                               ne= Ne_dic[aa], te= Te_dic[aa])
                tanh_ne_fit_dic[aa] = fit_tanh_dic['tanh_ne_fit']
                delta_dic[aa] = fit_tanh_dic['popt_ne'][1]
                
                
                dsa_cut = []
                an_cut = []
                for j in range(len(dsa_pol_loc_dic[aa])):
                    if dsa_pol_loc_dic[aa][j] <= delta_dic[aa] and dsa_pol_loc_dic[aa][j] >= -delta_dic[aa]:
                        dsa_cut.append(dsa_pol_loc_dic[aa][j])
                        an_cut.append(Nd_dic[aa][j])
                
                dsa_cut_dic[aa] = np.asarray(dsa_cut)
                an_cut = np.asarray(an_cut)
                
                fit_exp_dic = fm.exp_dsa_fit(dsa= dsa_cut_dic[aa], neuden= an_cut)
                exp_an_fit_dic[aa] = fit_exp_dic['exp_an_fit']
                efold_dic[aa] = 1/fit_exp_dic['popt_an'][1]
                
                
                opq_dic[aa] = 2*delta_dic[aa]/efold_dic[aa]
            
            
            result_dic = {'tanh_fit': tanh_ne_fit_dic, 'exp_fit': exp_an_fit_dic,
                          'pedestal_width': delta_dic, 'efold_length': efold_dic,
                          'dimensionless_opaqueness': opq_dic}
            
            self.data['opacity_study'] = result_dic
            
            
            log_flag = True
            ii = 1
            
            for i in self.data['dircomp']['multi_shift']:
                plt.figure(ii)
                x = [-efold_dic[i] + max(dsa_cut_dic[i]), max(dsa_cut_dic[i])]
                y = [max(exp_an_fit_dic[i]), max(exp_an_fit_dic[i])]
                plt.plot(dsa_pol_loc_dic[i], Nd_dic[i],'o-', color = 'green', label= 'solps neutral density')
                plt.plot(dsa_cut_dic[i], exp_an_fit_dic[i], color='r',lw= 5, label= 'exponential fit')
                plt.axvline(x=max(dsa_cut_dic[i]), color='orange',lw=3)
                plt.plot(x,y, color='orange', lw=3, label= 'Neutral penetration length [m]: $\lambda_{n_D}$')
                plt.axvline(x = -efold_dic[i] + max(dsa_cut_dic[i]), color='orange',lw=3)
                plt.axvline(x=delta_dic[i], color='black',lw=3, ls='--')
                plt.axvline(x=-delta_dic[i], color='black',lw=3, ls='--')
                plt.xlabel('Radial coordinate: $R- R_{sep}$')
                plt.ylabel(self.data['Parameter']['NeuDen'])
                shift_value = self.data['dircomp']['shift_dic'][i]
                plt.title('Modify {} m Neutral density with fits'.format(shift_value))
                # plt.title(plot_dic['an3da.last10'][0],fontdict={"family":"Calibri","size": 20})
                plt.legend()
                ii = ii + 1
                if log_flag:
                    plt.yscale('log')
                else:
                    pass
            
            # plt.plot(psi_xport, Ne,'o-', color = 'r', label= 'solps_electron density')
            for j in self.data['dircomp']['multi_shift']:
                plt.figure(ii)
                xd = [-delta_dic[j], delta_dic[j]]
                yd = [tanh_ne_fit_dic[j][20], tanh_ne_fit_dic[j][20]]
                plt.plot(dsa_pol_loc_dic[j], Ne_dic[j],'o-', label= 'solps electron density_{}'.format(j))
                plt.plot(dsa_pol_loc_dic[j], tanh_ne_fit_dic[j], color='r',lw= 3, label= 'exponential fit')
                plt.plot(xd, yd, color='black', lw=3, label= 'Pedestal width [m]: $\Delta n_e$')
                plt.axvline(x= delta_dic[j], color='black',lw=3)
                plt.axvline(x=-delta_dic[j], color='black',lw=3)
                plt.axvline(x=max(dsa_cut_dic[i]), color='orange',lw=3, ls='--')
                plt.axvline(x=-efold_dic[i] + max(dsa_cut_dic[i]), color='orange',lw=3, ls='--')
                plt.xlabel('Radial coordinate: $R- R_{sep}$')
                plt.ylabel(self.data['Parameter']['Ne'])
                plt.title('Electron density with fits')
                plt.legend()
                ii = ii + 1
            
            
            shift_ar = {}
            for k in self.data['dircomp']['multi_shift']:
                kk = float(self.data['dircomp']['shift_dic'][k])
                shift_ar[k] = kk
            
            
            plt.figure(ii)
            for j in self.data['dircomp']['multi_shift']:
                plt.plot(shift_ar[j], opq_dic[j],'o-', color= 'r', label= 'dimensionless opaqueness')
            plt.xlabel('shift: [m]')
            plt.title('opaqueness verses shift distance')
            # plt.title(plot_dic['ne3da.last10'][0],fontdict={"family":"Calibri","size": 20})
            
            ii = ii + 1
            plt.figure(ii)
            for j in self.data['dircomp']['multi_shift']:
                plt.plot(shift_ar[j], efold_dic[j],'o-', color= 'r', label= 'neutral penetration')
            plt.xlabel('shift: [m]')
            plt.ylabel('neutral penetration length: [m]')
            plt.title('neutral penetration length verses shift distance')
            
            ii = ii + 1
            plt.figure(ii)
            for p in self.data['dircomp']['multi_shift']:
                plt.plot(shift_ar[p], delta_dic[p],'o-', color= 'r', label= 'neutral penetration')
            plt.xlabel('shift: [m]')
            plt.ylabel('pedestal width: [m]')
            plt.title('pedestal width verses shift distance')
        
        
        elif self.withshift == False and self.withseries == True:
            self.load_output_data(param= 'NeuDen')
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
            
            efold_dic = {}
            delta_dic = {}
            opq_dic = {}
            
            dsa_pol_loc_dic = {}
            dsa_cut_dic = {}
            Nd_dic = {}
            Ne_dic = {}
            Te_dic = {}
            exp_an_fit_dic = {}
            tanh_ne_fit_dic = {}
            
            
            for aa in self.data['dircomp']['Attempt'].keys():
                psi = self.data['psi']['psi_{}_val'.format(pol_loc)]
                dsa_pol_loc_dic[aa] = self.data['psi']['dsa_{}_val'.format(pol_loc)]
                # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
                
                pol_index = int(pol_loc) + 1
                Nd_dic[aa] = self.data['outputdata']['NeuDen'][aa][:, pol_index]
                Ne_dic[aa] = self.data['outputdata']['Ne'][aa][:, pol_index]
                Te_dic[aa] = self.data['outputdata']['Te'][aa][:, pol_index]
                
                
                fit_tanh_dic = fm.tanh_dsa_fit(dsa= dsa_pol_loc_dic[aa], 
                                               ne= Ne_dic[aa], te= Te_dic[aa])
                tanh_ne_fit_dic[aa] = fit_tanh_dic['tanh_ne_fit']
                delta_dic[aa] = fit_tanh_dic['popt_ne'][1]
                
                
                dsa_cut = []
                an_cut = []
                for j in range(len(dsa_pol_loc_dic[aa])):
                    if dsa_pol_loc_dic[aa][j] <= delta_dic[aa] and dsa_pol_loc_dic[aa][j] >= -delta_dic[aa]:
                        dsa_cut.append(dsa_pol_loc_dic[aa][j])
                        an_cut.append(Nd_dic[aa][j])
                
                dsa_cut_dic[aa] = np.asarray(dsa_cut)
                an_cut = np.asarray(an_cut)
                
                fit_exp_dic = fm.exp_dsa_fit(dsa= dsa_cut_dic[aa], neuden= an_cut)
                exp_an_fit_dic[aa] = fit_exp_dic['exp_an_fit']
                efold_dic[aa] = 1/fit_exp_dic['popt_an'][1]
                
                
                opq_dic[aa] = 2*delta_dic[aa]/efold_dic[aa]
            
            
            result_dic = {'tanh_fit': tanh_ne_fit_dic, 'exp_fit': exp_an_fit_dic,
                          'pedestal_width': delta_dic, 'efold_length': efold_dic,
                          'dimensionless_opaqueness': opq_dic}
            
            self.data['opacity_study'] = result_dic
            
            density_dic = {}
            for k in self.data['dircomp']['Attempt'].keys():
                kk = float(k)*pow(10, 19)
                density_dic[k] = kk
            
            log_flag = True
            ii = 1
            
            for i in self.data['dircomp']['Attempt'].keys():
                plt.figure(ii)
                x = [-efold_dic[i] + max(dsa_cut_dic[i]), max(dsa_cut_dic[i])]
                y = [max(exp_an_fit_dic[i]), max(exp_an_fit_dic[i])]
                plt.plot(dsa_pol_loc_dic[i], Nd_dic[i],'o-', color = 'green', label= 'solps neutral density')
                plt.plot(dsa_cut_dic[i], exp_an_fit_dic[i], color='r',lw= 5, label= 'exponential fit')
                plt.axvline(x=max(dsa_cut_dic[i]), color='orange',lw=3)
                plt.plot(x,y, color='orange', lw=3, label= 'Neutral penetration length [m]: $\lambda_{n_D}$')
                plt.axvline(x = -efold_dic[i] + max(dsa_cut_dic[i]), color='orange',lw=3)
                plt.axvline(x=delta_dic[i], color='black',lw=3, ls='--')
                plt.axvline(x=-delta_dic[i], color='black',lw=3, ls='--')
                plt.xlabel('Radial coordinate: $R- R_{sep}$')
                plt.ylabel(self.data['Parameter']['NeuDen'])
                plt.title('Neutral density with fits and electron density is {}'.format(density_dic[i]))
                # plt.title(plot_dic['an3da.last10'][0],fontdict={"family":"Calibri","size": 20})
                plt.legend()
                ii = ii + 1
                if log_flag:
                    plt.yscale('log')
                else:
                    pass
            
            # plt.plot(psi_xport, Ne,'o-', color = 'r', label= 'solps_electron density')
            for j in self.data['dircomp']['Attempt'].keys():
                plt.figure(ii)
                xd = [-delta_dic[j], delta_dic[j]]
                yd = [tanh_ne_fit_dic[j][20], tanh_ne_fit_dic[j][20]]
                plt.plot(dsa_pol_loc_dic[j], Ne_dic[j],'o-', label= 'solps electron density_{}'.format(j))
                plt.plot(dsa_pol_loc_dic[j], tanh_ne_fit_dic[j], color='r',lw= 3, label= 'exponential fit')
                plt.plot(xd, yd, color='black', lw=3, label= 'Pedestal width [m]: $\Delta n_e$')
                plt.axvline(x= delta_dic[j], color='black',lw=3)
                plt.axvline(x=-delta_dic[j], color='black',lw=3)
                plt.axvline(x=max(dsa_cut_dic[i]), color='orange',lw=3, ls='--')
                plt.axvline(x=-efold_dic[i] + max(dsa_cut_dic[i]), color='orange',lw=3, ls='--')
                plt.xlabel('Radial coordinate: $R- R_{sep}$')
                plt.ylabel(self.data['Parameter']['Ne'])
                plt.title('Electron density with fits and core boundary is {}'.format(density_dic[i]))
                plt.legend()
                ii = ii + 1
            
            
            plt.figure(ii)
            for s in self.data['dircomp']['Attempt'].keys():
                plt.plot(density_dic[s], opq_dic[s],'o-', color= 'r', label= 'dimensionless opaqueness')
            plt.xlabel('core density: [$m^{-3}$]')
            plt.title('dimensionless opaqueness verses different core density')
            # plt.title(plot_dic['ne3da.last10'][0],fontdict={"family":"Calibri","size": 20})
            
            ii = ii + 1
            plt.figure(ii)
            for k in self.data['dircomp']['Attempt'].keys():
                plt.plot(density_dic[k], efold_dic[k],'o-', color= 'r', label= 'neutral penetration')
            plt.xlabel('core density: [$m^{-3}$]')
            plt.ylabel('neutral penetration length: [m]')
            plt.title('neutral penetration length verses different core density')
            
            ii = ii + 1
            plt.figure(ii)
            for p in self.data['dircomp']['Attempt'].keys():
                plt.plot(density_dic[p], delta_dic[p],'o-', color= 'r', label= 'neutral penetration')
            plt.xlabel('core density: [$m^{-3}$]')
            plt.ylabel('pedestal width: [m]')
            plt.title('pedestal width verses different core density')
        
        
        elif self.withshift == True and self.withseries == True:
            print('load_output_data is not there yet, to be continue...')
        
        else:
            print('plot_Ne_NeuDen_single function has a bug')
            
            
            
        
            
            
        
