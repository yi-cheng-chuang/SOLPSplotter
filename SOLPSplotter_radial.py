# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 21:00:45 2024

@author: ychuang
"""

from SOLPSplotter_fit import profile_fit
import opacity_plot_method as opm
import matplotlib.pyplot as plt
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate
from scipy.optimize import curve_fit
import numpy as np


class radial_plot(profile_fit):
    
    def __init__(self, DefaultSettings, loadDS):
        profile_fit.__init__(self, DefaultSettings, loadDS)
            
    
    def set_plot(self):
        if self.Publish == 'b2plottersetting':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 5, markersize= 9)
            plt.rcParams.update({'font.size': 16})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
  
        else:
            print('Publish setting is incorrect or add another setting')
    
    
    def opacity_radial_method(self, result_dic, SEP, x_coord, Nd, Ne, Te, 
                                          P, log_flag):
        
        tanh_ne_fit = result_dic['tanh_ne_fit']
        tanh_te_fit = result_dic['tanh_te_fit']
        exp_an_fit = result_dic['exp_fit']
        dn = result_dic['pedestal_width'][0]
        dtn = result_dic['temperature_pedestal_width']
        efold = result_dic['efold_length'][0]
        opq = result_dic['dimensionless_opaqueness']
        xcoord_cut = result_dic['x_coord_cut']
        sym_pt = result_dic['ne_symmetry_point']
        te_sym_pt = result_dic['te_symmetry_point']
        
                      
        x = [-efold + max(xcoord_cut), max(xcoord_cut)]
        y = [min(exp_an_fit), min(exp_an_fit)]
        xd = [-dn + sym_pt, dn + sym_pt]
        yd = [tanh_ne_fit[SEP] , tanh_ne_fit[SEP]]
        xt = [-dtn + te_sym_pt, dtn + te_sym_pt]
        yt = [tanh_te_fit[SEP], tanh_te_fit[SEP]]
        
        
        plt.figure(figsize=(7,7))
        if log_flag:
            plt.yscale('log')
        else:
            pass
        plt.plot(x_coord, Nd,'-', color = 'green', label= 'solps neutral density')
        # plt.plot(psi_RGI, Nd,'-', color = 'b', label= 'RGI_solps neutral density')
        plt.plot(xcoord_cut, exp_an_fit, color='r',lw= 5, ls='-', label= 'exponential fit')
        plt.axvline(x= max(xcoord_cut), color='orange',lw=3)
        plt.plot(x,y, color='orange', lw=3, ls='-', label= 'Neutral penetration length : $\lambda_{n_D}$')
        plt.axvline(x=-efold + max(xcoord_cut), color='orange',lw=3)
        plt.axvline(x= max(xcoord_cut), color='black',lw=3, ls='--', 
                    label= 'fit range : $\Delta n_e$')
        plt.axvline(x= min(xcoord_cut), color='black',lw=3, ls='--')
        # plt.axvline(x= x_m2[0], color='purple',lw=3, ls='--', 
        #             label= 'exp fitting width')
        # plt.axvline(x= x_m2[-1], color='purple',lw=3, ls='--')
        plt.xlabel('psiN')
        # plt.ylabel(P['NeuDen'])
        plt.title('Neutral density with fits')
        plt.legend()
            
        # plt.subplot(211, sharex= ax1)
        plt.figure(figsize=(7,7))
        if log_flag:
            plt.yscale('log')
        else:
            pass
        # plt.plot(psi_xport, Ne,'-', color = 'r', label= 'solps_electron density')
        plt.plot(x_coord, Ne,'-', color = 'b', label= 'solps electron density')
        # plt.plot(fitdsa, cutNe,'-', color = 'g', label= 'experiment electron density')
        plt.plot(x_coord, tanh_ne_fit, ls='-', color='r',lw= 3, label= 'tanh fit')
        plt.plot(xd, yd, color='black', ls='-', lw=3, label= 'Pedestal width : $\Delta n_e$')
        plt.axvline(x=dn + sym_pt, color='black',lw=3)
        plt.axvline(x=-dn + sym_pt, color='black',lw=3)
        plt.axvline(x=max(xcoord_cut), color='orange',lw=3, ls='--', 
                    label = 'Neutral penetration length : $\lambda_{n_D}$')
        plt.axvline(x=-efold + max(xcoord_cut), color='orange',lw=3, ls='--')
        plt.xlabel('psiN')
        # plt.ylabel(P['Ne'])
        plt.title('Electron density with fits')
        # plt.title(plot_dic['ne3da.last10'][0],fontdict={"family":"Calibri","size": 20})
        plt.legend()
        
        plt.figure(figsize=(7,7))
        if log_flag:
            plt.yscale('log')
        else:
            pass
        plt.plot(x_coord, Te, ls='-', color = 'b', label= 'solps electron tempurature')
        plt.plot(x_coord, tanh_te_fit, ls='-', color='r',lw= 3, label= 'tanh fit')
        plt.plot(xt, yt, color='black', ls='-', lw=3, label= 'temperature pedestal width [m]: $\Delta n_e$')
        plt.axvline(x=dtn + te_sym_pt, color='black', lw=3, ls='-')
        plt.axvline(x=-dtn + te_sym_pt, color='black', lw=3, ls='-')
        plt.xlabel('psiN')
        # plt.ylabel(P['Te'])
        plt.title('Electron temperature with fits')
        plt.legend()
    
    
    def Opacity_study_radial_plot(self, pol_loc):
        
        if self.withshift == False and self.withseries == False:
            
            pol_index = int(pol_loc[0])
            Nd = self.data['radial_fit_data']['NeuDen'][:, pol_index]
            Ne = self.data['radial_fit_data']['Ne'][:, pol_index]
            Te = self.data['radial_fit_data']['Te'][:, pol_index]
            SEP = int(self.data['DefaultSettings']['sep_index_dsa'])
            psi = self.data['psi']['psi_{}_val'.format(pol_loc[0])][:, 2]
            
            
            result_dic = self.data['radial_fit_data'] | self.data['opacity_poloidal']
            
            
            P = self.data['Parameter']
            self.opacity_radial_method(self, result_dic = result_dic, SEP = SEP, 
            x_coord = psi, Nd = Nd, Ne = Ne, Te = Te, P = P, log_flag = True)
        
        
        elif self.withshift == True and self.withseries == False:
            
            mix_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
            
                result_dic = self.data['radial_fit_data'][aa] | self.data['opacity_poloidal'][aa]
                mix_dic[aa] = result_dic
            
            self.data['mix_dic'] = mix_dic
            
            for aa in self.data['dircomp']['multi_shift']:
            
                pol_index = int(pol_loc[0])
                Nd = self.data['radial_fit_data'][aa]['NeuDen']
                Ne = self.data['radial_fit_data'][aa]['Ne']
                Te = self.data['radial_fit_data'][aa]['Te']
                SEP = int(self.data['DefaultSettings']['sep_index_dsa'][aa])
                psi = self.data['psi']['psi_{}_val'.format(pol_loc[0])][aa][:, 2]
                
                
                result_dic = self.data['radial_fit_data'][aa] | self.data['opacity_poloidal'][aa]
                mix_dic[aa] = result_dic
                
                P = self.data['Parameter']
                self.opacity_radial_method(result_dic = result_dic, SEP = SEP, 
                x_coord = psi, Nd = Nd, Ne = Ne, Te = Te, P = P, log_flag = True)
           
            
        
        elif self.withshift == False and self.withseries == True:
            
            
            for aa in self.data['dircomp']['Attempt'].keys():
            
                pol_index = int(pol_loc[0])
                Nd = self.data['radial_fit_data'][aa]['NeuDen'][:, pol_index]
                Ne = self.data['radial_fit_data'][aa]['Ne'][:, pol_index]
                Te = self.data['radial_fit_data'][aa]['Te'][:, pol_index]
                SEP = int(self.data['DefaultSettings'][aa]['sep_index_dsa'])
                psi = self.data['psi']['psi_{}_val'.format(pol_loc[0])][aa][:, 2]
                
                
                result_dic = self.data['radial_fit_data'][aa] | self.data['opacity_poloidal'][aa]
                
                
                P = self.data['Parameter']
                self.opacity_radial_method(result_dic = result_dic, SEP = SEP, 
                x_coord = psi, Nd = Nd, Ne = Ne, Te = Te, P = P, log_flag = True)
        
        elif self.withshift == True and self.withseries == True:
            print('Opacity_study_radial_plot_psi is not there yet, to be continue...')    
            
        else:
            print('Opacity_study_radial_plot_psi has a bug')
        
        
        
        
        
"""      
        
    
    def Opacity_study_radial_plot(self, pol_loc):
        self.load_output_data(param= 'NeuDen')
        self.load_output_data(param= 'Ne')
        self.load_output_data(param= 'Te')
        
        
        if self.withshift == False and self.withseries == False:
          
            psi = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 2]
            # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
            SEP = int(self.data['DefaultSettings']['sep_index_dsa'])
            
            pol_index = int(pol_loc)
            Nd = self.data['outputdata']['NeuDen'][:, pol_index]
            Ne = self.data['outputdata']['Ne'][:, pol_index]
            Te = self.data['outputdata']['Te'][:, pol_index]
            
            
            # result_dic = fm.Opacity_calculator(dsa_pol_loc, Ne, Te, Nd, psi)
            result_dic = fm.Opacity_calculator(x_coord = psi, 
                                           ne = Ne, te = Te, neuden = Nd)
            
            ped_index = result_dic['sep_index']
            flux_expand = self.calc_flux_expansion(pol_loc= pol_loc, 
                                    ped_index= ped_index, iter_index= None)
            
            
            flux_expand_dic = {'flux_expansion': flux_expand}
            opac_dic = result_dic | flux_expand_dic
            
            self.data['opacity_study'] = opac_dic
            
            
            P = self.data['Parameter']
            opm.opacity_radial_method_single(result_dic = result_dic, SEP = SEP, 
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
            sym_pt_dic = {}
            te_sym_pt_dic = {}
            
            
            psi_dic = {}
            dsa_pol_loc_dic = {}
            
            Nd_dic = {}
            Ne_dic = {}
            Te_dic = {}
            exp_an_fit_dic = {}
            tanh_ne_fit_dic = {}
            tanh_te_fit_dic = {}
            flux_expand_dic = {}
            
            
            
            for aa in self.data['dircomp']['multi_shift']:
                psi_dic[aa] = self.data['psi']['psi_{}_val'.format(pol_loc)][aa]
                # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
                
                
                pol_index = int(pol_loc)
                Nd_dic[aa] = self.data['outputdata']['NeuDen'][aa][:, pol_index]
                Ne_dic[aa] = self.data['outputdata']['Ne'][aa][:, pol_index]
                Te_dic[aa] = self.data['outputdata']['Te'][aa][:, pol_index]
                
                
                result_dic = fm.Opacity_calculator(x_coord = psi_dic[aa][:, 2], 
                            ne = Ne_dic[aa], te = Te_dic[aa], neuden = Nd_dic[aa])
                
                ped_index = result_dic['sep_index']
                
                flux_expand = self.calc_flux_expansion(pol_loc= pol_loc, 
                                    ped_index = ped_index, iter_index= aa)
                
                

                
                flux_expand_dic[aa] = flux_expand
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
                      'flux_expand': flux_expand_dic               
                      }
            
            
            
            rl = self.data['dircomp']['multi_shift']
            SEP = self.data['DefaultSettings']['sep_index_dsa']
            
            self.data['opacity_study'] = result
            
            char = {}
            char['withshift'] = self.withshift
            char['withseries'] = self.withseries
            char['series_flag'] = self.DefaultSettings['series_flag']
            
            shift_dic = {}
            for k in self.data['dircomp']['multi_shift']:
                p = str(self.data['dircomp']['shift_dic'][k])
                shift_dic[k] = p
            
            log_flag = True
            ii = 1
            
            P = self.data['Parameter']
            opm.opacity_radial_method_multi(result_dic = result, SEP = SEP, 
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
            
            Nd_dic = {}
            Ne_dic = {}
            Te_dic = {}
            exp_an_fit_dic = {}
            tanh_ne_fit_dic = {}
            tanh_te_fit_dic = {}
            flux_expand_dic = {}
            
            
            for aa in self.data['dircomp']['Attempt'].keys():
                psi_dic[aa] = self.data['psi']['psi_{}_val'.format(pol_loc)]
                # dsa_pol_loc_dic[aa] = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]
                # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
                
                
                pol_index = int(pol_loc)
                Nd_dic[aa] = self.data['outputdata']['NeuDen'][aa][:, pol_index]
                Ne_dic[aa] = self.data['outputdata']['Ne'][aa][:, pol_index]
                Te_dic[aa] = self.data['outputdata']['Te'][aa][:, pol_index]
                
                
                result_dic = fm.Opacity_calculator(x_coord = psi_dic[aa], 
                        ne = Ne_dic[aa], te = Te_dic[aa], neuden = Nd_dic[aa])
                
                ped_index = result_dic['sep_index']
                
                flux_expand = self.calc_flux_expansion(pol_loc= pol_loc, 
                                    ped_index = ped_index, iter_index= aa)
                
                flux_expand_dic[aa] = flux_expand           
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
                      'flux_expand': flux_expand_dic
                      
                      }
            
            SEP = self.data['DefaultSettings']['SEP']
            
            self.data['opacity_study'] = result
            
            char = {}
            char['withshift'] = self.withshift
            char['withseries'] = self.withseries
            char['series_flag'] = self.DefaultSettings['series_flag']
            
            density_dic = {}
            if self.DefaultSettings['series_flag'] == 'change_den':
                for k in self.data['dircomp']['Attempt'].keys():
                    kk = float(k)*pow(10, 19)
                    density_dic[k] = kk
            elif self.DefaultSettings['series_flag'] == 'eireneN':
                for k in self.data['dircomp']['Attempt'].keys():
                    density_dic[k] = k
                
            
            
            P = self.data['Parameter']
            opm.opacity_radial_method_multi(result_dic = result, 
                SEP = SEP, iter_list = self.data['dircomp']['Attempt'].keys(), 
                change_var_dic = density_dic, log_flag = True, char= char, P = P)
        
        
        elif self.withshift == True and self.withseries == True:
            print('Opacity_study_radial_plot_psi is not there yet, to be continue...')    
            
        else:
            print('Opacity_study_radial_plot_psi has a bug')
            
            
"""




    
    
    