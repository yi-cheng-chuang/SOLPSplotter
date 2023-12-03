# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 11:47:48 2023

@author: user
"""
from B2plotter_PRmap import RP_mapping
import B2plotter_set as b2s
import opacity_plot_method as opm
import matplotlib.pyplot as plt
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate
from scipy.optimize import curve_fit
import numpy as np
import B2TransportParser as b2tp


class Opacity_study(RP_mapping):
    def __init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters, Publish):
        RP_mapping.__init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters)
        
        self.Publish = Publish
        self.data['DefaultSettings']['Publish'] = self.Publish
    
    
    def set_plot(self):
        if self.Publish == 'b2plottersetting':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth=5, markersize=9)
            plt.rcParams.update({'font.size': 16})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
  
        else:
            print('Publish setting is incorrect or add another setting')
    

    def Opacity_study_poloidal_plot(self, pol_list):
        self.data['poloidal_index'] = pol_list
        
        for j in pol_list:
            self.calcpsi_1D(pol_loc= j)
            self.calc_dsa(pol_loc= j)
        
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
            result = self.opacity_data_method_single(pol_list = pol_list)
            
            self.data['opacity_poloidal'] = result
            
            unit = b2s.opacity_study_unit()
            
            char = {}
            char['withshift'] = self.withshift
            char['withseries'] = self.withseries
            
            log_flag = False
            
            # print(result.keys())
            # print(unit.keys())
            
            opm.opacity_plot(pol_loc = self.data['angle']['angle_list'], result_dic = result, unit_dic = unit,
                             log_flag = False, charactor= char,
                             iter_list = None, change_ver_dic = None,
                             xpoint_loc = self.data['angle']['xpoint_angle'])
            
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
                lb = np.asarray(result['dimensionless_opaqueness'][la])
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
            char['series_flag'] = self.DefaultSettings['series_flag']

            opm.data_reorder(iter_list = pol_list, change_var = shift_list,
                             data_collect = data_collect_opq, char = char)
            
            unit = b2s.opacity_study_unit()
                       
            shift_dic = {}
            for k in self.data['dircomp']['multi_shift']:
                p = str(self.data['dircomp']['shift_dic'][k])
                shift_dic[k] = p
            
            opm.opacity_plot(pol_loc = self.data['angle']['angle_list'], result_dic = result, unit_dic = unit,
                             log_flag = False, charactor= char,
                             iter_list = self.data['dircomp']['multi_shift'], 
                             change_ver_dic = shift_dic, 
                             xpoint_loc= self.data['angle']['xpoint_angle'])
        
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
            variable = 'dimensionless_opaqueness'
            for la in self.data['dircomp']['Attempt'].keys():
                lb = np.asarray(result[variable][la])
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
            char['series_flag'] = self.DefaultSettings['series_flag']
            char['variable'] = variable

            opm.data_reorder(iter_list = pol_list, change_var = density_list,
                             data_collect = data_collect_opq, char = char)
            
            unit = b2s.opacity_study_unit()
            
            
            density_dic = {}
            for k in self.data['dircomp']['Attempt'].keys():
                kk = float(k)*pow(10, 19)
                density_dic[k] = kk
            
            opm.opacity_plot(pol_loc = self.data['angle'], result_dic = result, unit_dic = unit,
                             log_flag = False, charactor= char,
                             iter_list = self.data['dircomp']['Attempt'].keys(), 
                             change_ver_dic = density_dic, xpoint_loc= None)
        
        elif self.withshift == True and self.withseries == True:
            print('Opacity_study_poloidal_plot is not there yet, to be continue...')
            
            
        else:
            print('more work need to be done')
            
            
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
        
            
            
            for k in pol_list:
                
                if self.withshift == True and self.withseries == False:
                    psi = self.data['psi']['psi_{}_val'.format(k)][aa]
                elif self.withshift == False and self.withseries == True:
                    psi = self.data['psi']['psi_{}_val'.format(k)]
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
        
        return result
          
            
    def Opacity_study_radial_plot(self, pol_loc):
        self.load_output_data(param= 'NeuDen')
        self.load_output_data(param= 'Ne')
        self.load_output_data(param= 'Te')
        
        if self.withshift == False and self.withseries == False:
          
            psi = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 2]
            # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
            SEP = int(self.data['DefaultSettings']['SEP'])
            
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
            SEP = self.data['DefaultSettings']['SEP']
            
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
    
    
    
    def transport_coe_align_plot(self, plot_transcoe):
        if self.withshift == True and self.withseries == False:
            trans_dic = {}
            jxa = self.data['b2mn']['org']['jxa']
            self.calcpsi_1D(pol_loc= jxa)
            for aa in self.data['dircomp']['multi_shift']:
                trans_file_dir = self.data['dirdata']['infolderdir'][aa]['simudir'] + '/b2.transport.inputfile'
                
                trans_list = b2tp.InputfileParser(trans_file_dir, plot= False)
                cod = trans_list['1'].T
                coki = trans_list['3'].T
                coke = trans_list['4'].T
                x= cod[:,0]  #the coordinate here is R-R_sep
                
                trans_dic[aa] = np.zeros([len(x), 4])
                trans_dic[aa][:, 0] = self.data['psi']['psi_{}_val'.format(jxa)][aa]
                trans_dic[aa][:, 1] = cod[:, 1]
                trans_dic[aa][:, 2] = coki[:, 1]
                trans_dic[aa][:, 3] = coke[:, 1]
    
   
            log_flag = False
            coe_label_dic = {'1': 'particle diffusivity', '2': 'ion thermal diffusivity'
                             ,'3': 'electron thermal diffusivity'}
            if plot_transcoe:
                for k in coe_label_dic.keys():
                    if log_flag:
                        plt.yscale('log')
                    plt.figure(figsize=(7,7))
                    color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                                 'dot7': 'blue', 'one': 'purple'}
                    for ab in self.data['dircomp']['multi_shift']: 
                        # plt.plot(trans_dic[ab][:, 0], trans_dic[ab][:, int(k)], 'o-', color= color_dic[ab],
                        #          label ='transport coefficient of modify {} m case'.format(self.data['dircomp']['shift_dic'][ab]))
                        plt.plot(trans_dic[ab][:, 0], trans_dic[ab][:, int(k)], 'o-', color= color_dic[ab] )
                        plt.xlabel('psiN')
                        plt.title('radial {} coefficient'.format(coe_label_dic[k]))
                        # plt.legend() 
                plt.show()
        else:
            print('transport_coe_align_plot is not there yet')
    
            

        
 