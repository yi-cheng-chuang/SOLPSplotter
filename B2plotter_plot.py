# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 11:47:48 2023

@author: user
"""
from B2plotter_class import load_data
import matplotlib.pyplot as plt
import load_mast_expdata_method as lmem
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate
from scipy.optimize import curve_fit
import numpy as np


class simple_plot(load_data):
    def __init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters, Publish):
        load_data.__init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters)
        
        self.Publish = Publish
        self.data['DefaultSettings']['Publish'] = self.Publish
    
    def set_plot(self):
        if self.Publish == 'b2plottersetting':
            plt.rcParams.update({'font.weight': 'bold'})
            plt.rc('lines', linewidth=2, markersize=11)
            plt.rcParams.update({'font.size': 16})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
            # if Publish != []:
            #     plt.rc('font',size=30)
            #     plt.rc('lines',linewidth=5,markersize=15)
            # else:
            #     plt.rc('font',size=14)
        else:
            print('Publish setting is incorrect or add another setting')
    

    
    def Opacity_study_poloidal_plot(self, pol_list):
        
        
        self.data['poloidal_index'] = pol_list
        
        for j in pol_list:
            self.calcpsi_1D(pol_loc= j)
        
        
        if self.withshift == False and self.withseries == False:
            ln = len(pol_list)
            efold = np.zeros(ln)
            delta = np.zeros(ln)
            opq = np.zeros(ln)
            pol_loc = np.zeros(ln)
            neu_den = np.zeros(ln)
            
            self.load_output_data(param= 'NeuDen')
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
            i = 0
            for k in self.data['poloidal_index']:
                psi = self.data['psi']['psi_{}_val'.format(k)][:, 2]
                dsa_pol_loc = self.data['psi']['dsa_{}_val'.format(k)]
                # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
                # SEP = int(self.data['DefaultSettings']['SEP'])
                
                pol_in = int(k) + 1
                Nd = self.data['outputdata']['NeuDen'][:, pol_in]
                Ne = self.data['outputdata']['Ne'][:, pol_in]
                Te = self.data['outputdata']['Te'][:, pol_in]
                
                rd = fm.Opacity_calculator_psi(psi, Ne, Te, Nd)
                
                efold[i] = rd['efold_length']
                delta[i] = rd['pedestal_width']
                opq[i] = rd['dimensionless_opaqueness']
                pol_loc[i] = int(k)
                neu_den[i] = max(rd['exp_fit'])
                i = i + 1
            
            log_flag = False
            
            plt.figure(1)
            if log_flag:
                plt.yscale('log')
            else:
                pass
            plt.plot(pol_loc, efold,'o-', color = 'green', label= 'Neutral penetration length')
            plt.xlabel('poloidal index')
            plt.ylabel('$\lambda_{n_D}$: [m]')
            plt.title('Neutral penetration length in different poloidal index')
            plt.legend()
            
            
            plt.figure(2)
            if log_flag:
                plt.yscale('log')
            else:
                pass
            plt.plot(pol_loc, opq,'o-', color = 'green', label= 'dimensionless opaqueness')
            plt.xlabel('poloidal index')
            plt.title('dimensionless opaqueness in different poloidal index')
            plt.legend()
            
            plt.figure(3)
            if log_flag:
                plt.yscale('log')
            else:
                pass
            plt.plot(pol_loc, neu_den,'o-', color = 'green', label= 'neutral density')
            plt.xlabel('poloidal index')
            plt.title('neutral density in different poloidal index')
            plt.legend()
            
                  
        elif self.withshift == True and self.withseries == False:
            self.load_output_data(param= 'NeuDen')
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
            
            efold_dic = {}
            delta_dic = {}
            opq_dic = {}
            neu_den_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                i = 0
                ln = len(pol_list)
                efold = np.zeros(ln)
                delta = np.zeros(ln)
                opq = np.zeros(ln)
                pol_loc = np.zeros(ln)
                neu_den = np.zeros(ln)
                for k in self.data['poloidal_index']:
                    psi = self.data['psi']['psi_{}_val'.format(k)][aa]
                    dsa_pol_loc = self.data['psi']['dsa_{}_val'.format(k)]
                    # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
                    # SEP = int(self.data['DefaultSettings']['SEP'])
                    
                    pol_in = int(k) + 1
                    Nd = self.data['outputdata']['NeuDen'][aa][:, pol_in]
                    Ne = self.data['outputdata']['Ne'][aa][:, pol_in]
                    Te = self.data['outputdata']['Te'][aa][:, pol_in]
                    
                    rd = fm.Opacity_calculator_psi(psi, Ne, Te, Nd)
                    
                    efold[i] = rd['efold_length']
                    delta[i] = rd['pedestal_width']
                    opq[i] = rd['dimensionless_opaqueness']
                    neu_den[i] = max(rd['exp_fit'])
                    pol_loc[i] = int(k)
                    i = i + 1
                efold_dic[aa] = efold
                delta_dic[aa] = delta
                opq_dic[aa] = opq
                neu_den_dic[aa] = neu_den
            
            result = {'efold_length': efold_dic, 'pedestal_width': delta_dic,
                      'dimensionless_opaqueness': opq_dic}
            self.data['opacity_poloidal'] = result
            
            
            log_flag = False
            
            plt.figure(1)
            if log_flag:
                plt.yscale('log')
            else:
                pass
            for aa in self.data['dircomp']['multi_shift']:
                p = str(self.data['dircomp']['shift_dic'][aa])
                plt.plot(pol_loc, efold_dic[aa],'o-', 
                         label= 'Neutral penetration length for modified {} m case'.format(p))
            plt.xlabel('poloidal index')
            plt.ylabel('$\lambda_{n_D}$: [m]')
            plt.title('Neutral penetration length in different poloidal index')
            plt.legend()
            
            plt.figure(2)
            if log_flag:
                plt.yscale('log')
            else:
                pass
            for aa in self.data['dircomp']['multi_shift']:
                q = str(self.data['dircomp']['shift_dic'][aa])
                plt.plot(pol_loc, opq_dic[aa],'o-',
                         label= 'dimensionless opaqueness for modified {} m case'.format(q))
            plt.xlabel('poloidal index')
            plt.title('dimensionless opaqueness in different poloidal index')
            plt.legend()
            
            plt.figure(3)
            if log_flag:
                plt.yscale('log')
            else:
                pass
            for ab in self.data['dircomp']['multi_shift']:
                r = str(self.data['dircomp']['shift_dic'][ab])
                plt.plot(pol_loc, neu_den_dic[ab],'o-', 
                         label= 'neutral density for modified {} m case'.format(r))
            plt.xlabel('poloidal index')
            plt.title('neutral density in different poloidal index')
            plt.legend()
        
        
        elif self.withshift == False and self.withseries == True:
            self.load_output_data(param= 'NeuDen')
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
            
            efold_dic = {}
            delta_dic = {}
            opq_dic = {}
            neu_den_dic = {}
            ne_ped_dic = {}
            
            for aa in self.data['dircomp']['Attempt'].keys():
                i = 0
                ln = len(pol_list)
                efold = np.zeros(ln)
                delta = np.zeros(ln)
                opq = np.zeros(ln)
                pol_loc = np.zeros(ln)
                neu_den = np.zeros(ln)
                ne_ped = np.zeros(ln)
                for k in self.data['poloidal_index']:
                    psi = self.data['psi']['psi_{}_val'.format(k)]
                    # dsa_pol_loc = self.data['psi']['dsa_{}_val'.format(k)]
                    # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
                    # SEP = int(self.data['DefaultSettings']['SEP'])
                    
                    pol_in = int(k) + 1
                    Nd = self.data['outputdata']['NeuDen'][aa][:, pol_in]
                    Ne = self.data['outputdata']['Ne'][aa][:, pol_in]
                    Te = self.data['outputdata']['Te'][aa][:, pol_in]
                    
                    rd = fm.Opacity_calculator_psi(psi, Ne, Te, Nd)
                    
                    efold[i] = rd['efold_length']
                    delta[i] = rd['pedestal_width']
                    opq[i] = rd['dimensionless_opaqueness']
                    neu_den[i] = max(rd['exp_fit_in_width'])
                    ne_ped[i] = rd['electron_density_pedestal']
                    pol_loc[i] = int(k)
                    i = i + 1
                efold_dic[aa] = efold
                delta_dic[aa] = delta
                opq_dic[aa] = opq
                neu_den_dic[aa] = neu_den
                ne_ped_dic[aa] = ne_ped
            
            result = {'efold_length': efold_dic, 'pedestal_width': delta_dic,
                      'dimensionless_opaqueness': opq_dic}
            self.data['opacity_poloidal'] = result
            
            
            log_flag = False
            
            density_dic = {}
            for k in self.data['dircomp']['Attempt'].keys():
                kk = float(k)*pow(10, 19)
                density_dic[k] = kk
            
            plt.figure(1)
            if log_flag:
                plt.yscale('log')
            else:
                pass
            for aa in self.data['dircomp']['Attempt'].keys():
                plt.plot(pol_loc, efold_dic[aa],'o-', 
                         label= 'Neutral penetration length and electron density is {}'.format(density_dic[aa]))
            if max(pol_loc) > 46:
                plt.axvline(x= 46, color='black',lw=3, ls='--')
            else:
                pass
            plt.xlabel('poloidal index')
            plt.ylabel('$\lambda_{n_D}$: [m]')
            plt.title('Neutral penetration length with different core electron density')
            plt.legend()
            
            plt.figure(2)
            if log_flag:
                plt.yscale('log')
            else:
                pass
            if max(pol_loc) > 46:
                plt.axvline(x= 46, color='black',lw=3, ls='--')
            else:
                pass
            for ab in self.data['dircomp']['Attempt'].keys():
                plt.plot(pol_loc, opq_dic[ab],'o-',
                         label= 'dimensionless opaqueness and electron density is {}'.format(density_dic[ab]))
            plt.xlabel('poloidal index')
            plt.title('dimensionless opaqueness with different core electron density')
            plt.legend()
            
            plt.figure(3)
            if log_flag:
                plt.yscale('log')
            else:
                pass
            if max(pol_loc) > 46:
                plt.axvline(x= 46, color='black',lw=3, ls='--')
            else:
                pass
            for ac in self.data['dircomp']['Attempt'].keys():
                plt.plot(pol_loc, neu_den_dic[ac],'o-', 
                         label= 'neutral density and electron density is {}'.format(density_dic[ac]))
            plt.xlabel('poloidal index')
            plt.title('neutral density with different core electron density')
            plt.legend()
            
            plt.figure(4)
            if log_flag:
                plt.yscale('log')
            else:
                pass
            for ad in self.data['dircomp']['Attempt'].keys():
                plt.plot(pol_loc, ne_ped_dic[ad],'o-', 
                         label= 'electron pedestal density and electron density is {}'.format(density_dic[ad]))
            plt.xlabel('poloidal index')
            plt.title('electron pedestal density with different core electron density')
            plt.legend()
            
            plt.figure(5)
            if log_flag:
                plt.yscale('log')
            else:
                pass
            for ae in self.data['dircomp']['Attempt'].keys():
                plt.plot(pol_loc, delta_dic[ae],'o-', 
                         label= 'pedestal width and electron density is {}'.format(density_dic[ae]))
            plt.xlabel('poloidal index')
            plt.title('pedestal width with different core electron density')
            plt.legend()
            
        
        else:
            print('more work needs to be done')
            
            
    def Opacity_study_radial_plot_psi(self, pol_loc):
        
        if self.withshift == False and self.withseries == False:
            self.load_output_data(param= 'NeuDen')
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
            
            psi = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 2]
            dsa_pol_loc = self.data['psi']['dsa_{}_val'.format(pol_loc)]
            # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
            SEP = int(self.data['DefaultSettings']['SEP'])
            
            pol_index = int(pol_loc) + 1
            Nd = self.data['outputdata']['NeuDen'][:, pol_index]
            Ne = self.data['outputdata']['Ne'][:, pol_index]
            Te = self.data['outputdata']['Te'][:, pol_index]
            
            
            result_dic = fm.Opacity_calculator_psi(psi, Ne, Te, Nd)
            
            tanh_ne_fit = result_dic['tanh_ne_fit']
            tanh_te_fit = result_dic['tanh_te_fit']
            exp_an_fit = result_dic['exp_fit_in_width']
            dn = result_dic['pedestal_width']
            dtn = result_dic['temperature_pedestal_width']
            efold = result_dic['efold_length']
            opq = result_dic['dimensionless_opaqueness']
            psi_cut = result_dic['psi_cut']
            
            
            
            self.data['opacity_study'] = result_dic
            
            
            x = [-efold + max(psi_cut), max(psi_cut)]
            y = [max(exp_an_fit), max(exp_an_fit)]
            xd = [-dn + 1, dn + 1]
            yd = [tanh_ne_fit[SEP], tanh_ne_fit[SEP]]
            xt = [-dtn + 1, dtn + 1]
            yt = [tanh_te_fit[SEP], tanh_te_fit[SEP]]
            
            
            plt.figure(1)
            plt.plot(psi, Nd,'o-', color = 'green', label= 'solps neutral density')
            # plt.plot(psi_RGI, Nd,'o-', color = 'b', label= 'RGI_solps neutral density')
            plt.plot(psi_cut, exp_an_fit, color='r',lw= 5, label= 'exponential fit')
            plt.axvline(x=max(psi_cut), color='orange',lw=3)
            plt.plot(x,y, color='orange', lw=3, label= 'Neutral penetration length [m]: $\lambda_{n_D}$')
            plt.axvline(x=-efold + max(psi_cut), color='orange',lw=3)
            plt.axvline(x=dn + 1, color='black',lw=3, ls='--')
            plt.axvline(x=-dn + 1, color='black',lw=3, ls='--')
            plt.text(0.1, 5*pow(10, 16), 'dimensionless opaqueness: {}'.format(opq))
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
            plt.ylabel(self.data['Parameter']['NeuDen'])
            plt.title('Neutral density with fits')
            # plt.title(plot_dic['an3da.last10'][0],fontdict={"family":"Calibri","size": 20})
            plt.legend()
                
            # plt.subplot(211, sharex= ax1)
            plt.figure(2)
            # plt.plot(psi_xport, Ne,'o-', color = 'r', label= 'solps_electron density')
            plt.plot(psi, Ne,'o-', color = 'b', label= 'solps electron density')
            # plt.plot(fitdsa, cutNe,'o-', color = 'g', label= 'experiment electron density')
            plt.plot(psi, tanh_ne_fit, color='r',lw= 3, label= 'tanh fit')
            plt.plot(xd, yd, color='black', lw=3, label= 'Pedestal width [m]: $\Delta n_e$')
            plt.axvline(x=dn + 1, color='black',lw=3)
            plt.axvline(x=-dn + 1, color='black',lw=3)
            plt.axvline(x=max(psi_cut), color='orange',lw=3, ls='--')
            plt.axvline(x=-efold + max(psi_cut), color='orange',lw=3, ls='--')
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
            plt.ylabel(self.data['Parameter']['Ne'])
            plt.title('Electron density with fits')
            # plt.title(plot_dic['ne3da.last10'][0],fontdict={"family":"Calibri","size": 20})
            plt.legend()
            
            plt.figure(3)
            # plt.plot(psi_xport, Ne,'o-', color = 'r', label= 'solps_electron density')
            plt.plot(psi, Te,'o-', color = 'b', label= 'solps electron tempurature')
            # plt.plot(fitdsa, cutNe,'o-', color = 'g', label= 'experiment electron density')
            plt.plot(psi, tanh_te_fit, color='r',lw= 3, label= 'tanh fit')
            plt.plot(xt, yt, color='black', lw=3, label= 'temperature pedestal width [m]: $\Delta n_e$')
            plt.axvline(x=dtn + 1, color='black',lw=3)
            plt.axvline(x=-dtn + 1, color='black',lw=3)
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
            plt.ylabel(self.data['Parameter']['Te'])
            plt.title('Electron temperature with fits')
            plt.legend()
            
            
            plt.show()
        
        else:
            print('More works for psi are needed')
        
        
    
    def Opacity_study_radial_plot(self, pol_loc):
           
        if self.withshift == False and self.withseries == False:
            self.load_output_data(param= 'NeuDen')
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
            
            psi = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 2]
            dsa_pol_loc = self.data['psi']['dsa_{}_val'.format(pol_loc)]
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
            
            
            fit_tanh_dic = fm.tanh_dsa_fit(dsa_pol_loc, Ne, Te)
            tanh_ne_fit = fit_tanh_dic['tanh_ne_fit']
            dn = fit_tanh_dic['popt_ne'][1]
            tanh_te_fit = fit_tanh_dic['tanh_te_fit']
            dtn = fit_tanh_dic['popt_te'][1]
            
            
            dsa_cut = []
            an_cut = []
            for j in range(len(dsa_pol_loc)):
                if dsa_pol_loc[j] <= dn and dsa_pol_loc[j] >= -dn:
                    dsa_cut.append(dsa_pol_loc[j])
                    an_cut.append(Nd[j])
            
            dsa_cut = np.asarray(dsa_cut)
            an_cut = np.asarray(an_cut)
            
            fit_exp_dic = fm.exp_dsa_fit(dsa_cut, an_cut)
            exp_an_fit = fit_exp_dic['exp_an_fit']
            efold = 1/fit_exp_dic['popt_an'][1]
            
            
            opq = 2*dn/efold
            print(opq)
            
            result_dic = {'tanh_fit': tanh_ne_fit, 'exp_fit': exp_an_fit,
                          'pedestal_width': dn, 'efold_length': efold,
                          'dimensionless_opaqueness': opq}
            
            self.data['opacity_study'] = result_dic
            
            
            x = [-efold + max(dsa_cut), max(dsa_cut)]
            y = [max(exp_an_fit), max(exp_an_fit)]
            xd = [-dn, dn]
            yd = [tanh_ne_fit[SEP], tanh_ne_fit[SEP]]
            xt = [-dtn, dtn]
            yt = [tanh_te_fit[SEP], tanh_te_fit[SEP]]
            
            
            plt.figure(1)
            plt.plot(dsa_pol_loc, Nd,'o-', color = 'green', label= 'solps neutral density')
            # plt.plot(psi_RGI, Nd,'o-', color = 'b', label= 'RGI_solps neutral density')
            plt.plot(dsa_cut, exp_an_fit, color='r',lw= 5, label= 'exponential fit')
            plt.axvline(x=max(dsa_cut), color='orange',lw=3)
            plt.plot(x,y, color='orange', lw=3, label= 'Neutral penetration length [m]: $\lambda_{n_D}$')
            plt.axvline(x=-efold + max(dsa_cut), color='orange',lw=3)
            plt.axvline(x=dn, color='black',lw=3, ls='--')
            plt.axvline(x=-dn, color='black',lw=3, ls='--')
            plt.text(0.1, 5*pow(10, 16), 'dimensionless opaqueness: {}'.format(opq))
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
            plt.ylabel(self.data['Parameter']['NeuDen'])
            plt.title('Neutral density with fits')
            # plt.title(plot_dic['an3da.last10'][0],fontdict={"family":"Calibri","size": 20})
            plt.legend()
                
            # plt.subplot(211, sharex= ax1)
            plt.figure(2)
            # plt.plot(psi_xport, Ne,'o-', color = 'r', label= 'solps_electron density')
            plt.plot(dsa_pol_loc, Ne,'o-', color = 'b', label= 'solps electron density')
            # plt.plot(fitdsa, cutNe,'o-', color = 'g', label= 'experiment electron density')
            plt.plot(dsa_pol_loc, tanh_ne_fit, color='r',lw= 3, label= 'tanh fit')
            plt.plot(xd, yd, color='black', lw=3, label= 'Pedestal width [m]: $\Delta n_e$')
            plt.axvline(x=dn, color='black',lw=3)
            plt.axvline(x=-dn, color='black',lw=3)
            plt.axvline(x=max(dsa_cut), color='orange',lw=3, ls='--')
            plt.axvline(x=-efold + max(dsa_cut), color='orange',lw=3, ls='--')
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
            plt.ylabel(self.data['Parameter']['Ne'])
            plt.title('Electron density with fits')
            # plt.title(plot_dic['ne3da.last10'][0],fontdict={"family":"Calibri","size": 20})
            plt.legend()
            
            plt.figure(3)
            # plt.plot(psi_xport, Ne,'o-', color = 'r', label= 'solps_electron density')
            plt.plot(dsa_pol_loc, Te,'o-', color = 'b', label= 'solps electron tempurature')
            # plt.plot(fitdsa, cutNe,'o-', color = 'g', label= 'experiment electron density')
            plt.plot(dsa_pol_loc, tanh_te_fit, color='r',lw= 3, label= 'tanh fit')
            plt.plot(xt, yt, color='black', lw=3, label= 'temperature pedestal width [m]: $\Delta n_e$')
            plt.axvline(x=dtn, color='black',lw=3)
            plt.axvline(x=-dtn, color='black',lw=3)
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
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
            
            
            
        
