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
    

    
    def plot_Ne_NeuDen_test(self, pol_loc):
        
        if self.withshift == False:
            self.load_output_data(param= 'NeuDen')
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
            
            psi = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 2]
            dsa_pol_loc = self.data['psi']['dsa_{}_val'.format(pol_loc)]
            # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
            
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
            
            dsa_cut = []
            an_cut = []
            for j in range(len(dsa_pol_loc)):
                if dsa_pol_loc[j] <= 0:
                    dsa_cut.append(dsa_pol_loc[j])
                    an_cut.append(Nd[j])
            
            dsa_cut = np.asarray(dsa_cut)
            an_cut = np.asarray(an_cut)
            
            
            fit_exp_dic = fm.exp_dsa_fit(dsa_cut, an_cut)
            exp_an_fit = fit_exp_dic['exp_an_fit']
            efold = 1/fit_exp_dic['popt_an'][1]
            
            
            fit_tanh_dic = fm.tanh_dsa_fit(dsa_pol_loc, Ne, Te)
            tanh_ne_fit = fit_tanh_dic['tanh_ne_fit']
            dn = fit_tanh_dic['popt_ne'][1]
            
            
            opq = 2*dn/efold
            print(opq)
            
            x = [-efold + max(dsa_cut), max(dsa_cut)]
            y = [max(exp_an_fit), max(exp_an_fit)]
            xd = [-dn, dn]
            yd = [cutNe[109], cutNe[109]]
            
            
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
            plt.plot(dsa_pol_loc, tanh_ne_fit, color='r',lw= 3, label= 'exponential fit')
            plt.plot(xd, yd, color='black', lw=3, label= 'Pedestal width [m]: $\Delta n_e$')
            plt.axvline(x=dn, color='black',lw=3)
            plt.axvline(x=-dn, color='black',lw=3)
            plt.axvline(x=0, color='orange',lw=3, ls='--')
            plt.axvline(x=-efold, color='orange',lw=3, ls='--')
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
            plt.ylabel(self.data['Parameter']['Ne'])
            plt.title('Electron density with fits')
            # plt.title(plot_dic['ne3da.last10'][0],fontdict={"family":"Calibri","size": 20})
            plt.legend()
            
            plt.show()
            
            return opq
        elif self.withshift == True:
            print('please use plot_Ne_NeuDen_withshift, this is for code testing')
            
        else:
            print('plot_Ne_NeuDen_test function has a bug')
        
        
        
    def plot_Ne_NeuDen_single(self, pol_loc):
           
        if self.withshift == False:
            self.load_output_data(param= 'NeuDen')
            self.load_output_data(param= 'Ne')
            self.load_output_data(param= 'Te')
            
            psi = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 2]
            dsa_pol_loc = self.data['psi']['dsa_{}_val'.format(pol_loc)]
            # psi_RGI = self.data['psi']['psi_{}_val'.format(pol_loc)][:, 0]
            
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
            
            x = [-efold + max(dsa_cut), max(dsa_cut)]
            y = [max(exp_an_fit), max(exp_an_fit)]
            xd = [-dn, dn]
            yd = [cutNe[109], cutNe[109]]
            
            
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
            plt.plot(dsa_pol_loc, tanh_ne_fit, color='r',lw= 3, label= 'exponential fit')
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
            
            plt.show()
            
            # dat_dic = {'tanh_fit': tanh_ne_fit}
            
            # self.data['fitplot'] = dat_dic
            
            return opq
        elif self.withshift == True:
            print('please use plot_Ne_NeuDen_withshift, this is for plotting single case')
            
        else:
            print('plot_Ne_NeuDen_single function has a bug')
            
            
    def plot_Ne_NeuDen_withshift(self, pol_loc):    
        
        if self.withshift == True:
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
                psi = self.data['psi']['psi_{}_val'.format(pol_loc)][aa][:, 2]
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
                # print(opq)
            
            
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
            
     
        elif self.withshift == False:
            print('please use plot_Ne_NeuDen_single, this is for shift cases')
            
        else:
            print('plot_Ne_NeuDen_withshift function has a bug')
