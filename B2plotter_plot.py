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
    def __init__(self, DEV, withshift, DefaultSettings, loadDS, Parameters, Publish):
        load_data.__init__(self, DEV, withshift, DefaultSettings, loadDS, Parameters)
        
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
            
            norm_an = an_cut/ max(an_cut)
            
            pn = [1, 250]
            
            popt_an, pcov_an = curve_fit(fm.expfit, dsa_cut, norm_an, pn)
            print(popt_an)
            
            exp_an_fit = fm.expfit(dsa_cut, popt_an[0]*max(an_cut), popt_an[1])
            
            efold = 1/popt_an[1]
            
            
            fit_dic = fm.tanh_dsa_fit(dsa_pol_loc, Ne, Te)
            tanh_ne_fit = fit_dic['tanh_ne_fit']*pow(10, 20)
            dn = fit_dic['popt_ne'][1]
            
            
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
                
                dsa_cut = []
                an_cut = []
                for j in range(len(dsa_pol_loc)):
                    if dsa_pol_loc[j] <= 0:
                        dsa_cut.append(dsa_pol_loc[j])
                        an_cut.append(Nd[j])
                
                dsa_cut = np.asarray(dsa_cut)
                an_cut = np.asarray(an_cut)
                
                norm_an = an_cut/ max(an_cut)
                
                pn = [1, 250]
                
                popt_an, pcov_an = curve_fit(fm.expfit, dsa_cut, norm_an, pn)
                print(popt_an)
                
                exp_an_fit = fm.expfit(dsa_cut, popt_an[0]*max(an_cut), popt_an[1])
                
                efold = 1/popt_an[1]
                
                
                fit_dic = fm.tanh_dsa_fit(dsa_pol_loc, Ne, Te)
                tanh_ne_fit = fit_dic['tanh_ne_fit']*pow(10, 20)
                dn = fit_dic['popt_ne'][1]
                
                
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
            
            
    def plot_Ne_NeuDen_single(self, pol_loc):    
        
        if self.withshift == False:
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
            
            # self.loadmastdata()
            # fitpsi = self.data['fitprofile']['psi_normal']
            # fitNe = self.data['fitprofile']['electron_density(m^(-3))']
            
                # psi_to_dsa_func = interpolate.interp1d(psi, dsa_pol_loc, fill_value = 'extrapolate')
                  
                # cutpsi = []
                # cutNe = []
                # n = len(fitpsi)
                # for i in range(n):
                #     if fitpsi[i] >= min(psi):
                #         cutpsi.append(fitpsi[i])
                #         cutNe.append(fitNe[i])
                
                # cutpsi = np.asarray(cutpsi)
                # cutNe = np.asarray(cutNe)
                
                # fitdsa = psi_to_dsa_func(cutpsi)
                
                dsa_cut = []
                an_cut = []
                for j in range(len(dsa_pol_loc_dic[aa])):
                    if dsa_pol_loc_dic[aa][j] <= 0:
                        dsa_cut.append(dsa_pol_loc_dic[aa][j])
                        an_cut.append(Nd_dic[aa][j])
                
                dsa_cut_dic[aa] = np.asarray(dsa_cut)
                an_cut = np.asarray(an_cut)
                
                norm_an = an_cut/ max(an_cut)
                
                pn = [1, 250]
                
                popt_an, pcov_an = curve_fit(fm.expfit, dsa_cut, norm_an, pn)
                print(popt_an)
                
                exp_an_fit_dic[aa] = fm.expfit(dsa_cut_dic[aa], popt_an[0]*max(an_cut), popt_an[1])
                
                efold_dic[aa] = 1/popt_an[1]
                
                
                fit_dic = fm.tanh_dsa_fit(dsa_pol_loc_dic[aa], Ne_dic[aa], Te_dic[aa])
                tanh_ne_fit_dic[aa] = fit_dic['tanh_ne_fit']*pow(10, 20)
                delta_dic[aa] = fit_dic['popt_ne'][1]
                
                
                opq_dic[aa] = 2*delta_dic[aa]/efold_dic[aa]
                # print(opq)
            
            
            plt.figure(1)
            for i in self.data['dircomp']['multi_shift']:
                plt.plot(dsa_pol_loc_dic[i], Nd_dic[i],'o-', label= 'solps neutral density_{}'.format(i))
                # plt.plot(psi_RGI, Nd,'o-', color = 'b', label= 'RGI_solps neutral density')
                plt.plot(dsa_cut_dic[i], exp_an_fit_dic[i], lw= 5, label= 'exponential fit_{}'.format(i))
            # plt.axvline(x=max(dsa_cut), color='orange',lw=3)
            # plt.plot(x,y, color='orange', lw=3, label= 'Neutral penetration length [m]: $\lambda_{n_D}$')
            # plt.axvline(x=-efold + max(dsa_cut), color='orange',lw=3)
            # plt.axvline(x=dn, color='black',lw=3, ls='--')
            # plt.axvline(x=-dn, color='black',lw=3, ls='--')
            # plt.text(0.1, 5*pow(10, 16), 'dimensionless opaqueness: {}'.format(opq))
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
            plt.ylabel(self.data['Parameter']['NeuDen'])
            plt.title('Neutral density with fits')
            # plt.title(plot_dic['an3da.last10'][0],fontdict={"family":"Calibri","size": 20})
            plt.legend()
                
            # plt.subplot(211, sharex= ax1)
            plt.figure(2)
            # plt.plot(psi_xport, Ne,'o-', color = 'r', label= 'solps_electron density')
            for j in self.data['dircomp']['multi_shift']:
                plt.plot(dsa_pol_loc_dic[i], Ne_dic[i],'o-', label= 'solps electron density_{}'.format(j))
            # plt.plot(fitdsa, cutNe,'o-', color = 'g', label= 'experiment electron density')
            # plt.plot(dsa_pol_loc_dic[i], tanh_ne_fit, color='r',lw= 3, label= 'exponential fit')
            # plt.plot(xd, yd, color='black', lw=3, label= 'Pedestal width [m]: $\Delta n_e$')
            # plt.axvline(x=dn, color='black',lw=3)
            # plt.axvline(x=-dn, color='black',lw=3)
            # plt.axvline(x=0, color='orange',lw=3, ls='--')
            # plt.axvline(x=-efold, color='orange',lw=3, ls='--')
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
            plt.ylabel(self.data['Parameter']['Ne'])
            plt.title('Electron density with fits')
            # plt.title(plot_dic['ne3da.last10'][0],fontdict={"family":"Calibri","size": 20})
            plt.legend()
            
            
            shift_ar = {}
            for k in self.data['dircomp']['multi_shift']:
                kk = float(self.data['dircomp']['shiftdic'][k])
                shift_ar[k] = kk
            
            plt.figure(3)
            for j in self.data['dircomp']['multi_shift']:
                plt.plot(shift_ar[j], opq_dic[j],'o-', color= 'r', label= 'dimensionless opaqueness')
            plt.xlabel('shift: [m]')
            plt.title('opaqueness verses shift distance')
            # plt.title(plot_dic['ne3da.last10'][0],fontdict={"family":"Calibri","size": 20})
            
            plt.figure(4)
            for j in self.data['dircomp']['multi_shift']:
                plt.plot(shift_ar[j], efold_dic[j],'o-', color= 'r', label= 'neutral penetration')
            plt.xlabel('shift: [m]')
            plt.ylabel('neutral penetration length: [m]')
            plt.title('neutral penetration length verses shift distance')
     
            
        else:
            print('We can go to sleep')