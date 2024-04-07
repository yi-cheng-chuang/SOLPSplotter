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
        
        self.Publish = DefaultSettings['Publish']
        self.data['DefaultSettings']['Publish'] = self.Publish
            
    
    def set_plot(self):
        if self.Publish == 'b2plottersetting':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 3, markersize= 7)
            plt.rcParams.update({'font.size': 20})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
            # plt.rcParams["text.usetex"] = True
  
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
        plt.plot(x_coord[1:37], Nd,'-', color = 'green', label= 'solps neutral density')
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
        plt.xlabel('Normalized flux coordinate $\psi_N$')
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
            # Nd = self.data['radial_fit_data']['NeuDen'][:, pol_index]
            # Ne = self.data['radial_fit_data']['Ne'][:, pol_index]
            # Te = self.data['radial_fit_data']['Te'][:, pol_index]
            Nd = self.data['radial_fit_data']['NeuDen']
            Ne = self.data['radial_fit_data']['Ne']
            Te = self.data['radial_fit_data']['Te']
            SEP = int(self.data['DefaultSettings']['sep_index_dsa'])
            psi = self.data['psi']['psi_{}_val'.format(pol_loc[0])][:, 2]
            
            
            result_dic = self.data['radial_fit_data'] | self.data['opacity_poloidal']
            
            
            P = self.data['Parameter']
            self.opacity_radial_method(result_dic = result_dic, SEP = SEP, 
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
        
    
    def plot_iout_radial_divertor(self, quant, log_scale):
        
        
        
        if self.withshift == False and self.withseries == False:
            
            print('plot_iout_radial function is in prepare ...')
        
        elif self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
                       
            div_side_list = ['inner target', 'outer target', 
                             'inner x point boundary', 'outer x point boundary']
            # qu = self.load_iout(filename = dname, simple_quant = quant)
            
            for side in div_side_list:
                
                plt.figure(figsize=(7,7))
                
                for aa in self.data['dircomp']['multi_shift']:
                    
                    if side == 'inner target':
                        psiN = self.data['psi']['psival'][aa][1:37, 1]
                        plot_data = self.data['iout_data'][quant][aa][:, 1]
                        
                    elif side == 'outer target':
                        
                        psiN = self.data['psi']['psival'][aa][1:37, -2]
                        plot_data = self.data['iout_data'][quant][aa][:, -2]
                    
                    elif side == 'inner x point boundary':
                        
                        inner_index = self.data['b2fgeo'][aa]['leftcut'][0] - 1
                        
                        psiN = self.data['psi']['psival'][aa][1:37, inner_index]
                        plot_data = self.data['iout_data'][quant][aa][:, inner_index]
                    
                    elif side == 'outer x point boundary':
                        
                        
                        outer_index = self.data['b2fgeo'][aa]['rightcut'][0] + 1
                        
                        psiN = self.data['psi']['psival'][aa][1:37, outer_index]
                        plot_data = self.data['iout_data'][quant][aa][:, outer_index]
                    
                    
                    
                    if log_scale:
                        plot_data = np.abs(plot_data)
                        plt.yscale('log')
                    else:
                        pass
                    
                        
                    plt.plot(psiN, plot_data, '-', color = color_dic[aa], 
                                 label = '{}'.format(A_dic[aa]))
                    plt.legend()
                
                plt.xlabel('psiN')
                plt.title('{} at {}'.format(quant, side))
                plt.show()
                
        
        elif self.withshift == False and self.withseries == True:
            
            print('plot_iout_radial function is in prepare ...')
    
    
    def plot_iout_radial_xpoint(self, quant, log_scale):
        
        
        
        if self.withshift == False and self.withseries == False:
            
            print('plot_iout_radial function is in prepare ...')
        
        elif self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
                       
            div_side_list = ['inner x top', 'outer x top', 
                             'inner x down', 'outer x down']
            # qu = self.load_iout(filename = dname, simple_quant = quant)
            
            for side in div_side_list:
                
                plt.figure(figsize=(7,7))
                
                for aa in self.data['dircomp']['multi_shift']:
                    
                    if side == 'inner x top':
                        
                        intop_index = self.data['b2fgeo'][aa]['leftcut'][0] + 2
                        
                        psiN = self.data['psi']['psival'][aa][1:37, intop_index]
                        plot_data = self.data['iout_data'][quant][aa][:, intop_index]
                        
                    elif side == 'outer x top':
                        
                        outop_index = self.data['b2fgeo'][aa]['rightcut'][0] - 2
                        
                        psiN = self.data['psi']['psival'][aa][1:37, outop_index]
                        plot_data = self.data['iout_data'][quant][aa][:, outop_index]
                    
                    elif side == 'inner x down':
                        
                        inner_index = self.data['b2fgeo'][aa]['leftcut'][0] - 1
                        
                        psiN = self.data['psi']['psival'][aa][1:37, inner_index]
                        plot_data = self.data['iout_data'][quant][aa][:, inner_index]
                    
                    elif side == 'outer x down':
                        
                        
                        outer_index = self.data['b2fgeo'][aa]['rightcut'][0] + 1
                        
                        psiN = self.data['psi']['psival'][aa][1:37, outer_index]
                        plot_data = self.data['iout_data'][quant][aa][:, outer_index]
                    
                    
                    
                    if log_scale:
                        plot_data = np.abs(plot_data)
                        plt.yscale('log')
                    else:
                        pass
                    
                        
                    plt.plot(psiN, plot_data, '-', color = color_dic[aa], 
                                 label = '{}'.format(A_dic[aa]))
                    plt.legend()
                
                plt.xlabel('psiN')
                plt.title('{} at {}'.format(quant, side))
                plt.show()
                
        
        elif self.withshift == False and self.withseries == True:
            
            print('plot_iout_radial function is in prepare ...')  
    
            
    def plot_radial_xpoint(self):
        
        
        
        if self.withshift == False and self.withseries == False:
            
            print('plot_iout_radial function is in prepare ...')
        
        elif self.withshift == True and self.withseries == False:
            
            color_dic = {'inner x top': 'red', 'outer x top': 'orange', 
                'inner x down': 'green', 'outer x down': 'blue'}
            
            
            
            div_side_list = ['inner x top', 'outer x top', 
                             'inner x down', 'outer x down']
            # qu = self.load_iout(filename = dname, simple_quant = quant)
            
            plt.figure(figsize=(7,7))
            
            for side in div_side_list:
                
                aa = 'org'
                    
                if side == 'inner x top':
                    
                    intop_index = self.data['b2fgeo'][aa]['leftcut'][0] + 2
                    
                    rloc = self.data['grid']['RadLoc'][aa][:, intop_index]
                    zloc = self.data['grid']['VertLoc'][aa][:, intop_index]
                    
                elif side == 'outer x top':
                    
                    outop_index = self.data['b2fgeo'][aa]['rightcut'][0] - 2
                    
                    rloc = self.data['grid']['RadLoc'][aa][:, outop_index]
                    zloc = self.data['grid']['VertLoc'][aa][:, outop_index]
                
                elif side == 'inner x down':
                    
                    inner_index = self.data['b2fgeo'][aa]['leftcut'][0] - 1
                    
                    rloc = self.data['grid']['RadLoc'][aa][:, inner_index]
                    zloc = self.data['grid']['VertLoc'][aa][:, inner_index]
                
                elif side == 'outer x down':
                    
                    
                    outer_index = self.data['b2fgeo'][aa]['rightcut'][0] + 1
                    
                    rloc = self.data['grid']['RadLoc'][aa][:, outer_index]
                    zloc = self.data['grid']['VertLoc'][aa][:, outer_index]
                
                
                    
                plt.plot(rloc, zloc, '-', color = color_dic[side], 
                             label = '{}'.format(side))
                plt.legend()
                
                plt.xlabel('R')
                plt.title('RZ location')
                plt.show()
                
        
        elif self.withshift == False and self.withseries == True:
            
            print('plot_iout_radial function is in prepare ...')
    
    
    def plot_iout_radial_percent(self, quant, log_scale):
        
               
        if self.withshift == False and self.withseries == False:
            
            print('plot_iout_radial function is in prepare ...')
        
        elif self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
                       
            div_side_list = ['inner target', 'outer target', 
                             'inner x point boundary', 'outer x point boundary']
            # qu = self.load_iout(filename = dname, simple_quant = quant)
            
            for side in div_side_list:
                
                plt.figure(figsize=(7,7))
                
                for aa in self.data['dircomp']['multi_shift']:
                    if aa == 'org':
                        pass
                    else:
                        
                        
                        if side == 'inner target':
                            psiN = self.data['psi']['psival'][aa][1:37, 1]
                            plot_data = self.data['iout_data'][quant][aa][:, 1]
                            std_data = self.data['iout_data'][quant]['org'][:, 1]
                            
                        elif side == 'outer target':
                            
                            psiN = self.data['psi']['psival'][aa][1:37, -2]
                            plot_data = self.data['iout_data'][quant][aa][:, -2]
                            std_data = self.data['iout_data'][quant]['org'][:, -2]
                        
                        elif side == 'inner x point boundary':
                            
                            inner_index = self.data['b2fgeo'][aa]['leftcut'][0] - 1
                            
                            psiN = self.data['psi']['psival'][aa][1:37, inner_index]
                            plot_data = self.data['iout_data'][quant][aa][:, inner_index]
                            std_data = self.data['iout_data'][quant]['org'][:, inner_index]
                        
                        elif side == 'outer x point boundary':
                            
                            
                            outer_index = self.data['b2fgeo'][aa]['rightcut'][0] + 1
                            
                            psiN = self.data['psi']['psival'][aa][1:37, outer_index]
                            plot_data = self.data['iout_data'][quant][aa][:, outer_index]
                            std_data = self.data['iout_data'][quant]['org'][:, outer_index]
                        
                        
                        if log_scale:
                            plot_data = np.abs(plot_data)
                            plt.yscale('log')
                        else:
                            pass
                        
                        
                        dat_diff = plot_data - std_data
                        dat_percent = np.divide(dat_diff, std_data)*100
                        
                            
                        plt.plot(psiN, dat_percent, '-', color = color_dic[aa], 
                                     label = '{}'.format(A_dic[aa]))
                        plt.legend()
                    
                
                plt.xlabel('psiN')
                plt.title('{} at {}'.format(quant, side))
                plt.show()
                        
                    
                    
                
        
        elif self.withshift == False and self.withseries == True:
            
            print('plot_iout_radial function is in prepare ...')
                  
        
                
    def divertor_te(self):  
        
        b2fstate = self.data['b2fstate']

        
        plt.figure(figsize=(7,7))

        for aa in self.data['dircomp']['multi_shift']:
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            rcood = self.data['psi']['psival'][aa][:, 1]
            Te_J = b2fstate[aa]['te'].transpose()
            ev = 1.6021766339999999 * pow(10, -19)
            te_pro = Te_J / ev
            te = te_pro[:, 0]
            plt.plot(rcood, te, '-', color = color_dic[aa], 
                     label= 'aspect ratio = {}'.format(A_dic[aa]))
            plt.title('inner target electron temperature')
            plt.xlabel('psiN')
            plt.legend()
      
        
        plt.figure(figsize=(7,7))

        for aa in self.data['dircomp']['multi_shift']:
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            rcood = self.data['psi']['psival'][aa][:, -2]
            Te_J = b2fstate[aa]['te'].transpose()
            ev = 1.6021766339999999 * pow(10, -19)
            te_pro = Te_J / ev
            te = te_pro[:, -1]
            plt.plot(rcood, te, '-', color = color_dic[aa], 
                     label= 'aspect ratio = {}'.format(A_dic[aa]))
            plt.title('outer target electron temperature')
            plt.xlabel('psiN')
            plt.legend()
      
      
            
    
    
    def plot_divertor_radial(self, r_coord, data, log_scale, quant, div_side):
        
        plt.figure(figsize=(7,7))
        if log_scale:
            plt.yscale('log')
        plt.plot(r_coord, data, '-', color = 'r')
        plt.xlabel('psiN')
        plt.title('{} at {}'.format(quant, div_side))
        plt.show()
    
    
    
    def plot_all_radial(self, separate):
        
        if self.withshift == False and self.withseries == False:
            
            """
            
            # if self.data['outputdata'].any() == None or self.data['outputdata']['Te'].any() == None:
            if 'Ne' and 'Te' and 'NeuDen' in self.data['outputdata']:
                pass
            else:
                self.load_output_data(param= 'Ne')
                self.load_output_data(param= 'Te')
                self.load_output_data(param= 'NeuDen')
            
            ne_pro = self.data['outputdata']['Ne']
            te_pro = self.data['outputdata']['Te']
            neu_pro = self.data['outputdata']['NeuDen']
            
            """
            
            b2fstate = self.data['b2fstate']
            
            ne_pro = b2fstate['ne'].transpose()
            Te_J = b2fstate['te'].transpose()
            ev = 1.6021766339999999 * pow(10, -19)
            te_pro = Te_J / ev
            
            data = self.data['ft44']['dab2']
            neu_pro = np.transpose(data[:, :, 0])
            
            
            core_ne_pro = ne_pro[:, 25:71]
            core_te_pro = te_pro[:, 25:71]
            core_neu_pro = neu_pro[:, 25:71]
            
            innerleg_ne = ne_pro[:, :25]
            innerleg_te = te_pro[:, :25]
            innerleg_neu = neu_pro[:, :25]
            
            outerleg_ne = ne_pro[:, 73:96]
            outerleg_te = te_pro[:, 73:96]
            outerleg_neu = neu_pro[:, 73:96]
            

            mean_core_ne = np.mean(core_ne_pro, axis=1)
            std_core_ne = np.std(core_ne_pro, axis=1)
            # print(std_core_ne)
            
            mean_core_te = np.mean(core_te_pro, axis=1)
            std_core_te = np.std(core_te_pro, axis=1)
            
            mean_core_neu = np.mean(core_neu_pro, axis=1)
            std_core_neu = np.std(core_neu_pro, axis=1)
            
            
            
            mean_innerleg_ne = np.mean(innerleg_ne, axis=1)
            std_innerleg_ne = np.std(innerleg_ne, axis=1)
            # print(std_innerleg_ne)
            mean_innerleg_te = np.mean(innerleg_te, axis=1)
            std_innerleg_te = np.std(innerleg_te, axis=1)
            
            mean_innerleg_neu = np.mean(innerleg_neu, axis=1)
            std_innerleg_neu = np.std(innerleg_neu, axis=1)
            
            
            
            mean_outerleg_ne = np.mean(outerleg_ne, axis=1)
            std_outerleg_ne = np.std(outerleg_ne, axis=1)
            # print(std_outerleg_ne)
            mean_outerleg_te = np.mean(outerleg_te, axis=1)
            std_outerleg_te = np.std(outerleg_te, axis=1)
            
            mean_outerleg_neu = np.mean(outerleg_neu, axis=1)
            std_outerleg_neu = np.std(outerleg_neu, axis=1)
            
            
            
            psiN = self.data['experimental_fit']['psiN']
            ne = self.data['experimental_fit']['ne']*pow(10, 20)
            te = self.data['experimental_fit']['te']*pow(10, 3)
            
            exp = self.data['ExpDict']
            psi = exp['psi_normal']
            
            
            psi = []
            exp_ne = []
            ne_er = []
            exp_te = []
            te_er = []
            for ep in range(len(exp['psi_normal'])):
                
                if  exp['psi_normal'][ep] >= min(psiN):
                    psi.append(exp['psi_normal'][ep])
                    exp_ne.append(exp['electron_density(10^20/m^3)'][ep]*pow(10, 20))
                    ne_er.append(exp['density error(10^20/m^3)'][ep]*pow(10, 20))
                    exp_te.append(exp['electron_temperature(KeV)'][ep]*pow(10, 3))
                    te_er.append(exp['temperature error(10^20/m^3)'][ep]*pow(10, 3))
                    
                    
            # exp_ne = exp['electron_density(10^20/m^3)']*pow(10, 20)
            # ne_er = exp['density error(10^20/m^3)']*pow(10, 20)
            # exp_te = exp['electron_temperature(KeV)']*pow(10, 3)
            # te_er = exp['temperature error(10^20/m^3)']*pow(10, 3)
            
            'core'
            
        
            
            if separate:
                plt.figure(figsize=(7,7))
                plt.yscale('log')
                plt.errorbar(psiN, mean_core_ne, yerr= std_core_ne, fmt = '-', color = 'g', label= 'ne_solps')
                plt.errorbar(psi, exp_ne, yerr= ne_er, fmt = 'o', color = 'b', label= 'ne_exp')
                # plt.plot(psiN, ne, 'o', color = 'r', label= 'ne_exp_fit')
                plt.xlabel('psiN')
                plt.title('electron density with experimental fit')
                plt.legend()
                
                
                plt.figure(figsize=(7,7))
                plt.yscale('log')
                plt.errorbar(psiN, mean_core_te, yerr= std_core_te, fmt = '-', color = 'g', label= 'te_solps')
                plt.errorbar(psi, exp_te, yerr= te_er, fmt = 'o', color = 'b', label= 'te_exp')
                # plt.plot(psiN, ne, 'o', color = 'r', label= 'ne_exp_fit')
                plt.xlabel('psiN')
                plt.title('electron temperature with experimental fit')
                plt.legend()
            else:
                fig, axs = plt.subplots(1, 2)
                
                axs[0].set_yscale('log')
                axs[0].errorbar(psiN, mean_core_ne, yerr= std_core_ne, fmt = '-', color = 'g', label= 'ne_solps')
                axs[0].errorbar(psi, exp_ne, yerr= ne_er, fmt = 'o', color = 'b', label= 'ne_exp')
                # plt.plot(psiN, ne, 'o', color = 'r', label= 'ne_exp_fit')
                axs[0].set_xlabel('Normalized flux coordinate $\psi_N$')
                axs[0].set_title('electron density with experimental fit')
                axs[0].legend()
                
                axs[1].set_yscale('log')
                axs[1].errorbar(psiN, mean_core_te, yerr= std_core_te, fmt = '-', color = 'g', label= 'te_solps')
                axs[1].errorbar(psi, exp_te, yerr= te_er, fmt = 'o', color = 'b', label= 'te_exp')
                # plt.plot(psiN, ne, 'o', color = 'r', label= 'ne_exp_fit')
                axs[1].set_xlabel('Normalized flux coordinate $\psi_N$')
                axs[1].set_title('electron temperature with experimental fit')
                axs[1].legend()
                
            
            
            
            print('the shape of mean_core_neu is {}'.format(np.shape(mean_core_neu)))
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN[1:37], mean_core_neu, yerr= std_core_neu, fmt = 'o', color = 'g', label= 'Neuden_solps')
            plt.xlabel('psiN')
            plt.title('Neutral density')
            plt.legend()
            
            'inner leg'
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN, mean_innerleg_ne, yerr= std_innerleg_ne, fmt = 'o', color = 'g', label= 'ne_solps')
            plt.xlabel('psiN')
            plt.title('inner leg electron density')
            plt.legend()
            
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN, mean_innerleg_te, yerr= std_innerleg_te, fmt = 'o', color = 'g', label= 'te_solps')
            plt.xlabel('psiN')
            plt.title('inner leg electron temperature')
            plt.legend()
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN[1:37], mean_innerleg_neu, yerr= std_innerleg_neu, fmt = 'o', color = 'g', label= 'neuden_solps')
            plt.xlabel('psiN')
            plt.title('inner leg neutral density')
            plt.legend()
            
            
            'outerleg'
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN, mean_outerleg_ne, yerr= std_outerleg_ne, fmt = 'o', color = 'g', label= 'ne_solps')
            plt.xlabel('psiN')
            plt.title('outer leg electron density')
            plt.legend()
            
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN, mean_outerleg_te, yerr= std_outerleg_te, fmt = 'o', color = 'g', label= 'te_solps')
            plt.xlabel('psiN')
            plt.title('outer leg electron temperature')
            plt.legend()
            
            plt.figure(figsize=(7,7))
            plt.yscale('log')
            plt.errorbar(psiN[1:37], mean_outerleg_neu, yerr= std_outerleg_neu, fmt = 'o', color = 'g', label= 'neuden_solps')
            plt.xlabel('psiN')
            plt.title('outer leg neutral density')
            plt.legend()
        
        

        



        
        
        
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




    
    
    