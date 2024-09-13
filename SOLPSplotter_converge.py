# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 14:24:08 2024

@author: ychuang
"""


from SOLPSplotter_fit import profile_fit
import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.offsetbox import AnchoredText
from scipy import interpolate


class SOLPS_converge(profile_fit):
    
    def __init__(self, DefaultSettings, loadDS):
        profile_fit.__init__(self, DefaultSettings, loadDS)
            
    
    def set_plot(self):
        
        plt.rcParams.update({'font.weight': 'normal'})
        plt.rc('lines', linewidth= 3, markersize= 5)
        plt.rcParams.update({'font.size': 14})
        plt.rcParams.update({'figure.facecolor':'w'})
        plt.rcParams.update({'mathtext.default': 'regular'})
    
    def load_numerical_midplane(self, data_size, plot):
        
    
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        b2fstate = self.data['b2fstate']
        
        if data_size == 'full':
            neu_pro = self.data['outputdata']['NeuDen']
            ne_pro = b2fstate['ne'].transpose()
            Te_J = b2fstate['te'].transpose()
            
        elif data_size == 'small':
            data = self.data['ft44']['dab2']
            neu_pro = np.transpose(data[:, :, 0])
            ne_pro = b2fstate['ne'][1:nx+1, 1:ny+1].transpose()
            Te_J = b2fstate['te'][1:nx+1, 1:ny+1].transpose()
            
        
        ev = 1.6021766339999999 * pow(10, -19)
        te_pro = Te_J / ev
        
        leftcut = self.data['b2fgeo']['leftcut'][0]
        rightcut = self.data['b2fgeo']['rightcut'][0]
        
        
        if data_size == 'full':
            psi_coord = self.data['midplane_calc']['psi_solps_mid']
            weight = self.data['midplane_calc']['weight']
            
        elif data_size == 'small':
            psi_coord = self.data['midplane_calc']['psi_solps_mid'][1:ny+1]
            weight = self.data['midplane_calc']['weight'][1:ny+1]
    
        
        weight_B = np.ones(len(weight))- weight
             
        mid_ne_pro = np.multiply(ne_pro[:, 58], weight) + np.multiply(ne_pro[:, 60], weight_B)
        mid_te_pro = np.multiply(te_pro[:, 58], weight) + np.multiply(te_pro[:, 60], weight_B)
        mid_neu_pro = np.multiply(neu_pro[:, 58], weight) + np.multiply(neu_pro[:, 60], weight_B)
        
        
        nete_mid_dic = {'midplane_psiN': psi_coord, 'midplane_ne': mid_ne_pro, 
                        'midplane_te': mid_te_pro}
        
        
        self.data['nete_midpro'] = nete_mid_dic
        
        if plot:
            fig, axs = plt.subplots(1, 2)
            
            anchored_text = AnchoredText('(a){}'.format('$n_e$ [$m^{-3}$]'), loc='upper right')
            axs[0].plot(psi_coord, mid_ne_pro, color = 'black', label= '$n_e$ midplane')
            axs[0].add_artist(anchored_text)
            axs[0].legend(loc='lower left', fontsize=10)
            
            anchored_text = AnchoredText('(b){}'.format('$t_e$ [$m^{-3}$]'), loc='upper right')
            axs[1].plot(psi_coord, mid_ne_pro, color = 'black', label= '$t_e$ midplane')
            axs[1].add_artist(anchored_text)
            axs[1].legend(loc='lower left', fontsize=10)
            
        
        
    
    def profile_combine(self, p_check):
        
        exp_psiN = self.data['ExpDict']['psi_normal']
        exp_ne = self.data['ExpDict']['electron_density(10^20/m^3)']
        exp_te = self.data['ExpDict']['electron_temperature(KeV)']
        
        fit_psiN = self.data['experimental_fit']['psiN']
        fit_ne = self.data['experimental_fit']['ne']
        fit_te = self.data['experimental_fit']['te']
        
        
        exp_sol_psiN = []
        exp_sol_ne = []
        exp_sol_te = []
        for ind, epN in enumerate(exp_psiN):
            if epN >= 1:
                exp_sol_psiN.append(epN)
                exp_sol_ne.append(exp_ne[ind])
                exp_sol_te.append(exp_te[ind])
        
        
        fit_edge_psiN = []
        fit_edge_ne = []
        fit_edge_te = []
        for ind, feN in enumerate(fit_psiN):
            if feN <= 1:
                fit_edge_psiN.append(feN)
                fit_edge_ne.append(fit_ne[ind])
                fit_edge_te.append(fit_te[ind])
        
        
        psiN_std = fit_edge_psiN + exp_sol_psiN
        ne_std = fit_edge_ne + exp_sol_ne
        te_std = fit_edge_te + exp_sol_te
        
        
        
        if p_check:
            
            print('this is the exp sol psiN')
            print(exp_sol_ne)
            print('this is the fit edge psiN')
            print(fit_edge_ne)
            
            k1 = min(exp_sol_psiN)
            k2 = max(fit_edge_psiN)
            weight = (1 - k2)/(k1 - k2)
            
            print(k1)
            print(k2)
            print(weight)
            
            print(psiN_std)
            print(ne_std)
            print(te_std)
        
        
        fig, axs = plt.subplots(1, 2)
        
        anchored_text = AnchoredText('(a){}'.format('$n_e$ [$m^{-3}$]'), loc='upper right')
        axs[0].plot(psiN_std, ne_std, color = 'black', label= '$n_e$ solution')
        axs[0].add_artist(anchored_text)
        axs[0].legend(loc='lower left', fontsize=10)
        
        anchored_text = AnchoredText('(b){}'.format('$t_e$ [$m^{-3}$]'), loc='upper right')
        axs[1].plot(psiN_std, te_std, color = 'black', label= '$t_e$ solution')
        axs[1].add_artist(anchored_text)
        axs[1].legend(loc='lower left', fontsize=10)
        
        
        std_dic = {'psiN_solution': psiN_std, 'ne_solution': ne_std, 
                   'te_solution': te_std}
        
        self.data['nete_ref'] = std_dic
    
    
    def solution_compare(self, plot):
        
        std_dic = self.data['nete_ref']
        nete_mid_dic = self.data['nete_midpro']
        
        psiN_std = std_dic['psiN_solution']
        ne = std_dic['ne_solution']
        te = std_dic['te_solution']
        
        ne_std = [i*pow(10, 20) for i in ne]
        te_std = [i*pow(10, 3) for i in te]
        
        
        mid_psiN = nete_mid_dic['midplane_psiN']
        mid_ne = nete_mid_dic['midplane_ne']
        mid_te = nete_mid_dic['midplane_te']
        
        
        f_ne = interpolate.interp1d(psiN_std, ne_std)
        f_te = interpolate.interp1d(psiN_std, te_std)
        
        
        
        mat_psiN = []
        mat_ne = []
        mat_te = []
        
        for ii, pN in enumerate(mid_psiN):
            
            if pN <= max(psiN_std) and pN >= min(psiN_std):
                
                mat_psiN.append(pN)
                mat_ne.append(mid_ne[ii])
                mat_te.append(mid_te[ii])
        
        
        
        ne_c = f_ne(mat_psiN)
        te_c = f_te(mat_psiN)
        
        if plot:
            
            fig, axs = plt.subplots(1, 2)
            
            anchored_text = AnchoredText('(a){}'.format('$n_e$ [$m^{-3}$]'), loc='upper right')
            axs[0].plot(psiN_std, ne_std, color = 'blue', label= '$n_e$ solution')
            axs[0].plot(mid_psiN, mid_ne, color = 'red', label= '$n_e$ numerical')
            axs[0].add_artist(anchored_text)
            axs[0].legend(loc='lower left', fontsize=10)
            
            anchored_text = AnchoredText('(b){}'.format('$t_e$ [$m^{-3}$]'), loc='upper right')
            axs[1].plot(psiN_std, te_std, color = 'blue', label= '$t_e$ solution')
            axs[1].plot(mid_psiN, mid_te, color = 'red', label= '$t_e$ numerical')
            axs[1].add_artist(anchored_text)
            axs[1].legend(loc='lower left', fontsize=10)
            
            
            fig, axs = plt.subplots(1, 2)
            
            anchored_text = AnchoredText('(a){}'.format('$n_e$ [$m^{-3}$]'), loc='upper right')
            axs[0].scatter(mat_psiN, ne_c, color = 'blue', label= '$n_e$ solution')
            axs[0].scatter(mat_psiN, mat_ne, color = 'red', label= '$n_e$ numerical')
            axs[0].add_artist(anchored_text)
            axs[0].legend(loc='lower left', fontsize=10)
            
            anchored_text = AnchoredText('(b){}'.format('$t_e$ [$m^{-3}$]'), loc='upper right')
            axs[1].scatter(mat_psiN, te_c, color = 'blue', label= '$t_e$ solution')
            axs[1].scatter(mat_psiN, mat_te, color = 'red', label= '$t_e$ numerical')
            axs[1].add_artist(anchored_text)
            axs[1].legend(loc='lower left', fontsize=10)
            
        
        
        
        
        prof_compare_dic = {'x_coord': mat_psiN, 'ne_std': ne_c, 'ne_numerical': mat_ne,
                            'te_std': te_c, 'te_numerical': mat_te }
        
        self.data['profile compare'] = prof_compare_dic
        
    
    
    def error_calculation(self):
        
        pc_dic = self.data['profile compare']
        
        psiN = pc_dic['x_coord']
        ne_std = pc_dic['ne_std']
        ne_num = pc_dic['ne_numerical']
        te_std = pc_dic['te_std']
        te_num = pc_dic['te_numerical']
        
        ne_qsum = 0
        te_qsum = 0
        
        for ii in range(len(psiN)):
            
            ne_qsum = ne_qsum + (ne_std[ii] - ne_num[ii])**2
            te_qsum = te_qsum + (te_std[ii] - te_num[ii])**2
        
        ne_err = np.sqrt(ne_qsum / len(psiN))
        te_err = np.sqrt(te_qsum / len(psiN))
        
        print('ne error')
        print(ne_err/ max(ne_num))
        print('te error')
        print(te_err/ max(te_num))
        
        
    
        
        
        
        
        
        
        