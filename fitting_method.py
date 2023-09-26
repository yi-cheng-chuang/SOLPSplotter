# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 13:28:56 2023

@author: user
"""

import numpy as np
from scipy.optimize import curve_fit
import load_mast_expdata_method as lmem
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from statistics import stdev

def tanh(r,r0,h,d,b,m):
    return b+(h/2)*(np.tanh((r0-r)/d)+1) + m*(r0-r-d)*np.heaviside(r0-r-d, 1)

def expfit(x,A,l):  #Removed vertical displacement variable B; seemed to cause 'overfitting'
    return A*np.exp(l*x)


def tanh_fit(x_choice, x_coord, ne, te):
    # mast_dat_dict = lmem.read_mastfile(mastloc)
    # psi = mast_dat_dict['psi_normal']
    # ne = mast_dat_dict['electron_density(10^20/m^3)']
    # te = mast_dat_dict['electron_temperature(KeV)']
    
    Ne = ne*pow(10, -19)
    Te = te*pow(10, -3)
    
    if x_choice == 'psiN':
        p0 = [1, 0.4, 0.009, 0.05, 0.5]
        p1 = [1, 0.3, 0.1, 0.001, 0.5]
        popt_ne, pcov_ne = curve_fit(tanh, x_coord, Ne, p0)
        # print(popt_ne)
        popt_te, pcov_te = curve_fit(tanh, x_coord, Te, p1)
        # print(popt_te) 
        tanh_ne_fit = tanh(x_coord, popt_ne[0], popt_ne[1], 
                               popt_ne[2], popt_ne[3], popt_ne[4])*pow(10, 19)
        tanh_te_fit = tanh(x_coord, popt_te[0], popt_te[1], 
                               popt_te[2], popt_te[3], popt_te[4])*pow(10, 3)
    
    
    elif x_choice == 'RRsep':
        p0 = [0.4, 0.005, 0.05, 0.5]
        p1 = [2.2, 0.005, 0.05, 3.5]
        popt_ne, pcov_ne = curve_fit(tanh_dsa, x_coord, Ne, p0)
        print(popt_ne)
        popt_te, pcov_te = curve_fit(tanh_dsa, x_coord, Te, p1)
        # print(popt_te) 
        tanh_ne_fit = tanh_dsa(x_coord, popt_ne[0], popt_ne[1], 
                               popt_ne[2], popt_ne[3])*max(ne)
        tanh_te_fit = tanh_dsa(x_coord, popt_te[0], popt_te[1], 
                               popt_te[2], popt_te[3])*max(te)
    
    
    
    fit_tanh_dic = {'tanh_ne_fit': tanh_ne_fit, 'tanh_te_fit': tanh_te_fit, 
               'popt_ne': popt_ne, 'popt_te': popt_te}
    
    
    
    
    return fit_tanh_dic


def exp_fit_sp(x_choice, x_coord, neuden):
    
    
    NeuDen = neuden/ neuden[0]
    
    if x_choice == 'psiN':
        # pn = [1, 1.015, 200.5]
        pn = [1.015, 200.5]
        x_sh = x_coord - x_coord[0]
        popt_an, pcov_an = curve_fit(expfit, x_sh, NeuDen, pn)
        # print(popt_an)
        
        exp_an_fit = expfit(x_sh, popt_an[0], popt_an[1])*neuden[0]
        
    elif x_choice == 'RRsep':
        pn = [1.015, 200.5]
        popt_an, pcov_an = curve_fit(expfit, x_coord, NeuDen, pn)
        print(popt_an)
        
        exp_an_fit = expfit(x_coord, popt_an[0], popt_an[1])*max(neuden)
        
    
    fit_exp_dic = {'exp_an_fit': exp_an_fit, 'popt_an': popt_an}
    
    return fit_exp_dic



def exp_fit(x_choice, x_coord, neuden):
    NeuDen = neuden/ neuden[0]
    
    if x_choice == 'psiN':
        # pn = [1, 1.015, 200.5]
        pn = [1.015, 200.5]
        x_sh = x_coord - x_coord[0]
        popt_an, pcov_an = curve_fit(expfit, x_sh, NeuDen, pn)
        # print(popt_an)
        
        exp_an_fit = expfit(x_sh, popt_an[0], popt_an[1])*neuden[0]
    
    fit_exp_dic = {'exp_an_fit': exp_an_fit, 'popt_an': popt_an}
    
    # plt.figure(figsize=(7,7))
    # plt.plot(dsa, neuden,'o-', color = 'green', label= 'solps neutral density')
    # plt.plot(dsa, exp_an_fit, color='r',lw= 5, label= 'exponential fit')
    # plt.xlabel('psiN')
    # plt.ylabel('Neutral Atom (D) Density $(m^{-3})$')
    # plt.title('Neutral density with fits')
    # plt.legend()
    
    return fit_exp_dic



def flux_expansion_fit(psi, dsa, flux_expansion_jxa):
    flux_expansion = np.polyfit( dsa, psi,1,cov=True)
    fluxpsn = flux_expansion[0][0]/flux_expansion_jxa[0][0]
    slope = flux_expansion[0][0]
    
    return slope
    
#     fluxpsnparam[PolPos] = np.polyfit(RR_exp,Psin_avg[PolPos][RR_i],1,cov=True)
#     fluxpsn[PolPos] = fluxpsnparam[PolPos][0][0]/fluxpsnparam[JXA][0][0]
#     fluxpsn_err[PolPos] = np.sqrt(np.diag(fluxpsnparam[PolPos][1]))[0]
#     efold_adj[PolPos] = fluxpsn[PolPos]*efold[PolPos]
#     efold_adj_err[PolPos] = efold_adj[PolPos]*np.sqrt((efold_err[PolPos]/efold[PolPos])**2 + (fluxpsn_err[PolPos]/fluxpsn[PolPos])**2)
            
# print('Exponential fit from r-r_sep={:.1e}m to r-r_sep={:.1e}m'.format(RR_exp[0],RR_exp[-1]))
# print('A0={:.3e}, lambda={:.1f}'.format(*eparam[PolPos]))
# print('Poloidal Slice {:0.0f}: Raw e-folding length={:.1f}mm, std error={:.1f}'.format(PolPos,efold[PolPos],efold_err[PolPos]))
# print('Slope={:.3f}, Flux Expansion={:.1f}, Adjusted e-folding length={:.1f}, std error={:.1f}'.format(fluxpsnparam[PolPos][0][0],fluxpsn[PolPos],efold_adj[PolPos],efold_adj_err[PolPos]))

def line(x, m, b):
    return m*x+b


"fit dsa with psi"
def dsa_psi_fit(dsa, psi, dsa_cut):
    
    dsa_psi_fitcoe = np.polyfit(dsa, psi, 1 , cov=True)
    # print(dsa_psi_fitcoe[1])
    dsa_psi_fitpoly = np.poly1d(dsa_psi_fitcoe[0])
    
    
    dsa_psi_fit = dsa_psi_fitpoly(dsa_cut)
   
    fit_dp_dic = {'dsa_psi_fit': dsa_psi_fit, 'dsa_psi_fitcoe': dsa_psi_fitcoe[0]}
    
    # plt.figure(figsize=(7,7))
    # plt.plot(dsa, psi,'o-', color = 'green', label= 'solps neutral density')
    # plt.plot(dsa_cut, psi_cut, color='r',lw= 5, label= 'exponential fit')
    # plt.xlabel('RR_sep')
    # plt.ylabel('psiN')
    # plt.title('dsa_psi_fit')
    # plt.legend()
    
    return fit_dp_dic



def dsa_psi_fit_sp(dsa, psi, dsa_cut, psi_cut):
    
    dsa_psi_fitcoe = np.polyfit( dsa_cut, psi_cut, 1 , cov=True)
    # print(dsa_psi_fitcoe[1])
    dsa_psi_fitpoly = np.poly1d(dsa_psi_fitcoe[0])
    
    
    dsa_psi_fit = dsa_psi_fitpoly(dsa_cut)
   
    fit_dp_dic = {'dsa_psi_fit': dsa_psi_fit, 'dsa_psi_fitcoe': dsa_psi_fitcoe[0]}
    
    # plt.figure(figsize=(7,7))
    # plt.plot(dsa, psi,'o-', color = 'green', label= 'solps neutral density')
    # plt.plot(dsa_cut, psi_cut, color='r',lw= 5, label= 'exponential fit')
    # plt.xlabel('RR_sep')
    # plt.ylabel('psiN')
    # plt.title('dsa_psi_fit')
    # plt.legend()
    
    return fit_dp_dic



def Opacity_calculator(x_choice, x_coord, ne, te, neuden):
    if x_choice == 'psiN':
        fit_tanh_dic = tanh_fit(x_choice, x_coord, ne, te)
        tanh_ne_fit = fit_tanh_dic['tanh_ne_fit']
        # dn = fit_tanh_dic['popt_ne'][2]
        dn = fit_tanh_dic['popt_ne'][2]
        tanh_te_fit = fit_tanh_dic['tanh_te_fit']
        # dtn = fit_tanh_dic['popt_te'][2]
        dtn = fit_tanh_dic['popt_te'][2]   
        ne_ped = (fit_tanh_dic['popt_ne'][1] + fit_tanh_dic['popt_ne'][3])*pow(10, 19)
        sym_pt_te = fit_tanh_dic['popt_te'][0]
        
        "fitting choice 1: in the width, symmetry point +- width/2"
        xcoord_exp = []
        an_cut = []
        sym_pt = fit_tanh_dic['popt_ne'][0]
        for j in range(len(x_coord)):
            if x_coord[j] <= sym_pt + dn and x_coord[j] >= sym_pt - dn:
                xcoord_exp.append(x_coord[j])
                an_cut.append(neuden[j])
        
        "plot to check the tanh fit result"
        # plt.figure(figsize=(7,7))
        # plt.plot(x_coord, ne,'o-', color = 'b', label= 'solps electron density')
        # plt.plot(x_coord, tanh_ne_fit, color='r',lw= 3, label= 'tanh fit')
        # plt.xlabel('Radial coordinate: $R- R_{sep}$')
        # plt.ylabel('Electron Density $n_e\;(m^{-3})$')
        # plt.title('Electron density with fits')
        # plt.legend()
        
        "print to check the exp fit result"
        # print(dn)
        # print(sym_pt)
        # print(xcoord_exp)
        # print(an_cut)
        
        "fitting choice 2: find the point before the kink"
        xcoord_k = []
        an_k_cut = []
        for af in range(len(x_coord)):
            if x_coord[af] <= sym_pt + dn:
                xcoord_k.append(x_coord[af])
                an_k_cut.append(neuden[af])
        
        pen_try = []
        
        "create try and error list"
        n_list = []
        # nn_i = len(xcoord_exp)
        nn_i = 4
        # print(len(xcoord_k))
        while nn_i < len(xcoord_k) + 1:
            n_list.append(nn_i)
            nn_i = nn_i + 1
        # print(n_list)
        
        an_tmp_dic = {}
        x_tmp_dic = {}
        # pen_try = []
        
        for ak in n_list:
            x_tmp = []
            an_tmp = []
            # print(ak)
            for ah in range(ak):
                x_tmp.append(xcoord_k[-ah - 1])
                an_tmp.append(an_k_cut[-ah - 1])
            x_tmp = np.asarray(x_tmp)
            an_tmp = np.asarray(an_tmp)
            x_tmp_dic[str(ak)] = x_tmp
            an_tmp_dic[str(ak)] = an_tmp   
            # print(ak)
            # print(x_tmp)
            # print(an_tmp)
            fitexp_tmp = exp_fit_sp(x_choice, x_tmp, an_tmp)
            pen_try.append(1/fitexp_tmp['popt_an'][1])
            
        diff_pen = np.diff(pen_try)
        
        
        "print to check fit method 2 result"
        # print(n_list)
        # print(n_list[-1])
        # print(x_tmp_dic[str(n_list[-1])])
        # print(pen_try)
        # print(diff_val)
        
        "plot the method 2 result"
        # plt.figure(figsize=(7,7))
        # plt.errorbar(x_tmp_dic.keys(), pen_try, yerr= stdev(pen_try), fmt = 'o-', color = 'b', label= 'efold_length')
        # plt.xlabel('number of fitting points')
        # # plt.ylabel('efold_length')
        # plt.title('efold length for fitting method 2')
        # plt.legend()
        
           
        dis_pen = []
        stable_pen = []
        for i_diff in diff_pen:
            dis_pen.append(abs(i_diff))
            if abs(i_diff) <= 0.001:
                stable_pen.append(abs(i_diff))
        
        # print(dis_pen)
        
        p = list(dis_pen).index(stable_pen[0])
        # print(diff_pen)
        ind = n_list[p]
        # if p + 1 < len(n_list):
        #     # print(len(n_list))
        #     # print(p)
        #     ind = n_list[p + 1]
        #     # print(ind)
        # elif p + 1 >= len(n_list):
        #     # print(len(n_list))
        #     # print(p)
        #     ind = n_list[p]
        #     # print(ind)
        x_m2 = np.asarray(x_tmp_dic[str(ind)])
        an_m2 = np.asarray(an_tmp_dic[str(ind)])
        
        # print(x_m2)
        
        
        fit_m2_dic = exp_fit_sp(x_choice, x_m2, an_m2)
        
        exp_fit_m2 = fit_m2_dic['exp_an_fit']
        efold_m2 = 1/fit_m2_dic['popt_an'][1]
        std_m2 = stdev(pen_try)/np.sqrt(len(pen_try))
        
        opq_m2 = 2*dn/efold_m2
        
        
        x_m3 = []
        an_m3 = []
        
        # i_nd = 1
        # for i in range(len(x_tmp_dic[str(n_list[-1])]) - 1):
        #     if i_nd >= ind - 1:
        #         x_m3.append(x_tmp_dic[str(n_list[-1])][i_nd])
        #         an_m3.append(an_tmp_dic[str(n_list[-1])][i_nd])
        #     i_nd = i_nd + 1
        
        # print(x_tmp_dic[str(n_list[-1])])
        
        diff_val = np.diff(np.log10(an_tmp_dic[str(n_list[-1])]))
        diff_val = np.insert(diff_val, 0, 0)
              
        dis_val = []
        for i_val in diff_val:
            dis_val.append(abs(i_val))
        
        "plot to check the result of fitting method 3"
        # plt.figure(figsize=(7,7))
        # plt.plot(x_tmp_dic[str(n_list[-1])], dis_val,'o-', color = 'b', label= 'neutral density')
        # plt.xlabel('magnetic flux coordinate: psiN')
        # # plt.ylabel('neutral density difference')
        # plt.title('log neutral density difference for method 3')
        # plt.legend()
        
                      
        qth = list(dis_val).index(max(dis_val))
        # qin = n_list[qth]
        # print(qin)
        
        
        if qth >= 4:
            x_m3 = np.asarray(x_tmp_dic[str(qth)])
            an_m3 = np.asarray(an_tmp_dic[str(qth)])
        elif qth < 4:
            x_m3 = np.asarray(x_tmp_dic[str(4)])
            an_m3 = np.asarray(an_tmp_dic[str(4)])
        
        
        fit_m3_dic = exp_fit_sp(x_choice, x_m3, an_m3)

        
        exp_fit_m3 = fit_m3_dic['exp_an_fit']
        efold_m3 = 1/fit_m3_dic['popt_an'][1]
        
        opq_m3 = 2*dn/efold_m3
        
        
        xcoord_exp = np.asarray(xcoord_exp)
        an_cut = np.asarray(an_cut)
        
        fit_exp_dic = exp_fit_sp(x_choice, xcoord_exp, an_cut)

        
        exp_an_fit = fit_exp_dic['exp_an_fit']
        efold = 1/fit_exp_dic['popt_an'][1]
        
        opq = 2*dn/efold
        
        "fit_dsa_psi"
        # xcoord_cut_dn = []
        # for j in range(len(xcoord)):
        #     if x_coord[j] <= fit_tanh_dic['popt_ne'][0] + dn and x_coord[j] >= fit_tanh_dic['popt_ne'][0] - dn:
        #         xcoord_cut_dn.append(xcoord[j])
        
        # psi_cut_dn = []
        # for j in range(len(psi)): 
        #     if psi[j] <= 2*dn + 1 and psi[j] >= -2*dn + 1:
        #         psi_cut_dn.append(psi[j])
        
        # # print(len(dsa_cut_dn))
        # # print(len(psi_cut_dn))
        # ir = min(len(dsa_cut_dn), len(psi_cut_dn))
        
        # psi_cut = np.zeros(ir)
        # dsa_cut = np.zeros(ir)
        # for i in range(ir):
        #     psi_cut[i] = psi_cut_dn[i]
        #     dsa_cut[i] = dsa_cut_dn[i]
        
        
        # dp_dic = dsa_psi_fit(dsa = dsa, psi = psi, dsa_cut = dsa_cut)
        # slope = dp_dic['dsa_psi_fitcoe'][0]
        # dp_fit = dp_dic['dsa_psi_fit']
        
        result_dic = {'tanh_ne_fit': tanh_ne_fit, 'exp_fit': exp_an_fit,
                      'electron_pedestal_density': ne_ped, 'x_coord_cut': xcoord_exp,
                      'tanh_te_fit': tanh_te_fit, 'pedestal_width': dn, 
                      'temperature_pedestal_width': dtn,
                      'efold_length': efold, 'dimensionless_opaqueness': opq,
                      'ne_symmetry_point': sym_pt, 'te_symmetry_point': sym_pt_te,
                      'penetration_length_change': pen_try, 'an_tmp_dic': an_tmp_dic,
                      'diff_pen': diff_pen,
                      
                      'exp_fit_m2': exp_fit_m2, 
                      'efold_length_method2': efold_m2, 
                      'dimensionless_opaqueness_method2': opq_m2,
                      'x_m2': x_m2, 'std_m2': std_m2,
                      
                      'exp_fit_m3': exp_fit_m3, 
                      'efold_length_method3': efold_m3, 
                      'dimensionless_opaqueness_method3': opq_m3, 
                      'x_m3': x_m3}
        
        
        return result_dic
        
    
    elif x_choice == 'RRsep':
        fit_tanh_dic = tanh_fit(x_choice, x_coord, ne, te)
        tanh_ne_fit = fit_tanh_dic['tanh_ne_fit']
        # dn = fit_tanh_dic['popt_ne'][2]
        dn = fit_tanh_dic['popt_ne'][1]
        tanh_te_fit = fit_tanh_dic['tanh_te_fit']
        # dtn = fit_tanh_dic['popt_te'][2]
        dtn = fit_tanh_dic['popt_te'][1]   
        ne_ped = (fit_tanh_dic['popt_ne'][0] + fit_tanh_dic['popt_ne'][2])*max(ne)
        
        "fitting choice 1: in the width, symmetry point +- width/2"
        xcoord_exp = []
        an_cut = []
        for j in range(len(x_coord)):
            if x_coord[j] <= fit_tanh_dic['popt_ne'][0] + dn and x_coord[j] >= fit_tanh_dic['popt_ne'][0] - dn:
                xcoord_exp.append(x_coord[j])
                an_cut.append(neuden[j])
        
        xcoord_exp = np.asarray(xcoord_exp)
        an_cut = np.asarray(an_cut)
        
        fit_exp_dic = exp_fit(xcoord_exp, an_cut)

    
    
    


#--------------------spare-zone--------------------------------------------
