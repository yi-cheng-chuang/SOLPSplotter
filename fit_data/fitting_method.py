# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 13:28:56 2023

@author: user
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from lmfit import Model



class fit_method_collection:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    def tanh(self, r,r0,h,d,b,m):
        return b+(h/2)*(np.tanh((r0-r)/d)+1) + m*(r0-r-d)*np.heaviside(r0-r-d, 1)

    def expfit(self, x,A,l):  #Removed vertical displacement variable B; seemed to cause 'overfitting'
        return A*np.exp(l*x)
    
    def tanh0(self, r,r0,h,d,b):
        return b+(h/2)*(np.tanh((r0-r)/d)+1)


    def tanh_fit(self, x_coord, ne, te):
        
        # Ne = ne*pow(10, -19)
        # Te = te*pow(10, -3)
        
        
        Ne = [x *pow(10, -19) for x in ne]
        Te = [x*pow(10, -3) for x in te]
          
        
        
        p0 = [1, 0.4, 0.009, 0.05, 0.5]
        p1 = [1, 0.3, 0.1, 0.001, 0.5]
        popt_ne, pcov_ne = curve_fit(self.tanh, x_coord, Ne, p0)
        # print(popt_ne)
        popt_te, pcov_te = curve_fit(self.tanh, x_coord, Te, p1)
        # print(popt_te) 
        tanh_ne_fit = self.tanh(x_coord, popt_ne[0], popt_ne[1], 
                               popt_ne[2], popt_ne[3], popt_ne[4])*pow(10, 19)
        tanh_te_fit = self.tanh(x_coord, popt_te[0], popt_te[1], 
                               popt_te[2], popt_te[3], popt_te[4])*pow(10, 3)
        
          
        fit_tanh_dic = {'tanh_ne_fit': tanh_ne_fit, 'tanh_te_fit': tanh_te_fit, 
                   'popt_ne': popt_ne, 'popt_te': popt_te}
           
        return fit_tanh_dic


    def exp_fit(self, x_coord, neuden):
        
        
        NeuDen = neuden/ neuden[0]
        

        pn = [1.015, 200.5]
        x_sh = x_coord - x_coord[0]
        popt_an, pcov_an = curve_fit(self.expfit, x_sh, NeuDen, pn)
        # print(popt_an)
        
        exp_an_fit = self.expfit(x_sh, popt_an[0], popt_an[1])*neuden[0]
        
        # plt.figure(figsize=(7,7))
        # plt.yscale('log')
        # plt.plot(x_coord, neuden,'o-', color = 'green', label= 'solps neutral density')
        # plt.plot(x_coord, exp_an_fit, color='r',lw= 5, label= 'exponential fit')
        # plt.xlabel('psiN')
        # plt.ylabel('Neutral Atom (D) Density $(m^{-3})$')
        # plt.title('Neutral density with fits')
        # plt.legend()
            
        
        fit_exp_dic = {'exp_an_fit': exp_an_fit, 'popt_an': popt_an}
        
        return fit_exp_dic



    def linear_fit(self, x_coord, neuden):
        
        NeuDen = np.log(neuden)
        
        
        
        x_sh = x_coord - x_coord[0]
        # print(x_coord[0])
        
        
        # x_sh = x_coord - 1
        
        ln_exp_fitcoe = np.polyfit(x_sh, NeuDen, 1 , cov=True)
        
        # print('exp coe')
        # print(ln_exp_fitcoe[0])
        ln_exp_fitpoly = np.poly1d(ln_exp_fitcoe[0])
        
        # print('exp poly')
        # print(ln_exp_fitpoly)
        exp_ln_fit = np.exp(ln_exp_fitpoly(x_sh))
        
        # plt.figure(figsize=(7,7))
        # plt.yscale('log')
        # plt.plot(x_coord, neuden,'o-', color = 'green', label= 'solps neutral density')
        # plt.plot(x_coord, exp_line_fit, color='r',lw= 5, label= 'exponential fit')
        # plt.xlabel('psiN')
        # plt.ylabel('Neutral Atom (D) Density $(m^{-3})$')
        # plt.title('Neutral density with fits')
        # plt.legend()
            
        fit_exp_dic = {'exp_line_fit': exp_ln_fit, 'popt_an': ln_exp_fitcoe[0]}
        
        return fit_exp_dic 



    "fit dsa with psi"
    def dsa_psi_fit(self, dsa, psi):
        dsa_psi_fitcoe = np.polyfit(psi, dsa, 1 , cov=True)
        # print(dsa_psi_fitcoe[1])
        dsa_psi_fitpoly = np.poly1d(dsa_psi_fitcoe[0])
        
        
        psi_dsa_fit = dsa_psi_fitpoly(psi)
       
        fit_dp_dic = {'dsa_psi_fit': psi_dsa_fit, 'dsa_psi_fitcoe': dsa_psi_fitcoe[0]}
        
        # plt.figure(figsize=(7,7))
        # plt.plot(psi, dsa,'o-', color = 'green', label= 'psi_dsa')
        # plt.plot(psi, psi_dsa_fit, color='r',lw= 5, label= 'psi_dsa_fit')
        # plt.xlabel('psiN')
        # plt.ylabel('R-Rsep')
        # plt.title('psi_dsa_fit')
        # plt.legend()
        
        # print(dsa_psi_fitcoe[0])

        return fit_dp_dic


    def flux_expand_fit(self, RRsep, arclength):
        flux_fitcoe = np.polyfit(RRsep, arclength, 1 , cov=True)
        # print(dsa_psi_fitcoe[1])
        flux_fitpoly = np.poly1d(flux_fitcoe[0])
              
        flux_fit = flux_fitpoly(RRsep)     
        flux_fit_dic = {'flux_fit': flux_fit, 'flux_fitcoe': flux_fitcoe[0]}
        
        # plt.figure(figsize=(7,7))
        # plt.plot(RRsep, arclength,'o-', color = 'green', label= 'psi_dsa')
        # plt.xlabel('R-Rsep')
        # plt.ylabel('arclength')
        # plt.title('flux_expansion_fit')
        # plt.legend()
        
        # print(flux_fitcoe[0])
            
        return flux_fit_dic


    def Opacity_calculator(self, x_coord, ne, te, neuden, check_ne, minor_rad):

            fit_tanh_dic = self.tanh_fit(x_coord, ne, te)
            tanh_ne_fit = fit_tanh_dic['tanh_ne_fit']
            dn = fit_tanh_dic['popt_ne'][2]
            tanh_te_fit = fit_tanh_dic['tanh_te_fit']
            dtn = fit_tanh_dic['popt_te'][2]  
            ne_ped = (fit_tanh_dic['popt_ne'][1] + fit_tanh_dic['popt_ne'][3])*pow(10, 19)
            te_ped = (fit_tanh_dic['popt_te'][1] + fit_tanh_dic['popt_te'][3])*pow(10, 3)
            sym_pt_te = fit_tanh_dic['popt_te'][0]
            sym_pt_ne = fit_tanh_dic['popt_ne'][0]
            # print(sym_pt_te + dn)
            
            popt_ne = fit_tanh_dic['popt_ne']
            popt_te = fit_tanh_dic['popt_te']
            
            ne_sep = self.tanh(1, popt_ne[0], popt_ne[1], 
                                   popt_ne[2], popt_ne[3], popt_ne[4])*pow(10, 19)
            
            te_sep = self.tanh(1, popt_te[0], popt_te[1], 
                                   popt_te[2], popt_te[3], popt_te[4])*pow(10, 3)
            
            # print(ne_sep)
            avg_ne = 0.5*(ne_ped + ne_sep)
            approxi_opq = avg_ne*minor_rad
            
            
            
            "fitting choice 1: in the width, symmetry point +- width/2"
            xcoord_exp = []
            an_cut = []
            index_cut = []
            sym_pt = fit_tanh_dic['popt_ne'][0]
            x_coord_rev = list(reversed(x_coord))
            # print(type(x_coord_rev))
            neuden_rev = list(reversed(neuden))
            # print(type(x_coord_rev))
            mod_dn = 0.5*np.log(2 + np.sqrt(3))*dn
            for j in range(len(x_coord_rev)):
                if x_coord_rev[j] <= sym_pt_ne + dn and x_coord_rev[j] >= sym_pt_ne - dn:
                    xcoord_exp.append(x_coord_rev[j])
                    an_cut.append(neuden_rev[j])
                    index_cut.append(j)
            
            # print(sym_pt, dn)
            # print(xcoord_exp)
            
            "plot to check the tanh fit result"
            
            if check_ne:
                
                plt.figure()
                plt.plot(x_coord, ne,'o-', color = 'b', label= 'solps electron density')
                plt.plot(x_coord, tanh_ne_fit, color='r',lw= 3, label= 'tanh fit')
                plt.scatter(1, ne_sep, color = 'g', label = 'separatrix density')
                plt.xlabel('Radial coordinate: $R- R_{sep}$')
                plt.ylabel('Electron Density $n_e\;(m^{-3})$')
                plt.title('Electron density with fits')
                plt.legend()
                
            

            
            "print to check the exp fit result"
            # print(dn)
            # print(sym_pt)
            # print(xcoord_exp)
            # print(an_cut)
            
            xcoord_exp = np.asarray(xcoord_exp)
            an_cut = np.asarray(an_cut)
            
            fit_exp_dic = self.linear_fit(xcoord_exp, an_cut)

            
            exp_an_fit = fit_exp_dic['exp_line_fit']
            efold = 1/fit_exp_dic['popt_an'][0]
            n_sep_fit = np.exp(fit_exp_dic['popt_an'][1])
            
            opq = 2*dn/efold
            
            
            result_dic = {'tanh_ne_fit': tanh_ne_fit, 'exp_fit': exp_an_fit,
                          'electron_pedestal_density': ne_ped,
                          'electron_pedestal_temperature': te_ped,
                          'x_coord_cut': xcoord_exp,
                          'tanh_te_fit': tanh_te_fit, 'pedestal_width': dn, 
                          'temperature_pedestal_width': dtn,
                          'efold_length': efold, 'dimensionless_opaqueness': opq,
                          'ne_symmetry_point': sym_pt, 'te_symmetry_point': sym_pt_te,
                          'n_sep_fit': n_sep_fit, 'sep_index': index_cut,
                          'opaqueness_approximation': approxi_opq,
                          'electron_density_separatrix': ne_sep,
                          'electron_temperature_separatrix': te_sep
                          }
            
            
            return result_dic
    
    
    def lmfit_tanh(self, TS_data, psi, plot):
        
        ne_sort = TS_data['ne']
        te_sort = TS_data['te']
        
        
        Ne = ne_sort
        Te = te_sort
        

        # Create the model
        model = Model(self.tanh)

        # Provide initial parameter guesses
        params = model.make_params(r0 = 1,h = 0.6,d = 0.01,b = 0.01,m = 3/14)
        # params['b'].set(min=0.1, max=10)

        # Fit to the data
        result_Te = model.fit(Te, params, r= psi)
        result_Ne = model.fit(Ne, params, r= psi)

        "Show results"
        print(result_Ne.fit_report())

        #Plot
        
        if plot:
            
            plt.figure()
            plt.plot(psi, Te, 'bo', label='data')
            plt.plot(psi, result_Te.best_fit, 'r-', label='fit')
            plt.title('check out Te lmfit')
            plt.legend()


            plt.figure()
            plt.plot(psi, Ne, 'bo', label='data')
            plt.plot(psi, result_Ne.best_fit, 'r-', label='fit')
            plt.title('check out Ne lmfit')
            plt.legend()
            
            
            
            n_tot = 200
            psi_pro = np.linspace(psi.min(), psi.max(), num= n_tot, dtype= float)
            
            
            
            # Evaluate the fitted model at x_new
            Te_new = result_Te.eval(r = psi_pro)
            Ne_new = result_Ne.eval(r = psi_pro)
            
            # Plot the new prediction
            plt.figure()
            plt.plot(psi, Te, 'bo', label='data')
            plt.plot(psi_pro, Te_new, 'g--', label='fit (extrapolated)')
            plt.legend()
            
            
            plt.figure()
            plt.plot(psi, Ne, 'bo', label='data')
            plt.plot(psi_pro, Ne_new, 'g--', label='fit (extrapolated)')
            plt.legend()
            plt.show()
            
            
            plt.show()





        


#--------------------spare-zone--------------------------------------------

# elif x_choice == 'RRsep':
#     fit_tanh_dic = tanh_fit(x_choice, x_coord, ne, te)
#     tanh_ne_fit = fit_tanh_dic['tanh_ne_fit']
#     # dn = fit_tanh_dic['popt_ne'][2]
#     dn = fit_tanh_dic['popt_ne'][1]
#     tanh_te_fit = fit_tanh_dic['tanh_te_fit']
#     # dtn = fit_tanh_dic['popt_te'][2]
#     dtn = fit_tanh_dic['popt_te'][1]   
#     ne_ped = (fit_tanh_dic['popt_ne'][0] + fit_tanh_dic['popt_ne'][2])*max(ne)
    
#     "fitting choice 1: in the width, symmetry point +- width/2"
#     xcoord_exp = []
#     an_cut = []
#     for j in range(len(x_coord)):
#         if x_coord[j] <= fit_tanh_dic['popt_ne'][0] + dn and x_coord[j] >= fit_tanh_dic['popt_ne'][0] - dn:
#             xcoord_exp.append(x_coord[j])
#             an_cut.append(neuden[j])
    
#     xcoord_exp = np.asarray(xcoord_exp)
#     an_cut = np.asarray(an_cut)
    
#     fit_exp_dic = exp_fit(xcoord_exp, an_cut)
