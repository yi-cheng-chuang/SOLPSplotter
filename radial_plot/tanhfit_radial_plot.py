# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 19:16:08 2025

@author: ychuang
"""



from SOLPS_input.header import *
import fitting_method as fm

class radial_tanhfitplot:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    def plot_tanh_fit(self, log_flag):
        
        if self.withshift == True and self.withseries == False:
            
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
            
            psiN = self.data['experimental_fit']['psiN']
            ne = self.data['experimental_fit']['ne']*pow(10, 20)
            te = self.data['experimental_fit']['te']*pow(10, 3)
            
            exp = self.data['ExpDict']
            # psi = exp['psi_normal']
            
            
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
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            
            fig, axs = plt.subplots()
            

            
            aa = 'org'
            
            b2fstate = self.data['b2fstate'][aa]
            
            ne_pro = b2fstate['ne'].transpose()
            Te_J = b2fstate['te'].transpose()
            ev = 1.6021766339999999 * pow(10, -19)
            te_pro = Te_J / ev
            
            data = self.data['ft44'][aa]['dab2']
            neu_pro = np.transpose(data[:, :, 0])
            leftcut = self.data['b2fgeo'][aa]['leftcut'][0]
            rightcut = self.data['b2fgeo'][aa]['rightcut'][0]
            
            weight = self.data['midplane_calc'][aa]['weight']
            weight_B = np.ones(len(weight))- weight
            
            
            mid_ne_pro = np.multiply(ne_pro[:, 58], weight) + np.multiply(ne_pro[:, 60], weight_B)
            mid_te_pro = np.multiply(te_pro[:, 58], weight) + np.multiply(te_pro[:, 60], weight_B)
            
            
        
            psi_coord = self.data['midplane_calc'][aa]['psi_solps_mid']
            
            fit_tanh_dic = fm.tanh_fit(x_coord = psi_coord, ne = mid_ne_pro, te = mid_te_pro)
        
            
            # axs[0].set_yscale('log')
            
            psi_cut = []
            ne_cut = []
            fit_cut = []
            
            for ik, psi_it in enumerate(psi_coord):
                
                if psi_it >= 0.92:
                    psi_cut.append(psi_it)
                    ne_cut.append(mid_ne_pro[ik])
                    fit_cut.append(fit_tanh_dic['tanh_ne_fit'][ik])
            
            if log_flag:
                plt.yscale('log')
            else:
                pass
            
            
            axs.plot(psi_cut, ne_cut, color = 'green', 
                        label= 'SOLPS simulation profile')
            
            axs.plot(psi_cut, fit_cut, color = 'red', 
                        label= 'simulation tanh fit')
            
            
            dn = fit_tanh_dic['popt_ne'][2]
            sym_pt = fit_tanh_dic['popt_ne'][0]
            
            
            axs.axvline(x=dn + sym_pt, ls= '--', color='black',lw=3, 
                        label= 'pedestal width')
            axs.axvline(x=-dn + sym_pt, ls= '--', color='black',lw=3)
            
            simtanh_text = AnchoredText('{}'.format('$n_e$ [$m^{-3}$]'), loc='upper right')
            
            axs.add_artist(simtanh_text)
            axs.set_xlabel('$\psi_N$')
            
            
            axs.legend(loc= 'lower left', fontsize=10)






