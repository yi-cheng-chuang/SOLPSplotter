# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 19:13:47 2025

@author: ychuang
"""



from SOLPS_input.header import *
import fitting_method as fm




class radial_plot_targets:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
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
            
            
            fig, axs = plt.subplots(4, 1)
            
            
            
            for ii, side in enumerate(div_side_list):
                
                # axs.figure(figsize=(7,7))
                
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
                        axs.yscale('log')
                    else:
                        pass
                    
                        
                    axs[ii].plot(psiN, plot_data, '-', color = color_dic[aa], 
                                 label = '{}'.format(A_dic[aa]))
                    
                axs[3].legend()
                plt.xlabel('psiN')
                axs[0].set_title('{}'.format(quant))
                # plt.show()
            
            plt.subplots_adjust(hspace=.0)
                
        
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
                
                plt.figure()
                
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
                  
                        
    def divertor_te(self, sep_plot):  
        
        b2fstate = self.data['b2fstate']

        if self.withshift == True and self.withseries == False:
            
            if sep_plot:
                
                plt.figure()
                
                color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                             'dot7': 'blue', 'one': 'purple'}
                A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                          'dot7': '2.8', 'one': '3.4'}
        
                for aa in self.data['dircomp']['multi_shift']:
                    
                    
                    
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
              
                
                plt.figure()
        
                for aa in self.data['dircomp']['multi_shift']:
                    

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
                
            else:
                
                fig, axs = plt.subplots(1, 2)
                anchored_text_in = AnchoredText('(a){}'.format('inner target'), loc=2)
                anchored_text_out = AnchoredText('(b){}'.format('outer target'), loc=2)
                
                for aa in self.data['dircomp']['multi_shift']:

                    

                    rcood = self.data['psi']['psival'][aa][:, 1]
                    Te_J = b2fstate[aa]['te'].transpose()
                    ev = 1.6021766339999999 * pow(10, -19)
                    te_pro = Te_J / ev
                    te_in = te_pro[:, 0]
                    te_out = te_pro[:, -1]
                    axs[0].plot(rcood, te_in, '-', color = color_dic[aa])
                    axs[1].plot(rcood, te_out, '-', color = color_dic[aa])
                
                
                axs[0].add_artist(anchored_text_in)
                axs[0].set_xlabel('$\psi_N$')
                axs[1].add_artist(anchored_text_out)
                axs[1].set_xlabel('$\psi_N$')
                
         
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
                
                anchored_text = AnchoredText('(a){}'.format('Electron density'), loc=2)
                axs[0].set_yscale('log')
                axs[0].errorbar(psiN, mean_core_ne, yerr= std_core_ne, fmt = '-', color = 'g', label= 'ne_solps')
                axs[0].errorbar(psi, exp_ne, yerr= ne_er, fmt = 'o', color = 'b', label= 'ne_exp')
                # plt.plot(psiN, ne, 'o', color = 'r', label= 'ne_exp_fit')
                axs[0].set_xlabel('Normalized flux coordinate $\psi_N$')
                axs[0].set_title('(a)Electron density')
                # axs[0].add_artist(anchored_text)
                
                
                anchored_text2 = AnchoredText('(b){}'.format('Electron temperature'), loc=2)
                axs[1].set_yscale('log')
                axs[1].errorbar(psiN, mean_core_te, yerr= std_core_te, fmt = '-', color = 'g', label= 'te_solps')
                axs[1].errorbar(psi, exp_te, yerr= te_er, fmt = 'o', color = 'b', label= 'te_exp')
                # plt.plot(psiN, ne, 'o', color = 'r', label= 'ne_exp_fit')
                axs[1].set_xlabel('Normalized flux coordinate $\psi_N$')
                axs[1].set_title('(b)Electron temperature')
                # axs[1].add_artist(anchored_text2)
                
            
            
            
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
        
        
        elif self.withshift == True and self.withseries == False:
            
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
            
            fig, axs = plt.subplots(2, 1, figsize= (7, 7))
            
            anchored_text = AnchoredText('(a){}'.format('Electron density'), loc=3)
            axs[0].errorbar(psi, exp_ne, yerr= ne_er, fmt = 'o', color = 'purple', label= 'ne_exp')
            # plt.plot(psiN, ne, 'o', color = 'r', label= 'ne_exp_fit')
            # axs[0].set_xlabel('Normalized flux coordinate $\psi_N$')
            # axs[0].set_title('(a)Electron density')
            axs[0].add_artist(anchored_text)
            axs[0].legend(loc='center left')
            
            
            anchored_text2 = AnchoredText('(b){}'.format('Electron temperature'), loc=3)
            axs[1].errorbar(psi, exp_te, yerr= te_er, fmt = 'o', color = 'purple', label= 'te_exp')
            # plt.plot(psiN, ne, 'o', color = 'r', label= 'ne_exp_fit')
            axs[1].set_xlabel('Normalized flux coordinate $\psi_N$')
            # axs[1].set_title('(b)Electron temperature')
            axs[1].add_artist(anchored_text2)
            axs[1].legend(loc = 'center left')
            
            plt.subplots_adjust(hspace=.0)
            
            
            for aa in self.data['dircomp']['Attempt'].keys():
            
            
                b2fstate = self.data['b2fstate'][aa]
                
                ne_pro = b2fstate['ne'].transpose()
                Te_J = b2fstate['te'].transpose()
                ev = 1.6021766339999999 * pow(10, -19)
                te_pro = Te_J / ev
                
                data = self.data['ft44'][aa]['dab2']
                neu_pro = np.transpose(data[:, :, 0])
                leftcut = self.data['b2fgeo']['leftcut'][0]
                rightcut = self.data['b2fgeo']['rightcut'][0]
                
                core_ne_pro = ne_pro[:, int(leftcut)+1:int(rightcut)-1]
                core_te_pro = te_pro[:, int(leftcut)+1:int(rightcut)-1]
                core_neu_pro = neu_pro[:, int(leftcut)+1:int(rightcut)-1]
            
                mean_core_ne = np.mean(core_ne_pro, axis=1)
                std_core_ne = np.std(core_ne_pro, axis=1)
                # print(std_core_ne)
                
                mean_core_te = np.mean(core_te_pro, axis=1)
                std_core_te = np.std(core_te_pro, axis=1)
                
                mean_core_neu = np.mean(core_neu_pro, axis=1)
                std_core_neu = np.std(core_neu_pro, axis=1)
            
                psi_coord = self.data['psi']['psi_59_val'][aa][:, 1]
            
            
                # axs[0].set_yscale('log')
                axs[0].errorbar(psi_coord, mean_core_ne, yerr= std_core_ne, fmt = '-', color = color_dic[aa], label= 'ne_solps')
                
                # axs[1].set_yscale('log')
                axs[1].errorbar(psi_coord, mean_core_te, yerr= std_core_te, fmt = '-', color = color_dic[aa], label= 'te_solps')
                
                # axs[1].add_artist(anchored_text2)





