# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 10:50:23 2025

@author: ychuang
"""

from matplotlib.offsetbox import AnchoredText
from modify_transport_coefficient.transport_coefficient_adjust_method import transcoe_method
from load_directory.grab_attempt_number import grab_aptn_method
from SOLPS_input.input_setting import set_figdir
from load_coordinate.SOLPSplotter_geo import load_geometry
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np



class transport_coefficient_alignment:
    
    def __init__(self, DF, data, geo: load_geometry, gam: grab_aptn_method, tm: transcoe_method):
        
        self.DF = DF
        self.data = data
        self.geo = geo
        self.gam = gam
        self.tm = tm
    
    
    
    def transport_coe_align_plot(self, plot_transcoe, paper_transcoe, save_eps):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == True and withseries == False:
            trans_dic = {}
            jxa = self.data['b2mn']['org']['jxa']
            self.geo.calcpsi_1D(pol_loc= jxa, no_coord_avg_check = False)
            for aa in self.data['dircomp']['multi_shift']:
                trans_file_dir = self.data['dirdata']['simudir'][aa] + '/b2.transport.inputfile'
                
                trans_list = self.tm.load_transcoefile_method(trans_file_dir, plot= False)
                cod = trans_list['1'].T
                coki = trans_list['3'].T
                coke = trans_list['4'].T
                x= cod[:,0]  #the coordinate here is R-R_sep
                
                trans_dic[aa] = np.zeros([len(x), 4])
                trans_dic[aa][:, 0] = self.data['psi']['psi_{}_val'.format(jxa)][aa][:, 1]
                trans_dic[aa][:, 1] = cod[:, 1]
                trans_dic[aa][:, 2] = coki[:, 1]
                trans_dic[aa][:, 3] = coke[:, 1]
    
   
            log_flag = False
            coe_label_dic = {'1': 'particle transport coefficient', '2': 'ion thermal diffusivity'
                             ,'3': 'electron heat transport coefficient'}
            
            
            
            if plot_transcoe:
                for k in coe_label_dic.keys():
                    if log_flag:
                        plt.yscale('log')
                    plt.figure(figsize=(7,7))
                    color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green', 'dot7': 'blue'}
                    A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4', 'dot7': '2.8'}
                    for ab in self.data['dircomp']['multi_shift']: 
                        # plt.plot(trans_dic[ab][:, 0], trans_dic[ab][:, int(k)], 'o-', color= color_dic[ab],
                        #          label ='transport coefficient of modify {} m case'.format(self.data['dircomp']['shift_dic'][ab]))
                        plt.plot(trans_dic[ab][:, 0], trans_dic[ab][:, int(k)], 'o-', color= color_dic[ab], 
                                 label = 'A = {}'.format(A_dic[ab]))
                        plt.xlabel('$\Psi_N$')
                        plt.title('Radial particle diffusion coefficient at [$m^{2}/s$] outer midplane')
                        plt.legend() 
                    
                    if save_eps:
                        
                        fig_dir  = set_figdir()
                        plt.savefig('{}/{}.pdf'.format(fig_dir, coe_label_dic[k]), format='pdf')
                    
                    
                plt.show()
            else:
                pass
            
            
            if paper_transcoe:
                
                alphabat_list = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
                
                fig, axs = plt.subplots()
                
                label_list = ['1', '3']
                
                for i, k in enumerate(label_list):
                    if log_flag:
                        plt.yscale('log')
                    color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green', 'dot7': 'blue'}
                    A_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4', 'dot7': '2.8'}
                    
                    coe = coe_label_dic[k]
                    po = alphabat_list[i]
                    anchored_text = AnchoredText('{}{}'.format(po, coe), loc= 'upper left')
                    
                    if i == 0:
                        for ab in self.data['dircomp']['multi_shift']:
                            
                            
                            axs.plot(trans_dic[ab][:, 0], trans_dic[ab][:, int(k)],
                                     '-', color= color_dic[ab], label = 'A = {}'.format(A_dic[ab]))
                    else:
                        for ab in self.data['dircomp']['multi_shift']:
                            
                            
                            axs.plot(trans_dic[ab][:, 0], trans_dic[ab][:, int(k)],
                                     '--', color= color_dic[ab])
                        
                    
                    
                    
                    # axs[i].add_artist(anchored_text)
                    
                    
                    
                    
                    
                axs.set_xlabel('$\psi_N$')
                axs.legend(loc= 'lower left')
                axs.set_yscale('log')
                
                # axs[2, 1].set_xlabel('poloidal angle')
                
                plt.subplots_adjust(hspace=.0)
                if save_eps:
                    
                    fig_dir  = set_figdir()
                    plt.savefig('{}/{}.eps'.format(fig_dir, coe_label_dic[k]), format='eps')
                
            
        else:
            print('transport_coe_align_plot is not there yet')
    
    
    
    def transport_coe_compare_plot(self, file_loc_list, plot_compare):
        trans_dic = {}
        psi_1d_dic = {}
        for fl in file_loc_list:
            fname = fl.rsplit("/",1)[1]
            psi_1d_dic[fname] = self.geo.calcpsi_block_method(file_loc = fl, shift = 0)
            
        
        self.data['psi_1d'] = psi_1d_dic    
            
            
        for an in file_loc_list:
            aa = an.rsplit("/",1)[1]
            trans_file_dir = an + '/b2.transport.inputfile'
            
            trans_list = self.tm.load_transcoefile_method(trans_file_dir, plot= False)
            cod = trans_list['1'].T
            coki = trans_list['3'].T
            coke = trans_list['4'].T
            x= cod[:,0]  #the coordinate here is R-R_sep
            
            trans_dic[aa] = np.zeros([len(x), 4])
            trans_dic[aa][:, 0] = self.geo.calcpsi_block_method(file_loc = an, shift = 0)
            trans_dic[aa][:, 1] = cod[:, 1]
            trans_dic[aa][:, 2] = coki[:, 1]
            trans_dic[aa][:, 3] = coke[:, 1]
    
   
        log_flag = False
        coe_label_dic = {'1': 'particle transport coefficient', '2': 'ion thermal diffusivity'
                         ,'3': 'electron heat transport coefficient'}
        
        note_list = ['leakage', 'decay length']
        
        
        if plot_compare:
            for k in coe_label_dic.keys():
                if log_flag:
                    plt.yscale('log')
                plt.figure()
                for ab in trans_dic.keys(): 
                    # plt.plot(trans_dic[ab][:, 0], trans_dic[ab][:, int(k)], 'o-', color= color_dic[ab],
                    #          label ='transport coefficient of modify {} m case'.format(self.data['dircomp']['shift_dic'][ab]))
                    plt.plot(trans_dic[ab][:, 0], trans_dic[ab][:, int(k)], 'o-', 
                             label= '{}'.format(ab))
                    plt.xlabel('psiN')
                    plt.title('radial {} coefficient at outer midplane'.format(coe_label_dic[k]))
                    plt.legend() 
            plt.show()
            
            
                
        else:
            print('transport_coe_align_plot is not there yet')
    
    
    
    def align_transcoe_method(self, std_file_loc, input_file_loc, itername, plot_align, log_flag):
        
        # one_list, n , filename_a, shift_a = tl.mast_tranco_dir('one')
        # print(one_list)
        # org_list = tl.mast_std_dir()

        std_trans = self.tm.load_transcoefile_method(std_file_loc, plot= False)
        cod = std_trans['1'].T
        coki = std_trans['3'].T
        coke = std_trans['4'].T
        x= cod[:,0]  #the coordinate here is R-R_sep
        yd= cod[:,1]
        yki = coki[:,1]
        yke = coke[:,1]

        input_trans = self.tm.load_transcoefile_method(input_file_loc, plot= False)
        ond = input_trans['1'].T
        onki = input_trans['3'].T
        onke = input_trans['4'].T
        onx= ond[:,0]  #the coordinate here is R-R_sep
        fd = ond[:,1]
        fki = onki[:,1]
        fke = onke[:,1]

        d_func = interpolate.interp1d(x, yd, fill_value = 'extrapolate')
        ond[:,1] = d_func(onx)
        ki_func = interpolate.interp1d(x, yki, fill_value = 'extrapolate')
        onki[:,1] = ki_func(onx)
        ke_func = interpolate.interp1d(x, yke, fill_value = 'extrapolate')
        onke[:,1] = ke_func(onx)

        n = str(self.gam.s_number(input_file_loc)[0])
        simu_dir = input_file_loc.rsplit("/",2)[0]
        filename = self.gam.s_number(std_file_loc)[1]
        char = filename.split('_')[2]
        print(char)
        # k = str(ss.s_number(file_loc, series_flag= None))
        print(simu_dir)

        
        # tcam.Generate_transcoefile_method(cod, CoeffID=1, SpeciesID=2, M=[1])
        self.tm.Write_transcoefile_method(file='{}/b2.transport.inputfile_align_{}_{}_{}'.format(simu_dir, char, itername, n), 
                                       points= input_trans ,M_1 = True, M=[1])

        # log_flag = True
        specieslist = ['1','3','4']
        d = self.tm.transport_coe_unit()
        
        if plot_align:
            for k in specieslist:
                if log_flag:
                    plt.yscale('log')
                else:
                    pass
                
                plt.figure()
                plt.plot(std_trans[k][0,:], std_trans[k][1,:], 'o-',color = 'blue', label ='orgin_case transport coefficient')
                plt.plot(input_trans[k][0,:], input_trans[k][1,:], 'o-', color = 'orange', label ='{}_case transport coefficient'.format(itername))
                plt.xlabel('Radial coordinate: $R- R_{sep}$')
                # plt.ylabel(d[k][1])
                plt.title(d[k][0])
                plt.legend()
                

            plt.show()
        


    def align_transco(self, plot_align, log_flag):
        # folder = '72_n100000_m12n8e3_nts5_a'
        
        if self.DF.withshift == True and self.DF.withseries == False:
            simudir = self.data['dirdata']['simudir']['org']
            stdfileloc = '{}/b2.transport.inputfile'.format(simudir)
            for shiftname in self.data['dircomp']['multi_shift']:
                if shiftname == 'org':
                    pass
                else:
                    input_simudir = self.data['dirdata']['simudir'][shiftname]
                    inpfileloc = '{}/b2.transport.inputfile'.format(input_simudir)
                    
                    self.align_transcoe_method(std_file_loc = stdfileloc, 
                                    input_file_loc = inpfileloc, itername = shiftname,
                                    plot_align = plot_align, log_flag = log_flag)
                    
                
        
        else:
            print('align_transco is not there yet!')