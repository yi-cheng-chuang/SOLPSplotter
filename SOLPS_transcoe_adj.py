# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:00:38 2024

@author: user
"""
from SOLPSplotter_geo import load_geometry
from matplotlib.offsetbox import AnchoredText
import SOLPS_set as ss
import transport_coefficient_adjust_method as tcam
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np



class transport_coefficient_adjustment(load_geometry):
    
    def __init__(self, DefaultSettings):
        load_geometry.__init__(self, DefaultSettings)
    
    
    def set_plot(self):
        
        plt.rcParams.update({'font.weight': 'normal'})
        plt.rc('lines', linewidth= 5, markersize= 9)
        plt.rcParams.update({'font.size': 16})
        plt.rcParams.update({'figure.facecolor':'w'})
        plt.rcParams.update({'mathtext.default': 'regular'})
  
    
    
    
    def mod_transco_method(self,file_loc, withmod, de_SOL, ki_SOL, ke_SOL, log_flag):
        
        trans_list = tcam.load_transcoefile_method(file_loc, plot= False)
        cod = trans_list['1'].T
        coki = trans_list['3'].T
        coke = trans_list['4'].T
        x= cod[:,0]  #the coordinate here is R-R_sep
        yd= cod[:,1]
        yki = coki[:,1]
        yke = coke[:,1]

        m = len(yd)
        if withmod:
            mod_y = np.zeros(m)
            for j in range(m):
                if j<= de_SOL:
                    mod_y[j] = cod[j,1]
                else:
                    mod_y[j] = 12.0
            cod[:,1] = mod_y

            mod_yki = np.zeros(m)
            for j in range(m):
                if j<= ki_SOL:
                    mod_yki[j] = coke[j,1]  
                else:
                    mod_yki[j] = 18.0
            coki[:,1] = mod_yki

            mod_yke = np.zeros(m)
            for j in range(m):
                if j<= ke_SOL:
                    mod_yke[j] = coke[j,1]  
                else:
                    mod_yke[j] = 18.0
            coke[:,1] = mod_yke
        else:
            pass


        tcam.Generate_transcoefile_method(cod, CoeffID=1, SpeciesID=2, M=[1])
        
        shift = 'org'
        n = str(ss.s_number(file_loc, series_flag= None)[0])
        simu_dir = file_loc.rsplit("/",2)[0]
        # k = str(ss.s_number(file_loc, series_flag= None))
        print(simu_dir)
        
        tcam.Write_transcoefile_method(file = '{}/b2.transport.inputfile_mod_{}{}'.format(simu_dir, shift, n), points= trans_list ,M_1 = True, M=[1])

        log_flag = False
        specieslist = ['1','3','4']
        transcoe_unit = tcam.transport_coe_unit()

        for k in specieslist:
            if log_flag:
                plt.yscale('log')
            else:
                pass
            plt.figure(figsize=(7,7))
            plt.plot(trans_list[k][0,:], trans_list[k][1,:], 'o-', color = 'orange')
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
            # plt.ylabel(transcoe_unit[k][1])
            plt.title(transcoe_unit[k][0])

        plt.show()
        
        
        
        
    def mod_transco(self, withmod, de_SOL, ki_SOL, ke_SOL, log_flag):
        # folder = '72_n100000_m12n8e3_nts5_a'
        if self.withshift == False and self.withseries == False:
            simudir = self.data['dirdata']['simudir']
            fileloc = '{}/b2.transport.inputfile_new'.format(simudir)
            self.mod_transco_method(file_loc = fileloc, withmod = withmod, de_SOL = de_SOL, 
                                    ki_SOL = ki_SOL, ke_SOL = ke_SOL, log_flag = log_flag)
        
        elif self.withshift == True and self.withseries == False:
            simudir = self.data['dirdata']['simudir']['org']
            fileloc = '{}/b2.transport.inputfile_new'.format(simudir)
            self.mod_transco_method(file_loc = fileloc, withmod = withmod, de_SOL = de_SOL, 
                                    ki_SOL = ki_SOL, ke_SOL = ke_SOL, log_flag = log_flag)
        
        
        
    def transport_coe_align_plot(self, plot_transcoe, paper_transcoe, save_eps):
        if self.withshift == True and self.withseries == False:
            trans_dic = {}
            jxa = self.data['b2mn']['org']['jxa']
            self.calcpsi_1D(pol_loc= jxa, no_coord_avg_check = False)
            for aa in self.data['dircomp']['multi_shift']:
                trans_file_dir = self.data['dirdata']['simudir'][aa] + '/b2.transport.inputfile'
                
                trans_list = tcam.load_transcoefile_method(trans_file_dir, plot= False)
                cod = trans_list['1'].T
                coki = trans_list['3'].T
                coke = trans_list['4'].T
                x= cod[:,0]  #the coordinate here is R-R_sep
                
                trans_dic[aa] = np.zeros([len(x), 4])
                trans_dic[aa][:, 0] = self.data['psi']['psi_{}_val'.format(jxa)][aa][:, 2]
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
                        
                        fig_dir  = ss.set_figdir()
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
                    
                    fig_dir  = ss.set_figdir()
                    plt.savefig('{}/{}.eps'.format(fig_dir, coe_label_dic[k]), format='eps')
                
            
        else:
            print('transport_coe_align_plot is not there yet')
    




    def transport_coe_compare_plot(self, file_loc_list, plot_compare):
        trans_dic = {}
        psi_1d_dic = {}
        for fl in file_loc_list:
            fname = fl.rsplit("/",1)[1]
            psi_1d_dic[fname] = self.calcpsi_block_method(file_loc = fl, 
                                                                   shift = 0)
            
        
        self.data['psi_1d'] = psi_1d_dic    
            
            
        for an in file_loc_list:
            aa = an.rsplit("/",1)[1]
            trans_file_dir = an + '/b2.transport.inputfile'
            
            trans_list = tcam.load_transcoefile_method(trans_file_dir, plot= False)
            cod = trans_list['1'].T
            coki = trans_list['3'].T
            coke = trans_list['4'].T
            x= cod[:,0]  #the coordinate here is R-R_sep
            
            trans_dic[aa] = np.zeros([len(x), 4])
            trans_dic[aa][:, 0] = self.calcpsi_block_method(file_loc = an, shift = 0)
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
                plt.figure(figsize=(7,7))
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

        std_trans = tcam.load_transcoefile_method(std_file_loc, plot= False)
        cod = std_trans['1'].T
        coki = std_trans['3'].T
        coke = std_trans['4'].T
        x= cod[:,0]  #the coordinate here is R-R_sep
        yd= cod[:,1]
        yki = coki[:,1]
        yke = coke[:,1]

        input_trans = tcam.load_transcoefile_method(input_file_loc, plot= False)
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

        n = str(ss.s_number(input_file_loc, series_flag= None)[0])
        simu_dir = input_file_loc.rsplit("/",2)[0]
        filename = ss.s_number(std_file_loc, series_flag= None)[1]
        char = filename.split('_')[2]
        print(char)
        # k = str(ss.s_number(file_loc, series_flag= None))
        print(simu_dir)

        
        # tcam.Generate_transcoefile_method(cod, CoeffID=1, SpeciesID=2, M=[1])
        tcam.Write_transcoefile_method(file='{}/b2.transport.inputfile_align_{}_{}_{}'.format(simu_dir, char, itername, n), 
                                       points= input_trans ,M_1 = True, M=[1])

        # log_flag = True
        specieslist = ['1','3','4']
        d = tcam.transport_coe_unit()
        
        if plot_align:
            for k in specieslist:
                if log_flag:
                    plt.yscale('log')
                else:
                    pass
                plt.figure(figsize=(7,7))
                plt.plot(std_trans[k][0,:], std_trans[k][1,:], 'o-',color = 'blue', label ='orgin_case transport coefficient')
                plt.plot(input_trans[k][0,:], input_trans[k][1,:], 'o-', color = 'orange', label ='{}_case transport coefficient'.format(itername))
                plt.xlabel('Radial coordinate: $R- R_{sep}$')
                # plt.ylabel(d[k][1])
                plt.title(d[k][0])
                plt.legend()
                

            plt.show()
        


    def align_transco(self, plot_align, log_flag):
        # folder = '72_n100000_m12n8e3_nts5_a'
        
        if self.withshift == True and self.withseries == False:
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



    
    
    