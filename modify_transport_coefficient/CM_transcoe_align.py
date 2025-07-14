# -*- coding: utf-8 -*-
"""
Created on Sun Jul  6 22:35:18 2025

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



class CrossMachine_transcoe_align:
    
    def __init__(self, DF, data, geo: load_geometry, gam: grab_aptn_method, tm: transcoe_method):
        
        self.DF = DF
        self.data = data
        self.geo = geo
        self.gam = gam
        self.tm = tm
    
    
    def CM_align_transcoe_method(self, std_file_loc, input_file_loc, device, plot_align, log_flag):
        
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
        char = self.data['{}_dircomp'.format(device)]['Shot']
        print(char)
        # k = str(ss.s_number(file_loc, series_flag= None))
        print(simu_dir)

        
        # tcam.Generate_transcoefile_method(cod, CoeffID=1, SpeciesID=2, M=[1])
        self.tm.Write_transcoefile_method(file='{}/b2.transport.inputfile_align_{}_{}'.format(simu_dir, char, n), 
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
                plt.plot(input_trans[k][0,:], input_trans[k][1,:], 'o-', color = 'orange', label ='{}_case transport coefficient'.format(device))
                plt.xlabel('Radial coordinate: $R- R_{sep}$')
                # plt.ylabel(d[k][1])
                plt.title(d[k][0])
                plt.legend()
                

            plt.show()
    
    def CM_align_transco(self, plot_align, log_flag):

        

        if self.DF.withshift == False and self.DF.withseries == False:
                                    
            simudir = self.data['mast_dirdata']['simudir']
            stdfileloc = '{}/b2.transport.inputfile'.format(simudir)
                        
            input_simudir = self.data['mastu_dirdata']['simudir']
            inpfileloc = '{}/b2.transport.inputfile'.format(input_simudir)
            
            
            self.CM_align_transcoe_method(std_file_loc = stdfileloc, input_file_loc = inpfileloc,
                            device = 'mastu', plot_align = plot_align, log_flag = log_flag)
                
                    
        else:
            print('align_transco is not there yet!')
    
    
    
    def CM_transcoe_align_plot(self, plot_transcoe, paper_transcoe, save_eps):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == False:
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