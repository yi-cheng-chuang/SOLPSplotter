# -*- coding: utf-8 -*-
"""
Created on Thu May 23 12:45:06 2024

@author: ychuang
"""

import SOLPS_set as sps
import matplotlib.pyplot as plt
from SOLPSplotter_contour import PlotContour
import SOLPS_transcoe_adj as sta
import numpy as np
from matplotlib.offsetbox import AnchoredText
import scipy.stats as stats
from matplotlib.colors import LogNorm
from matplotlib import cm
from numpy import ma
from scipy.optimize import curve_fit


class Plotiout(PlotContour):
    def __init__(self, DefaultSettings, loadDS):
        PlotContour.__init__(self, DefaultSettings, loadDS)
        
        self.Publish = DefaultSettings['Publish']
        self.data['DefaultSettings']['Publish'] = self.Publish
    
    
    def set_plot(self):
        if self.Publish == 'b2plottersetting':
            plt.rcParams.update({'font.weight': 'normal'})
            plt.rc('lines', linewidth= 3, markersize= 6)
            plt.rcParams.update({'font.size': 10})
            plt.rcParams.update({'figure.facecolor':'w'})
            plt.rcParams.update({'mathtext.default': 'regular'})
  
        else:
            print('Publish setting is incorrect or add another setting')


    def weight_generater(self, pol_list, psi_dic, input_dat, itername, nnp):
        
        sk = int(pol_list[0])
        sd = int(pol_list[-1]) + 1
            
        # psi_st = 17
        # psi_ed = 38
        
        # psi_dic = {'st': psi_st, 'ed': psi_ed}
        # np = 1
        
        psi_ref = self.data['psi']['psival']['org'][psi_dic['st']:psi_dic['ed'], sk:sd]
        avg_psi = []
        
        print('len pol_list is :{}, len psi ref is: {}'.format(len(pol_list), len(psi_ref[0, :])))
        
        for xa in range(len(pol_list)):
            avg_psi.append(0.5*(float(psi_ref[nnp + 1, xa]) + float(psi_ref[nnp + 2, xa])))
        # avg_psi = 0.5*(psi_ref[:, np + 1] + psi_ref[:, np + 2])
        
        psi_w_dic = {}
        
        for aa in self.data['dircomp']['multi_shift']:
            
            psi_dat = self.data['psi']['psival'][aa][psi_dic['st']:psi_dic['ed'], sk:sd]
            
            if aa == 'org':
                
                # ln = int(len(pol_list))
                # print('look at type: {}'.format(type(ln)))
                w_psi = []
                for xa in range(len(pol_list)):
                    w_psi.append(0.5)
                
                
            else:
                
                # w_psi = np.zeros(len(pol_list))
                # print('avg shape: {}'.format(avg_psi.shape))
                w_psi = []
                for xa in range(len(pol_list)):
                    w_psi.append((avg_psi[xa] - psi_dat[nnp + 1, xa])/(psi_dat[nnp + 2, xa] - psi_dat[nnp + 1, xa]))
                
            psi_w_dic[aa] = w_psi
            
        
        
        if input_dat.shape == (98, 38):
            
            plot_dat = input_dat[sk:sd, psi_dic['st']:psi_dic['ed']]
            
            wei = psi_w_dic[itername]
            avg_dat = []
            
            for wa in range(len(pol_list)):
                avg_dat.append((1- wei[wa])*plot_dat[wa, nnp + 1] + wei[wa]*plot_dat[wa, nnp + 2])
                
            
        
        elif input_dat.shape == (36, 96):
            
            wei = psi_w_dic[itername]
            avg_dat = []
            
            
            plot_dat = input_dat[psi_dic['st']:psi_dic['ed'], sk:sd]
            
            
            for wa in range(len(pol_list)):
                avg_dat.append((1- wei[wa])*plot_dat[nnp + 0, wa] + wei[wa]*plot_dat[nnp + 1, wa])
        
        
        elif input_dat.shape == (38, 98):
            
            wei = psi_w_dic[itername]
            avg_dat = []
            
            
            plot_dat = input_dat[psi_dic['st']:psi_dic['ed'], sk:sd]
            
            
            for wa in range(len(pol_list)):
                avg_dat.append((1- wei[wa])*plot_dat[nnp + 1, wa] + wei[wa]*plot_dat[nnp + 2, wa])
        
        return avg_dat
    
    
    
    
    
    
    def flux_poloidal_sub(self, input_dat, itername, pol_list, psi_dic, ang_list, art_text,
                          axs, color_dic, A_dic, no_A_label, input_ls, nnp):
    
            
            # np, psi_st, psi_ed, sk, sd, psi_w_dic = weight_generater(pol_list)
            
            
            avg_dat = self.weight_generater(pol_list = pol_list, input_dat = input_dat, 
                                itername = itername, nnp = nnp, psi_dic = psi_dic)
            
            axs.add_artist(art_text)
            
            
            if no_A_label:
                
                axs.plot(ang_list, avg_dat, linestyle = input_ls, 
                    color= color_dic[itername])
            else:
                
                axs.plot(ang_list, avg_dat, linestyle = input_ls, 
                    color= color_dic[itername], label= 'A = {}'.format(A_dic[itername]))
                
    
    def general_iout_loader(self, tpl_list):
        
        qu_list = []
        
        for tpl_qu in tpl_list:
            
            quant = self.load_iout(filename = tpl_qu[0], simple_quant = tpl_qu[1])
            print(quant)
            # print(flux_qu.split('_'))
            qu_list.append(quant)
        
        return qu_list
    
    
    def derive_no_geo_quant(self, qu_list, gcoe_list, pol_name, rad_name):
        
        res_qu_list = []
        
        for item in qu_list:
            
            res_qu_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                if item.split('_')[1] == 'x':
                    data, res_qu = self.load_iout_multi(name1 = gcoe_list[0], name2 = item, 
                                        input_name= pol_name, itername = aa)
                    print('flux check: {}'.format(gcoe_list[0]))
                    print('flux check: {}'.format(item))
                    res_qu_dic[aa] = data
                
                elif item.split('_')[1] == 'y':
                    data, res_qu = self.load_iout_multi(name1 = gcoe_list[1], name2 = item,
                                        input_name= rad_name, itername = aa)
                    
                    print('flux check: {}'.format(gcoe_list[1]))
                    print('flux check: {}'.format(item))
                    
                    res_qu_dic[aa] = data
          
                else:
                    print('a bug! flag is {}'.format(item.split('_')[1]))
            
            self.data['iout_data'][res_qu] = res_qu_dic
            
        
            print(res_qu)
            res_qu_list.append(res_qu)
            
            
        
        return res_qu_list
    
    
    def flux_iout_loader(self):
        
        
        f_tuple5 = [('hx.dat', 'hx'),('vol.dat', 'sqrt_g')]
        f_tuple6 = [('hy.dat', 'hy'), ('vol.dat', 'sqrt_g')]
        
        
            
        f_list2 = [f_tuple5, f_tuple6]
        fcoe_list = []
        
        
        for tuples in f_list2:
            
            qu_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                data, qu = self.load_iout_ratio(file_tuple = tuples, itername = aa)
                qu_dic[aa] = data
            
            self.data['iout_data'][qu] = qu_dic
                
            print(qu)
            fcoe_list.append(qu)
        
        # flux_qu_list = self.general_iout_loader(tpl_list = flux_tuple)
        
        # res_qu_list = self.derive_no_geo_quant(qu_list = flux_qu_list, gcoe_list = fcoe_list, 
        #                             pol_name = 'poloidal_flux', rad_name = 'radial_flux')
        
        
        ioutq_dic = {'hz': [('hz.dat', 'hz')], 'mag': [('bbx.dat', 'bx'), ('bbz.dat', 'bz'), ('bb.dat', 'B')],
                  'flux_no_psch': [('b2tfnb_fnbx001.dat', 'pol_flux_no_psch'), 
                                        ('b2tfnb_fnby001.dat', 'rad_flux_no_psch')],
                  'ni': [('b2npc11_na001.dat', 'ni')], 
                  'vp and p': [('b2npmo_ua001.dat', 'vp'), ('b2nppo_po.dat', 'e_potential')],
                  'tcoe_density': [('b2tqna_dna0001.dat', 'tcoe_n')],
                  'source': [('b2npc11_sna001.dat', 'source')]}
        
        
        qu_list_dic = {}
        
        for dickey in ioutq_dic:
            
            tuple_list = ioutq_dic[dickey]
            quantlist = self.general_iout_loader(tpl_list = tuple_list)
            
            qu_list_dic[dickey] = quantlist
            
            
        
        quwithder_dic = {'flux': [('b2npc11_fnax001.dat', 'flux_x_0'), 
                                ('b2npc11_fnay001.dat', 'flux_y_0')],
            'psch_flux': [('b2tfnb_fnbPSchx001.dat', 'psch_x'), 
                           ('b2tfnb_fnbPSchx001.dat','psch_y')],
            'vpara': [('b2tfnb_bxuanax001.dat', 'vpara_x')],
            'gradn': [('b2tfnb_dPat_mdf_gradnax001.dat', 'gradn_x'), 
                       ('b2tfnb_dPat_mdf_gradnay001.dat', 'gradn_y')],
            'tcoe': [('b2trno_cdnax001.dat', 'tcoe_x'), 
                     ('b2trno_cdnay001.dat', 'tcoe_y')],
            'corrdpc': [('b2tfnb_dpccornax001.dat', 'corrdpc_x')],
            
            }
        
        
        quwithder_list_dic = {}
        
        for dickey in quwithder_dic:
            
            tuple_list = quwithder_dic[dickey]
            quantlist = self.general_iout_loader(tpl_list = tuple_list)
            derivequant_list = self.derive_no_geo_quant(qu_list = quantlist, gcoe_list = fcoe_list, 
                                        pol_name = 'poloidal_{}'.format(dickey), 
                                        rad_name = 'radial_{}'.format(dickey))
            
            quwithder_list_dic[dickey] = quantlist
            quwithder_list_dic['derive {}'.format(dickey)] = derivequant_list
        
        
            
        qu_dic = {}
        for aa in self.data['dircomp']['multi_shift']:
            data, qu = self.load_iout_divide(name1 = 'gradn_x', 
                        name2 = 'tcoe_x', itername = aa, 
                        input_name = 'pnpx')
            qu_dic[aa] = data
            
        
        self.data['iout_data'][qu] = qu_dic
        
        
        qy_dic = {}
        for aa in self.data['dircomp']['multi_shift']:
            data_y, qy = self.load_iout_divide(name1 = 'gradn_y', 
                        name2 = 'tcoe_y', itername = aa, 
                        input_name = 'pnpy')
            
            qy_dic[aa] = data_y
    
        self.data['iout_data'][qy] = qy_dic
    
    
        geoscale_list_dic = {'geo_coe': fcoe_list}
        
        
        
        qunamelist_dic = qu_list_dic | quwithder_list_dic
        ioutlist_dic = qunamelist_dic | geoscale_list_dic
        
        
        self.data['iout_load_quant'] = ioutlist_dic
















