# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 20:09:09 2024

@author: ychuang
"""

import matplotlib.pyplot as plt
from SOLPSplotter_iout_dataprocess import iout_process
import numpy as np
from matplotlib.offsetbox import AnchoredText


class iout_flux(iout_process):
    def __init__(self, DefaultSettings, loadDS):
        iout_process.__init__(self, DefaultSettings, loadDS)
        
    
    
    def load_b2fstate_fna(self):
        
        if self.withshift == False and self.withseries == False:
            
            fna_data = self.data['b2fstate']['fna']
            fna_a = fna_data[:, :, 0, 0]
            fna_b = fna_data[:, :, 0, 1]
            fna_c = fna_data[:, :, 1, 0]
            fna_d = fna_data[:, :, 1, 1]
            
            fna_dic = {'fna_a': fna_a, 'fna_b': fna_b, 'fna_c': fna_c, 'fna_d': fna_d}
            
            self.data['fna_check'] = fna_dic
        
        elif self.withshift == False and self.withseries == True:
            
            if self.series_flag == 'two_compare':
                
                fna_seriesdic = {}
                
                for aa in list(self.data['dircomp']['Attempt'].keys()):
                    
                    fna_data = self.data['b2fstate'][aa]['fna']
                    fna_a = fna_data[:, :, 0, 0]
                    fna_b = fna_data[:, :, 0, 1]
                    fna_c = fna_data[:, :, 1, 0]
                    fna_d = fna_data[:, :, 1, 1]
                    
                    fna_dic = {'fna_a': fna_a, 'fna_b': fna_b, 'fna_c': fna_c, 'fna_d': fna_d}
        
                    fna_seriesdic[aa] = fna_dic
                
                self.data['fna_check'] = fna_seriesdic
                
                
                
        else:
            print('load_b2fstate_fna function is not there yet!')
    
    def iout_geo_ratio(self, iterlist):
        
        f_tuple5 = [('hx.dat', 'hx'),('vol.dat', 'sqrt_g')]
        f_tuple6 = [('hy.dat', 'hy'), ('vol.dat', 'sqrt_g')]
        f_tuple7 = [('hy1.dat', 'hy1'), ('vol.dat', 'sqrt_g')]
        
        
            
        f_list2 = [f_tuple5, f_tuple6, f_tuple7]
        
        qu_dic = {}
        
        for tuples in f_list2:

            data, qu = self.load_iout_ratio(file_tuple = tuples, itername = iterlist)
            
            qu_dic[qu] = data
        
        return qu_dic
    
    
    
    def load_fluxes_iout(self):
        
        
        
        fcoe_list = []
        
        if self.withshift == False and self.withseries == False:
            
            qu_dic = self.iout_geo_ratio(iterlist = None)
        
            for aa in list(qu_dic.keys()):
                
                self.data['iout_data'][aa] = qu_dic[aa]
                
                print(aa)
                fcoe_list.append(aa)
        
        elif self.withshift == False and self.withseries == True:
            
            scan = list(self.data['dircomp']['Attempt'].keys())
            
            qu_seriesdic = {}
            
            for tt in scan:
                
                qu_dic = self.iout_geo_ratio(iterlist = tt)
                
                qu_seriesdic[tt] = qu_dic
            
            
            iout_datadic = {}
            
            
            for aa in list(qu_dic.keys()):
                
                self.data['iout_data'][aa] = {}
                
                for tt in scan:
                    
                    self.data['iout_data'][aa][tt] = qu_seriesdic[tt][aa]
                    
                    # print(tt)
                
                
            fcoe_list = list(qu_dic.keys())
        
            

        ioutq_dic = {'hz': [('hz.dat', 'hz')], 'hy1': [('hy1.dat', 'hy1')], 'mag': [('bbx.dat', 'bx'), ('bbz.dat', 'bz'), ('bb.dat', 'B')],
                  'flux_no_psch': [('b2tfnb_fnbx001.dat', 'pol_flux_no_psch'), 
                                  ('b2tfnb_fnby001.dat', 'rad_flux_no_psch')]}
        
        
        qu_list_dic = {}
        
        for dickey in ioutq_dic:
            
            tuple_list = ioutq_dic[dickey]
            quantlist = self.general_iout_loader(tpl_list = tuple_list)
            
            qu_list_dic[dickey] = quantlist
            
            
        
        quwithder_dic = {'particle_flux': [('b2npc11_fnax001.dat', 'flux_x_0'), 
                                ('b2npc11_fnay001.dat', 'flux_y_0')],
             'heat_flux': [('b2nph9_fhex.dat', 'heatflux_x_0'), 
                                     ('b2nph9_fhey.dat', 'heatflux_y_0')],
             'ion_heat_flux': [('b2nph9_fhix.dat', 'ionheatflux_x_0'), 
                                     ('b2nph9_fhiy.dat', 'ionheatflux_y_0')]
            }
        
        
        quwithder_list_dic = {}
        
        for dickey in quwithder_dic:
            
            tuple_list = quwithder_dic[dickey]
            quantlist = self.general_iout_loader(tpl_list = tuple_list)
            
            
            
            if self.withshift == False and self.withseries == False:
                
                derivequant_list = self.derive_no_geo_quant(qu_list = quantlist, gcoe_list = fcoe_list, 
            pol_name = 'poloidal_{}'.format(dickey), rad_name = 'radial_{}'.format(dickey),
            iterlist = None)
            
            elif self.withshift == True and self.withseries == False:
                
                scan = self.data['dircomp']['multi_shift']
                
                derivequant_list = self.derive_no_geo_quant(qu_list = quantlist, gcoe_list = fcoe_list, 
    pol_name = 'poloidal_{}'.format(dickey), rad_name = 'radial_{}'.format(dickey),
                            iterlist = scan)
             
            elif self.withshift == False and self.withseries == True:
                
                    
                scan = list(self.data['dircomp']['Attempt'].keys())
                
                derivequant_list = self.derive_no_geo_quant(qu_list = quantlist, gcoe_list = fcoe_list, 
    pol_name = 'poloidal_{}'.format(dickey), rad_name = 'radial_{}'.format(dickey),
                            iterlist = scan)
                
            else:
                print('derive_no_geo_quant function is not there yet!')
                
            
            
            quwithder_list_dic[dickey] = quantlist
            quwithder_list_dic['derive {}'.format(dickey)] = derivequant_list
        
    
    
        geoscale_list_dic = {'geo_coe': fcoe_list}
        
        
        
        qunamelist_dic = qu_list_dic | quwithder_list_dic
        ioutlist_dic = qunamelist_dic | geoscale_list_dic
        
        
        self.data['iout_load_quant'] = ioutlist_dic
        
    
    def radial_cell_plot(self):
        
        sy = [0.000210521024697835, 0.0215047650857415, 0.0224393452735885, 0.023427682988235998,
    0.0244982167368385, 0.0256863306887055, 0.027034308347024, 0.028594232849626, 
    0.0304292536771715, 0.0326077559177885, 0.0351687338448725, 0.038036349055041496, 
    0.040985645438948004, 0.043680071833872, 0.045754398282541, 0.046901682164282, 
    0.0469318152180315, 0.045738931234544505, 0.043194431777765996, 0.039365793458189496, 
    0.034870733756757505, 0.030799151423874502, 0.02817372339416, 0.0268549705697585, 
    0.026970947747741503, 0.0141682235400465, 0.0497515894897005, 0.11138558777931501, 
    0.16960155730197501, 0.21825005723858, 0.24982056924682, 0.26505290834583, 
    0.270309320500165, 0.27057965173766496, 0.26934981397658997, 0.26718203095713, 
    0.26579232894353, 0.26486031458229, 0.263402926912145, 0.261527233794115, 
    0.25738305745592, 0.250095507768165, 0.23862941769967, 0.22403946457153, 
    0.21140725030527, 0.208297586514855, 0.22264897797628, 0.255393724934, 
    0.30266952243837497, 0.35968586986788, 0.41718709322201497, 0.467052462742065, 
    0.503368255414555, 0.5216879082206, 0.52390807024146, 0.5225831606086999, 
    0.537495829085225, 0.573888872291905, 0.61786966379477, 0.64316299021126, 
    0.63698191831342, 0.60998434653954, 0.57976106911706, 0.556645029323235, 
    0.531283208523925, 0.49006833490595497, 0.42401653612185497, 
    0.332874957194765, 0.233728089392925, 0.14811940660067, 0.08507361787591, 
    0.034563124058114, 0.0145969519046295, 0.029759570372418, 0.032589263502433005, 
    0.0378737114435935, 0.04729745439016, 0.0584218839119795, 0.06882144003815549, 
    0.0781077463538455, 0.08662059163542199, 0.094330255962412, 0.1005858779193275, 
    0.105420972910645, 0.108845355187335, 0.11039781322575, 0.10955158548803, 
    0.106248730403215, 0.10066176331530399, 0.0942752425688385, 0.08940494957696149, 
    0.0877127023159335, 0.0888393030788365, 0.09153622092792749, 0.09474816858198151, 
    0.097944968008587, 0.1007866712845905, 0.0010142869512735002]
        
        print('the length of sy list is:')
        print(len(sy))
        
        ioutgcoey = self.data['iout_data']['hy_divide_sqrt_g']
        ioutgcoey1 = self.data['iout_data']['hy_divide_sqrt_g']
        
        ioutgcoey_core = self.data['iout_data']['hy_divide_sqrt_g'][0, :]
        ioutgcoey1_core = self.data['iout_data']['hy1_divide_sqrt_g'][0, :]
        
        
        
        ioutsy = []
        for ii in ioutgcoey_core:
            ioutsy.append(1/ii)
        
        ioutsy1 = []
        for ib in ioutgcoey1_core:
            ioutsy1.append(1/ib)
        
        pol_list_a = []
        for i in range(96):
            pol_list_a.append('{}'.format(1 + i))
        
        
        
        fig, axs = plt.subplots()
        
        anchored_text = AnchoredText('(a){}'.format('rad cell area [$m^2$]'), loc='upper right')
        axs.plot(pol_list_a, sy[1:97], color = 'black', label= 'b2pl_radcell_area')
        axs.plot(pol_list_a, ioutsy, color = 'red', label= 'iout_radcell_area')
        axs.plot(pol_list_a, ioutsy1, color = 'blue', label= 'iouty1_radcell_area')
        axs.add_artist(anchored_text)
        axs.legend(loc='lower left', fontsize=10)
        
    
    def plot_radial_particleflux(self):
        
        b2fna = self.data['fna_check']['fna_d']
        outputfna = self.data['outputdata']['IonFlx']['D_1']
        ioutfna = self.data['iout_data']['flux_y_0']
        nogeo_ioutfna = self.data['iout_data']['radial_particle_flux']
        fpfna = self.data['b2fplasmf']['fna'][:, :, 1, 1]
        
        
        pol_list_a = []
        for i in range(48):
            pol_list_a.append('{}'.format(24 + i))
        
        
        b2fna_core = b2fna[25:73, 1]
        outputfna_core = outputfna[0, 25:73]
        ioutfna_core = ioutfna[0, 24:72]
        nogeo_ioutfna_core = nogeo_ioutfna[0, 24:72]
        fpfna_core = fpfna[25:73, 1]
        
        tot_particle_flux = sum(fpfna_core)
        print('the total particle flux is:')
        print(tot_particle_flux)
        
        
        
        b2fhe = self.data['b2fstate']['fhe'][:, :, 1]
        ioutfhe = self.data['iout_data']['heatflux_y_0']
        nogeo_ioutfhe = self.data['iout_data']['radial_heat_flux']
        fpfhe = self.data['b2fplasmf']['fhe'][:, :, 1]
        
        b2fhe_core = b2fhe[25:73, 1]
        ioutfhe_core = ioutfhe[0, 24:72]
        nogeo_ioutfhe_core = nogeo_ioutfhe[0, 24:72]
        fpfhe_core = fpfhe[25:73, 1]
        
        
        tot_heat_flux = sum(fpfhe_core)
        print('the total heat flux is:')
        print(tot_heat_flux)
        
        
        "b2pl data from 76_n900000_leakbsol_nts5_a"
        
        core_qe_dat = [ 201.593169786325, 768.94945454343, 1674.3401617572, 2590.42657722835, 
 3467.3839488935, 4216.5630555642, 4788.32494200055, 5203.034446893949, 5470.766505588201, 5637.1060473104, 
   5728.4740157359, 5767.9826797243, 5750.58373425685, 5650.6642311816995, 5489.1289623938, 
   5269.3740836949, 4997.775881337549, 4651.60466730775, 4271.76740635155, 3989.6376144459, 
   3940.22126647365, 4276.10699358925, 5057.88108788085, 6296.6001412634, 7994.5061247461, 
   10026.25960977345, 12116.028298046, 13913.43390966, 15105.659592105, 15704.231412956, 
   16127.768081341, 17131.4080702885, 19193.5477367535, 21778.522268839, 23531.9768739735, 
   23379.470467462, 21717.375552681, 19642.431893807, 17601.550906162498, 15292.473764461, 
   12434.183485932499, 9204.005185334649, 6156.62976261505, 3845.52155425525, 
   2304.2740817407002, 1299.8920843014498, 540.649026585805, 205.86765175645002]
        
        
        b2fhi = self.data['b2fstate']['fhi'][:, :, 1]
        ioutfhi = self.data['iout_data']['ionheatflux_y_0']
        nogeo_ioutfhi = self.data['iout_data']['radial_ion_heat_flux']
        fpfhi = self.data['b2fplasmf']['fhi'][:, :, 1]
        
        b2fhi_core = b2fhi[25:73, 1]
        ioutfhi_core = ioutfhi[0, 24:72]
        nogeo_ioutfhi_core = nogeo_ioutfhi[0, 24:72]
        fpfhi_core = fpfhi[25:73, 1]
        
        
        tot_ionheat_flux = sum(fpfhi_core)
        print('the total ion heat flux is:')
        print(tot_ionheat_flux)
        
        
        
        
        fig, axs = plt.subplots()
        
        anchored_text = AnchoredText('(a){}'.format('Particle flux [$s^{-1}$]'), loc='upper right')
        axs.plot(pol_list_a, b2fna_core, color = 'black', label= 'b2fstate_radial flux')
        axs.plot(pol_list_a, outputfna_core, color = 'red', label= 'output_radial flux')
        axs.plot(pol_list_a, ioutfna_core, color = 'blue', label= 'iout_radial flux')
        axs.plot(pol_list_a, fpfna_core, color = 'green', label= 'fplasmf_radial flux')
        axs.add_artist(anchored_text)
        axs.legend(loc='lower left', fontsize=10)
        
        
        fig, axs = plt.subplots()
        
        anchored_text = AnchoredText('(b){}'.format('Particle flux [$s^{-1}$]'), loc='upper right')
        axs.plot(pol_list_a, nogeo_ioutfna_core, color = 'black', label= 'iout_radial flux')
        axs.add_artist(anchored_text)
        axs.legend(loc='lower left', fontsize=10)
        
        
        fig, axs = plt.subplots()
        
        anchored_text = AnchoredText('(c){}'.format('Heat flux [$s^{-1}$]'), loc='upper right')
        axs.plot(pol_list_a, ioutfhe_core, color = 'red', label= 'iout_radheatflux')
        axs.plot(pol_list_a, b2fhe_core, color = 'black', label= 'b2fstate_radheatflux')
        axs.plot(pol_list_a, core_qe_dat, color = 'blue', label= 'b2pl_radheatflux')
        axs.plot(pol_list_a, fpfhe_core, color = 'green', label= 'fplasmf_radheatflux')
        axs.add_artist(anchored_text)
        axs.legend(loc='lower left', fontsize=10)
        
        
        b2plot_core_qe = [14228.542429227044, 15455.776638103527, 15031.928233610624, 
    15273.601365676755, 15887.207512175497, 16878.36621410578, 18065.54386399316, 
    19248.446325367364, 20218.69150342566, 20928.56855583475, 21440.34160985567, 
    21701.087847985866, 21711.760568305865, 21452.549132327647, 20988.747071423804, 
    20472.88634993912, 19983.46921916893, 19493.0059845433, 19067.030956002327, 
    18871.810728747016, 18916.30782862022, 19205.599021634884, 19804.249650957285, 
    20803.548670961474, 22226.355813428665, 24033.005269503297, 25941.471814349927, 
    27640.66617232631, 28955.35693673285, 29975.166073922464, 30861.62987447879, 
    31872.634434103034, 33444.71144753408, 35247.7610489595, 36587.89021153712, 
    36703.50726652556, 35603.168631923654, 33880.218835185326, 31620.781609354053, 
    28784.033673769576, 25372.346263340278, 21706.712831335375, 18495.322731691133, 
    16452.971331958594, 15556.868168889067, 15279.614488683168, 15642.365709672673, 
    14103.468525587044]
        
        
        
        
        fig, axs = plt.subplots()
        
        anchored_text = AnchoredText('(d){}'.format('Heat flux [$s^{-1}$]'), loc='upper right')
        axs.plot(pol_list_a, nogeo_ioutfhe_core, color = 'black', label= 'iout_radheatflux')
        axs.plot(pol_list_a, b2plot_core_qe, color = 'green', label= 'iout_radheatflux')
        axs.add_artist(anchored_text)
        axs.legend(loc='lower left', fontsize=10)
        
        
        fig, axs = plt.subplots()
        
        anchored_text = AnchoredText('(e){}'.format('Ion Heat flux [$s^{-1}$]'), loc='upper right')
        axs.plot(pol_list_a, b2fhi_core, color = 'black', label= 'b2fstate_ionradheatflux')
        axs.plot(pol_list_a, ioutfhi_core, color = 'blue', label= 'iout_ionradheatflux')
        axs.plot(pol_list_a, fpfhi_core, color = 'green', label= 'fplasmf_ionradheatflux')
        axs.add_artist(anchored_text)
        axs.legend(loc='lower left', fontsize=10)
    
    
    
    def plot_flux_compare(self):
        
        
        
        
        # outputfna = self.data['outputdata']['IonFlx']['D_1']
        # ioutfna = self.data['iout_data']['flux_y_0']
        # nogeo_ioutfna = self.data['iout_data']['radial_particle_flux']
        # fpfna = self.data['b2fplasmf']['fna'][:, :, 1, 1]
        
        
        pol_list_a = []
        for i in range(48):
            pol_list_a.append('{}'.format(24 + i))
            
        
        fig, axs = plt.subplots()
        
        anchored_text = AnchoredText('(a){}'.format('Particle flux [$s^{-1}$]'), loc='upper right')
        
        for aa in list(self.data['dircomp']['Attempt'].keys()):
            
            b2fna = self.data['fna_check'][aa]['fna_d']
            b2fna_core = b2fna[25:73, 1]
                   
            axs.plot(pol_list_a, b2fna_core, label= '{}_particleflux'.format(aa))
        
        axs.add_artist(anchored_text)
        axs.legend(loc='lower left', fontsize=10)
        
        
        
    

    
    
   