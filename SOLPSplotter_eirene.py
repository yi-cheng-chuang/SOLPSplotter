# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 18:55:31 2024

@author: ychuang
"""

from SOLPSplotter_NTplot import NT_plot
import matplotlib.pyplot as plt
import SOLPS_set as ss
from matplotlib import colors, cm, ticker
import matplotlib.tri as tri
from matplotlib.colors import LogNorm
import fitting_method as fm
import numpy as np
from numpy import ma
from matplotlib.offsetbox import AnchoredText




class Eirene_contour(NT_plot):
    
    def __init__(self, DefaultSettings, loadDS):
        NT_plot.__init__(self, DefaultSettings, loadDS)
        
        
    
    def eirene_contourplot_method(self, simudir, data, plot33):
        
        # base_start_data = np.log10(np.nanmin(data))
        # base_end_data = np.log10(np.nanmax(data))
        # log_level_data = np.logspace(base_start_data, base_end_data, num = 20, base= 10)
        # print(log_level_data)
        
        
        # NORM_data = colors.LogNorm(np.nanmin(data), np.nanmax(data))
        
        # NORM_data = colors.Normalize(np.nanmin(data), np.nanmax(data))
        
        data_mask = ma.masked_where(data <= 0, data)
        
        datamap = np.abs(data_mask)
        CMAP = cm.viridis
        NORM_data = colors.LogNorm(np.nanmin(datamap), np.nanmax(datamap))
        # Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())
        
        
        if plot33:
            
            Nodes=np.fromfile('{}/fort.33'.format(simudir),sep=' ') #Alternatively use fort.33
            NN=int(Nodes[0])
            XNodes=Nodes[1:NN+1]
            YNodes=Nodes[NN+1:]
        
        
            numberlist = np.zeros(NN)
            for i in range(NN):
                numberlist[i] = i
        
            plt.figure()
            plt.scatter(XNodes[:500], YNodes[:500])
            
         
        Triangles = np.loadtxt('{}/fort.34'.format(simudir), 
            skiprows=1, usecols=(1,2,3)) #Alternatively use fort.34
        # print(Triangles -1)
    
        TP = tri.Triangulation(XNodes, YNodes, triangles= (Triangles -1))
    
        # CMAP = cm.viridis
    
        plt.figure()
        plt.tripcolor(TP, data, shading = 'flat', cmap = CMAP, norm = NORM_data)
        plt.title('Neutral density contour plot')
        # plt.title('Neutral density contour plot outerleg')
    
        SM_data= cm.ScalarMappable(NORM_data, CMAP)    
        plt.colorbar(SM_data)
    
    
    
    def twcontourplot_prep(self, ta, keylist_b, scan_style):
     
         if self.withshift == False and self.withseries == True:
             
             if self.series_flag == 'twin_scan':
                 
                 color_list = ['red', 'orange', 'green', 'blue', 'purple']
                 
                 color_dic = self.pair_dic(keys = keylist_b, values = color_list)
                 
                 iter_key = []
                 
                 
                 for tb in keylist_b:
                     
                     if scan_style == 'tempscan':
                         
                         it_in = (ta, tb)
                     
                     elif scan_style == 'denscan':
                         
                         it_in = (tb, ta)
                     
                     else:
                         print('twinscan_plot_method, please check the scan_style!')
                     
                     
                     nx = self.data['b2fgeo']['nx']
                     ny = self.data['b2fgeo']['ny']
                     
                     iter_key.append(it_in)
             
             
                 return iter_key, color_dic


    def twcontourplot(self, scan_style, plot_option, format_option):
        
        if self.withshift == False and self.withseries == True:
            
            # series_flag = self.DefaultSettings['series_flag']
            
            
            if self.series_flag == 'twin_scan':
                
                dircomp = self.data['dircomp']
                
                if scan_style == 'tempscan':
                    
                    key_a = 'denscan_list'
                    key_b = 'tempscan_list'
                
                elif scan_style == 'denscan':
                    
                    key_a = 'tempscan_list'
                    key_b = 'denscan_list'
                
                else:
                    print('twinscan_plot_method, please check the scan_style!')
                
                keylist_a = []
                
                # pol_list_a = []
                # for i in range(48):
                #     pol_list_a.append('{}'.format(25 + i))
                
                
                # self.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
                # ang_list = self.data['angle']['angle_list']
                
                for x in dircomp[key_a]:
                    keylist_a.append('{:.3f}'.format(x))
                
                for ta in keylist_a:
                    
                    keylist_b = []
                    
                    for x in dircomp[key_b]:
                        keylist_b.append('{:.3f}'.format(x))
                    
                    
                    iter_key, color_dic= self.twcontourplot_prep(ta = ta, 
                    keylist_b = keylist_b, scan_style = scan_style)
                    
                    
                    print('check:')
                    print(iter_key)
                    print(color_dic)
                    # print(label_dic)
                    
                    
                    self.twcontourplot_method(iterlist = iter_key, cl_dic = color_dic, 
                                scan_style = scan_style, 
                    plot_option = plot_option, format_option = format_option)
                    
             
            else:
                print('neteTS_plot, please check the series flag')
    
    
    def twcontourplot_method(self, iterlist, cl_dic, scan_style, plot_option, format_option):
        
        
        
        
        
        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']
        
        
        plt.subplots_adjust(hspace=.0)
        anchored_text_1 = AnchoredText('{}'.format('density pedestal width [mm]'), 
                                              loc='lower center')
        anchored_text_2 = AnchoredText('{}'.format('neutral penetration length [mm]'), 
                                              loc='lower right')
        anchored_text_3 = AnchoredText('{}'.format('dimensionless opaqueness'), 
                                              loc='lower right')
        anchored_text_4 = AnchoredText('{}'.format('neutral density'), 
                                              loc='lower right')
        
        for aa in iterlist:
            
            if format_option == '1x1':
                
                fig, axs = plt.subplots()
            
            # psi_coord, mid_ne_pro, mid_te_pro, mid_neu_pro = self.nete_midprof(itername = aa, 
            #                                         data_struc = dat_struc)
            if plot_option == 'Neuden contour':
                simu_dir = self.data['dirdata']['simudir'][aa[0]][aa[1]]
                dat = self.data['ft46'][aa[0]][aa[1]]['pdena'][:, 0]
                
                
            
            if self.series_flag == 'twin_scan':
                
                if scan_style == 'tempscan':
                    
                    ad = aa[1]
                    ap = aa[0]
                    
                
                elif scan_style == 'denscan':
                    
                    ad = aa[0]
                    ap = aa[1]
                
                else:
                    print('neteTSplot_method, please check scan_style')
            
            else:
                ad = aa
                
            plot_num = 0

            if scan_style == 'denscan':
                
                title_ap = float(ap)*pow(10, 5)
                
                if format_option == '1x1':
                    
                    Nodes=np.fromfile('{}/fort.33'.format(simu_dir),sep=' ') #Alternatively use fort.33
                    NN=int(Nodes[0])
                    XNodes=Nodes[1:NN+1]
                    YNodes=Nodes[NN+1:]
                    
                    
                    data_mask = ma.masked_where(dat <= 0, dat)
                    
                    datamap = np.abs(data_mask)
                    CMAP = cm.viridis
                    NORM_data = colors.LogNorm(np.nanmin(datamap), np.nanmax(datamap))
                    # Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())
                
                    Triangles = np.loadtxt('{}/fort.34'.format(simu_dir), 
                        skiprows=1, usecols=(1,2,3)) #Alternatively use fort.34
                    # print(Triangles -1)
                
                    TP = tri.Triangulation(XNodes, YNodes, triangles= (Triangles -1))
                
                    # CMAP = cm.viridis
                    
                    
                    
                    axs.tripcolor(TP, dat, shading = 'flat', cmap = CMAP, norm = NORM_data)
                    axs.set_title('Neutral density contour plot {}_{}'.format(aa[0], aa[1]))
                    # plt.title('Neutral density contour plot outerleg')
                
                    SM_data= cm.ScalarMappable(NORM_data, CMAP)
                    
                    plt.colorbar(SM_data)
                    
            elif scan_style == 'tempscan':
                
                title_ap = float(ap)*pow(10, 20)
                
                if format_option == '1x1':
                    
                    Nodes=np.fromfile('{}/fort.33'.format(simu_dir),sep=' ') #Alternatively use fort.33
                    NN=int(Nodes[0])
                    XNodes=Nodes[1:NN+1]
                    YNodes=Nodes[NN+1:]
                    
                    
                    data_mask = ma.masked_where(dat <= 0, dat)
                    
                    datamap = np.abs(data_mask)
                    CMAP = cm.viridis
                    NORM_data = colors.LogNorm(np.nanmin(datamap), np.nanmax(datamap))
                    # Lnorm = LogNorm(vmax = datamap.max(), vmin = datamap.min())
                
                    Triangles = np.loadtxt('{}/fort.34'.format(simu_dir), 
                        skiprows=1, usecols=(1,2,3)) #Alternatively use fort.34
                    # print(Triangles -1)
                
                    TP = tri.Triangulation(XNodes, YNodes, triangles= (Triangles -1))
                
                    # CMAP = cm.viridis
                    
                    
                    
                    axs.tripcolor(TP, dat, shading = 'flat', cmap = CMAP, norm = NORM_data)
                    axs.set_title('Neutral density contour plot {}_{}'.format(aa[0], aa[1]))
                    # plt.title('Neutral density contour plot outerleg')
                
                    SM_data= cm.ScalarMappable(NORM_data, CMAP)
                    
                    plt.colorbar(SM_data)
                    
                    
            
            else:
                print('neudenplot_method, please check the scan parameter')
        
            plot_num = plot_num + 1
            plt.tight_layout()
    
    
    
    def eirene_contour_plot(self):
        
        
        if self.withshift == False and self.withseries == False:
            
            simu_dir = self.data['dirdata']['simudir']
            dat = self.data['ft46']['pdena'][:, 0]
            abs_dat = np.absolute(dat)
            
            self.eirene_contourplot_method(simudir = simu_dir, data = abs_dat)
            
        
        
        elif self.withshift == False and self.withseries == True:
            
            print('we are working on it!')
            
            
            
        else:
            print('load_b2fplasmf function is not there yet!')
            
        
        
        
        

        

"""

base_start_eiratom = np.log10(np.nanmin(Eiratom))
base_end_eiratom = np.log10(np.nanmax(Eiratom))
log_level_eiratom = np.logspace(base_start_eiratom, base_end_eiratom, num=20, base= 10)
# print(log_level_eiratom)

NORM_neu_eiratom = colors.LogNorm(np.nanmin(Eiratom), np.nanmax(Eiratom))


Nodes=np.fromfile('{}/fort.33'.format(simudir),sep=' ') #Alternatively use fort.33
NN=int(Nodes[0])
XNodes=Nodes[1:NN+1]
YNodes=Nodes[NN+1:]


numberlist = np.zeros(NN)
for i in range(NN):
    numberlist[i] = i

plt.figure(figsize=(7,7))
plt.scatter(XNodes[:500], YNodes[:500])



Triangles = np.loadtxt('{}/fort.34'.format(simudir), 
    skiprows=1, usecols=(1,2,3)) #Alternatively use fort.34
# print(Triangles -1)

TP = tri.Triangulation(XNodes, YNodes, triangles= (Triangles -1))


plt.figure(figsize=(6,12))
plt.tripcolor(TP, Eiratom, shading='flat', cmap= CMAP, norm= NORM_neu_eiratom)
plt.title('Neutral density contour plot')
# plt.title('Neutral density contour plot outerleg')

SM_neu_eiratom= cm.ScalarMappable(NORM_neu_eiratom, CMAP)    
plt.colorbar(SM_neu_eiratom)


"""      
        
        








    

