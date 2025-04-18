# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 18:55:31 2024

@author: ychuang
"""


import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from numpy import ma
from matplotlib.offsetbox import AnchoredText
from matplotlib import colors, cm
from twscan_module.twinscan_prepare import twscan_assist
from contour_plot.contourplot_toolbox import contour_plot_method_collect




class Eirene_contour:


    
    def __init__(self, DF, data, twa: twscan_assist, cpmc: contour_plot_method_collect):
        
        self.DF = DF
        self.data = data
        self.twa = twa
        self.cpmc = cpmc


    
    def twscan_eirene_contourplot_method(self, iterlist, cl_dic, scan_style, plot_option, format_option, norm_type):
        
        
        
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
            
            
            simu_dir = self.data['dirdata']['simudir'][aa[0]][aa[1]]
            
            
            if plot_option == 'Neuden':
                
                dat = self.data['ft46'][aa[0]][aa[1]]['pdena'][:, 0]
                datname = 'Atomic neutral density'
            
            
            elif plot_option == 'testparticleden':

                dat = self.data['ft46'][aa[0]][aa[1]]['pdeni'][:, 0]
                datname = 'test particle density'
            
            
            elif plot_option == 'moleculeden':

                dat = self.data['ft46'][aa[0]][aa[1]]['pdenm'][:, 0]
                datname = 'Molecule neutral density'
            
            elif plot_option == 'ndvxden':
                
                dat = self.data['ft46'][aa[0]][aa[1]]['vxdena'][:, 0]
                datname = 'Atomic neutral x momentum density'
            
            elif plot_option == 'atomic energyden':

                dat = self.data['ft46'][aa[0]][aa[1]]['edena'][:, 0]
                datname = 'Atomic neutral energy density'
            
            
            
            
          
            if self.DF.series_flag == 'twin_scan':
                
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

            if scan_style == 'denscan':
                
                title_ap = float(ap)*pow(10, 5)
                
                if format_option == '1x1':
                    
                    
                    self.cpmc.eirene_contourplot_method(simu_direc = simu_dir, data = dat, 
                                                   itername = aa, axs = axs, datname = datname, norm_type = norm_type)
                   
            elif scan_style == 'tempscan':
                
                title_ap = float(ap)*pow(10, 20)
                
                if format_option == '1x1':
                    
                    self.cpmc.eirene_contourplot_method(simu_direc = simu_dir, data = dat, 
                                                   itername = aa, axs = axs, datname = datname, norm_type = norm_type)
          
            else:
                print('neudenplot_method, please check the scan parameter')
        
            plt.tight_layout()
            plt.xlabel('R: [cm]')
            plt.ylabel('Z: [cm]')



    def twscan_eirene_contourplot(self, scan_style, plot_option, format_option, norm_type):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
               
        if withshift == False and withseries == True:
            
            # series_flag = self.DefaultSettings['series_flag']
            
            
            if self.DF.series_flag == 'twin_scan':
                
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
                
                for x in dircomp[key_a]:
                    keylist_a.append('{:.3f}'.format(x))
                
                for ta in keylist_a:
                    
                    keylist_b = []
                    
                    for x in dircomp[key_b]:
                        keylist_b.append('{:.3f}'.format(x))
                    
                    
                    iter_key, color_dic= self.twa.twscan_plot_prep(ta = ta, 
                    keylist_b = keylist_b, scan_style = scan_style)
                    
                    
                    print('check:')
                    print(iter_key)
                    print(color_dic)
                    # print(label_dic)
                    
                    
                    self.twscan_eirene_contourplot_method(iterlist = iter_key, cl_dic = color_dic, 
                                scan_style = scan_style, norm_type = norm_type,
                    plot_option = plot_option, format_option = format_option)
                    
                    
                    
             
            else:
                print('neteTS_plot, please check the series flag')
    

    





"""



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
        
        








    

