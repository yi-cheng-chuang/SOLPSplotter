# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 18:55:31 2024

@author: ychuang
"""













def eirene_contour_plot(self):
    
       
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
    

