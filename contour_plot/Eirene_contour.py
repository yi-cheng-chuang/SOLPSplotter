# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 11:20:13 2024

@author: ychuang
"""


from R_diff_calc import Diff_R_calc
from SOLPSplotter_contour import PlotContour
from matplotlib.offsetbox import AnchoredText
import load_B2_data_method as lBdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker
import matplotlib.tri as tri
from matplotlib.colors import LogNorm
from numpy import ma


class eirene_contour(PlotContour):
    def __init__(self, DefaultSettings, loadDS):
        PlotContour.__init__(self, DefaultSettings, loadDS)
        
    
    
    def plot_eireneoutput(self):
            
        
        if self.withshift == True and self.withseries == False:
            
            fig, axs = plt.subplots(2, 2, sharey= True)
            triangle_dic = {}
            
            for ii, aa in enumerate(self.data['dircomp']['multi_shift']):
                
                Attempt = self.data['dircomp']['Attempt'][aa]
                DRT = self.data['dirdata']['outputdir'][aa]['Output']
                Eirout = self.data['dirdata']['outputdir'][aa]['EirOutput']
                simudir = self.data['dirdata']['simudir'][aa]
                simutop = self.data['dirdata']['simutop'][aa]
                XDIM = self.data['b2fgeo'][aa]['nx'] + 2
                YDIM = self.data['b2fgeo'][aa]['ny'] + 2
                
                
                RadLoc = self.data['grid']['RadLoc'][aa]
                VertLoc = self.data['grid']['VertLoc'][aa]
                
                
                tz=np.loadtxt('{}/TriangVertLoc{}'.format(Eirout, str(Attempt)), 
                              usecols = (2))
                              
                tr=np.loadtxt('{}/TriangRadLoc{}'.format(Eirout, str(Attempt)), 
                              usecols = (2))
                
                Eiratom =np.loadtxt('{}/EirAtom{}'.format(Eirout, str(Attempt)), 
                              usecols = (2))

                VVFILE = np.loadtxt('{}/baserun/vvfile.ogr'.format(simutop))

                # Nodes=np.fromfile('{}/{}.tria.{}.nodes'.format(BASEDRT,Device,MeshID),sep=' ') #Alternatively use fort.33
                Nodes=np.fromfile('{}/fort.33'.format(simudir),sep=' ') #Alternatively use fort.33
                NN=int(Nodes[0])
                XNodes=Nodes[1:NN+1]/100
                YNodes=Nodes[NN+1:]/100

                Triangles = np.loadtxt('{}/fort.34'.format(simudir),skiprows=1, usecols=(1,2,3)) #Alternatively use fort.34
                
                
                triangle_dic[aa] = {'tz': tz, 'tr': tr, 'Eiratom': Eiratom, 
                                    'Nodes': Nodes, 'Triangles': Triangles}
                
                
                
                TP = tri.Triangulation(XNodes,YNodes,triangles=(Triangles-1))
                
                triangle_dic[aa] = {'tz': tz, 'tr': tr, 'Eiratom': Eiratom, 
                                    'Nodes': Nodes, 'Triangles': Triangles}
                
                CPB = cm.viridis
                Lnorm = LogNorm(vmax = pow(10, 19), vmin = pow(10, 14))
                
                
                
                if ii <2:
                    
                    axs[ii, 0].plot(VVFILE[:,0]/1000,VVFILE[:,1]/1000,'k-')
                    
                    
                    # smap = cm.ScalarMappable(Lnorm, CPB)
                    
                    
                    IM= axs[ii, 0].tripcolor(TP , Eiratom, norm = Lnorm, cmap = CPB)
                
                else:
                    ik = ii % 2
                    
                    axs[ik, 1].plot(VVFILE[:,0]/1000,VVFILE[:,1]/1000,'k-')
                    
                    # CPB = cm.viridis
                    # Lnorm = LogNorm(vmax = Eiratom.max(), vmin = pow(10, 14))
                    # smap = cm.ScalarMappable(Lnorm, CPB)
                    
                    IM= axs[ik, 1].tripcolor(TP , Eiratom, norm = Lnorm, cmap = CPB)
                    
                
                
                
                if aa == 'org':
                        
                    axs[ii, 0].set_xlim(0, 2)
                    
                    
                elif aa == 'dot3':
                    
                    axs[ii, 0].set_xlim(0.3, 2.3)
                
                elif aa == 'dot5':
                    
                    ik = ii % 2
                    
                    axs[ik, 1].set_xlim(0.5, 2.5)
                
                elif aa == 'dot7':
                    
                    ik = ii % 2
                    
                    axs[ik, 1].set_xlim(0.7, 2.7)
                
                    
                # axs[ii, 0].add_artist(text_list[ii])
                if ii <2:
                    
                    axs[ii, 0].set_xlabel('R [m]')              
                    axs[ii, 0].set_ylim(-2, -0.7)
                    
                else:
                    
                    ik = ii % 2
                    
                    axs[ik, 1].set_xlabel('R [m]')              
                    axs[ik, 1].set_ylim(-2, -0.8)
                
                
                self.data['tridata'] = triangle_dic
                
                fig.supylabel('Z [m]')
                fig.subplots_adjust(right=0.8)
                cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
                plt.suptitle('Atomic neutral density contour plot')
                fig.colorbar(IM, cax= cbar_ax)
        
        
        
    def eirene_contour(self):
            
        
        if self.withshift == True and self.withseries == False:
            
            fig, axs = plt.subplots(1, 4, sharey= True)
            triangle_dic = {}
            
            for ii, aa in enumerate(self.data['dircomp']['multi_shift']):
                
                Attempt = self.data['dircomp']['Attempt'][aa]
                DRT = self.data['dirdata']['outputdir'][aa]['Output']
                Eirout = self.data['dirdata']['outputdir'][aa]['EirOutput']
                simudir = self.data['dirdata']['simudir'][aa]
                simutop = self.data['dirdata']['simutop'][aa]
                XDIM = self.data['b2fgeo'][aa]['nx'] + 2
                YDIM = self.data['b2fgeo'][aa]['ny'] + 2
                
                
                RadLoc = self.data['grid']['RadLoc'][aa]
                VertLoc = self.data['grid']['VertLoc'][aa]
                
                
                tz=np.loadtxt('{}/TriangVertLoc{}'.format(Eirout, str(Attempt)), 
                              usecols = (2))
                              
                tr=np.loadtxt('{}/TriangRadLoc{}'.format(Eirout, str(Attempt)), 
                              usecols = (2))
                
                Eiratom =np.loadtxt('{}/EirAtom{}'.format(Eirout, str(Attempt)), 
                              usecols = (2))

                VVFILE = np.loadtxt('{}/baserun/vvfile.ogr'.format(simutop))

                # Nodes=np.fromfile('{}/{}.tria.{}.nodes'.format(BASEDRT,Device,MeshID),sep=' ') #Alternatively use fort.33
                Nodes=np.fromfile('{}/fort.33'.format(simudir),sep=' ') #Alternatively use fort.33
                NN=int(Nodes[0])
                XNodes=Nodes[1:NN+1]/100
                YNodes=Nodes[NN+1:]/100

                Triangles = np.loadtxt('{}/fort.34'.format(simudir),skiprows=1, usecols=(1,2,3)) #Alternatively use fort.34
                
                
                triangle_dic[aa] = {'tz': tz, 'tr': tr, 'Eiratom': Eiratom, 
                                    'Nodes': Nodes, 'Triangles': Triangles}
                
                
                
                TP = tri.Triangulation(XNodes,YNodes,triangles=(Triangles-1))
                
                
                

                    
                axs[ii].plot(VVFILE[:,0]/1000,VVFILE[:,1]/1000,'k-')
                
                CPB = cm.viridis
                Lnorm = LogNorm(vmax = pow(10, 19), vmin = pow(10, 14))
                smap = cm.ScalarMappable(Lnorm, CPB)
                
                # IM = axs[ii, 0].tripcolor(tr, tz, Eiratom, 
                #              cmap = CPB, norm = Lnorm)
                
                IM= axs[ii].tripcolor(TP , Eiratom, norm = Lnorm, cmap = CPB)
                
                    
                # if aa == 'org':
                        
                #     axs[ii].set_xlim(0, 2)
                    
                    
                # elif aa == 'dot3':
                    
                #     axs[ii].set_xlim(0.3, 2.3)
                
                # elif aa == 'dot5':
                    
                #     ik = ii % 2
                    
                #     axs[ii].set_xlim(0.5, 2.5)
                
                # elif aa == 'dot7':
                    
                #     ik = ii % 2
                    
                #     axs[ii].set_xlim(0.7, 2.7)
                
                axs[ii].set_xlabel('R [m]')              
                # axs[ii].set_ylim(-2, -0.7)
                    
               
                
                
                self.data['tridata'] = triangle_dic
                
                fig.supylabel('Z [m]')
                plt.suptitle('Atomic neutral density contour plot')
                fig.subplots_adjust(right=0.8)
                cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
                
                fig.colorbar(smap, cax= cbar_ax)
        
        
        



