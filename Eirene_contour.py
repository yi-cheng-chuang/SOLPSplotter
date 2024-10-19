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


class target_contour(Diff_R_calc, PlotContour):
    def __init__(self, DefaultSettings, loadDS):
        Diff_R_calc.__init__(self, DefaultSettings, loadDS)
        PlotContour.__init__(self, DefaultSettings, loadDS)
        
    
    
    def plot_all_radial(self):
            

        
        Attempt = self.data['dircomp']['Attempt']
        DRT = self.data['dirdata']['outputdir']['Output']
        Eirout = self.data['dirdata']['outputdir']['EirOutput']
        simudir = self.data['dirdata']['simudir']
        XDIM = self.data['b2fgeo']['nx'] + 2
        YDIM = self.data['b2fgeo']['ny'] + 2
        
        RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                    usecols = (3)).reshape((YDIM, XDIM))
        VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                      usecols = (3)).reshape((YDIM,XDIM))
        
        
        tz=np.loadtxt('{}/TriangVertLoc{}'.format(Eirout, str(Attempt)), 
                      usecols = (2))
                      
        tr=np.loadtxt('{}/TriangRadLoc{}'.format(Eirout, str(Attempt)), 
                      usecols = (2))
        
        Eiratom =np.loadtxt('{}/EirAtom{}'.format(Eirout, str(Attempt)), 
                      usecols = (2))



        VVFILE = np.loadtxt('{}/{}/{}/{}/baserun/vvfile.ogr'.format(basedrt, dev, shot, shift))

        # Nodes=np.fromfile('{}/{}.tria.{}.nodes'.format(BASEDRT,Device,MeshID),sep=' ') #Alternatively use fort.33
        Nodes=np.fromfile('{}/fort.33'.format(newbase),sep=' ') #Alternatively use fort.33
        NN=int(Nodes[0])
        XNodes=Nodes[1:NN+1]/100
        YNodes=Nodes[NN+1:]/100

        Triangles=np.loadtxt('{}/fort.34'.format(newbase),skiprows=1, usecols=(1,2,3)) #Alternatively use fort.34

        TP = tri.Triangulation(XNodes,YNodes,triangles=(Triangles-1))
        
        
        Contour.plot(VVFILE[:,0]/1000,VVFILE[:,1]/1000,'k-')
        IM=Contour.tripcolor(TP,Data)
        Contour.set_xlabel('Radial Position R (m)')
        Contour.set_ylabel('Vertical Position Z (m)')
        Contour.set_xlim(0, 2)
        Contour.set_ylim(-2, -0.5)
        plt.colorbar(IM,ax=Contour)



