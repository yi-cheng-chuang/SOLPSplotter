# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 11:20:13 2024

@author: ychuang
"""

 def plot_all_radial(self):
         
     'contour plot for ne & te'
     
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