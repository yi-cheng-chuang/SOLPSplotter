# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 23:00:37 2024

@author: ychuang
"""

#! /usr/bin/env python
import netCDF4
import os
import matplotlib.pyplot as plt
import sys
import numpy
import re
from SOLPSplotter_PB import PB_plot





class EB_plot(PB_plot):
    
    def __init__(self, DefaultSettings, loadDS):
        PB_plot.__init__(self, DefaultSettings, loadDS)


    def plot_EB(self, scan_style, dat_size):
        
        
        if self.withshift == False and self.withseries == True:
            
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
                
                
                for x in dircomp[key_a]:
                    keylist_a.append('{:.3f}'.format(x))
                
                for ta in keylist_a:
                    
                    keylist_b = []
                    
                    for x in dircomp[key_b]:
                        keylist_b.append('{:.3f}'.format(x))
                    
                        
                    iter_key, color_dic, scan_title, label_dic = self.twinscan_prep(ta = ta, 
                    keylist_b = keylist_b, scan_style = scan_style, dat_size = dat_size)
                    
                    # print('check:')
                    # print(iter_key)
                    # print(color_dic)
                    # print(label_dic)
                                      
                    self.plot_PB_method(iterlist = iter_key, scandetail = scan_title,
                            cl_dic = color_dic, A_dic = label_dic,
                            scan_style = scan_style, dat_size = dat_size)
            
    
    
    
    def plot_PB_method(self, iterlist, cl_dic, A_dic, scan_style, 
                          scandetail, dat_size):
        
        
        
        if self.withshift == False and self.withseries == True:
            
            
            if self.series_flag == 'twin_scan':
            
                fig, axs = plt.subplots()
                
                
                nx = self.data['b2fgeo']['nx']
                ny = self.data['b2fgeo']['ny']
                
                
                if dat_size == 'full':
        
                    dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
                
                elif dat_size == 'small':
                    dat_struc = {'size': dat_size, 'nx': nx, 'ny': ny}
                
                
                
                for aa in iterlist:
                    
                    
                    if scan_style == 'tempscan':
                        
                        ad = aa[1]
                        axs.set_title('{} Temperature scan with Ne = {}'.format('Normalised energy error', scandetail))
                    
                    elif scan_style == 'denscan':
                        
                        ad = aa[0]
                        axs.set_title('{} Density scan with Te = {} eV'.format('Normalised energy error', scandetail))
                    
                    else:
                        print('neteTSplot_method, please check scan_style')
                    
                    
                    sim_dir = self.data['dirdata']['simudir']
                    f_loc = sim_dir[aa[0]][aa[1]]

                    
                    f = self.load_b2tallies(file_loc = f_loc)
                    
                    vreg = f.dimensions['vreg'].size
                    xreg = f.dimensions['xreg'].size
                    yreg = f.dimensions['yreg'].size
                    ns = f.dimensions['ns'].size
                    time = f.dimensions['time'].size
                    times = f.variables['times']
                    species_names=f.variables['species']
                    species=[b''.join(species_names[i,:]).strip().decode('utf-8') for i in range(species_names.shape[0])]
                    elements=[re.sub('[-+0-9]','',species[i]) for i in range(len(species))]
                    mask=[re.match('[a-zA-Z]+0$',species[i])!=None for i in range(len(species))]
                    s=0
                    bounds=[]
                    for i in range(mask.count(True)-1):
                        e=mask.index(True,s+1)
                        bounds+=[[s,e]]
                        s=e
                    bounds+=[[s,len(mask)]]


                    print(bounds)


                    "ion index"

                    fhixreg = f.variables['fhixreg']
                    fhiyreg = f.variables['fhiyreg']

                    "electron index"


                    fhexreg = f.variables['fhexreg']
                    fheyreg = f.variables['fheyreg']


                    "volume index"

                    b2divua = f.variables['b2divua']

                    b2exba = f.variables['b2exba']

                    b2stbr_shi_reg = f.variables['b2stbr_shi_reg']

                    b2sext_shi_reg = f.variables['b2sext_shi_reg']


                    b2divue = f.variables['b2divue']

                    b2exbe = f.variables['b2exbe']

                    b2stbr_she_reg = f.variables['b2stbr_she_reg']

                    b2sext_she_reg = f.variables['b2sext_she_reg']


                    rsanareg = f.variables['rsanareg']

                    rranareg = f.variables['rranareg']

                    rsahireg = f.variables['rsahireg']

                    rrahireg = f.variables['rrahireg']

                    rqahereg = f.variables['rqahereg']

                    rcxhireg = f.variables['rcxhireg']

                    b2visa = f.variables['b2visa']

                    b2joule = f.variables['b2joule']

                    b2fraa = f.variables['b2fraa']

                    b2wrong1 = f.variables['b2wrong1']

                    b2wrong2 = f.variables['b2wrong2']

                    b2wrong3 = f.variables['b2wrong3']

                    rdneureg = f.variables['rdneureg']


                    poreg = f.variables['poreg']

                    b2stbr_sna_reg = f.variables['b2stbr_sna_reg']




                    if vreg == 5:
                        FULL_X = numpy.array([0,1,0,0,-1,0,0])
                        FULL_Y = numpy.array([0,1,1,1,0,-1,-1,-1])
                        # pot_V = numpy.array([1, 2, 3, 4, 0]*1.602176634*pow(10, -19))
                    elif vreg ==2:
                        FULL_X = numpy.array([0,1,-1])
                        FULL_Y = numpy.array([0,1,-1])
                    else:
                        raise ValueError('Value of vreg=%s not currently coded' % vreg)


                    for i in range(len(bounds)):
                        
                        
                        fhex = (fhexreg[:,:]*FULL_X).sum(axis=1)
                        fhix = (fhixreg[:,:]*FULL_X).sum(axis=1)

                        fhx_dat = fhexreg[:,:] + fhixreg[:,:]

                        fhx = fhex + fhix

                        # print(numpy.shape(fhexreg[:,:]*FULL_X))


                        fhey = (fheyreg[:,:]*FULL_Y).sum(axis=1)
                        fhiy = (fhiyreg[:,:]*FULL_Y).sum(axis=1)

                        fhy_dat = fheyreg[:,:] + fhiyreg[:,:]

                        fhy = fhex + fhiy

                        
                        b2bnrp = (b2stbr_sna_reg[:, 0, 1:]*1.602176634*pow(10, -19)).sum(axis=1)


                        rsap = (rsanareg[:, 0, 1:].copy()*1.602176634*pow(10, -19)).sum(axis=1)

                        rrap = (rranareg[:, 0, 1:].copy()*1.602176634*pow(10, -19)).sum(axis=1)

                        
                        rdneureg = rdneureg[:,1:].copy()

                        rdneu = rdneureg.sum(axis=1)

                        b2stbr_e = b2stbr_she_reg[:, 2].copy()

                        print(numpy.shape(b2stbr_e))

                        b2stbr_i = b2stbr_shi_reg[:, 2].copy()

                        b2stbr = b2stbr_e + b2stbr_i

                        # print(b2stbr)

                        b2sext_e = b2sext_she_reg[:,1:].copy()

                        b2sext_i = b2sext_shi_reg[:,1:].copy()

                        b2sext = b2sext_e.sum(axis=1) + b2sext_i.sum(axis=1)
                        
                        b2divu_a = b2divua[:,2].copy()
                        
                        b2divu_e = b2divue[:,2].copy()
                        
                        print('the shape of b2divu is:')
                        print(numpy.shape(b2divu_e))
                        
                        b2divu = b2divu_e + b2divu_a
                        

                        b2exb_a = b2exba[:,2].copy()
                        
                        b2exb_e = b2exbe[:,2].copy()
                        
                        b2exb = b2exb_e + b2exb_a
                        
                        b2visa = b2visa[:, 2].copy()

                        b2joule = b2joule[:, 2].copy()

                        b2fraa = b2fraa[:, 2].copy()
                        
                        b2wrong1 = b2wrong1[:, 2].copy()
                        
                        b2wrong2 = b2wrong2[:, 2].copy()
                        
                        b2wrong3 = b2wrong3[:, 2].copy()
                        
                        
                          
                        rsahireg = ((rsahireg[:,bounds[i][0]:bounds[i][1] ,1:].copy()).sum(axis=1)).sum(axis=1)
                        
                        print('the shape of rrahireg is:')
                        print(numpy.shape(rsahireg))
                        
                        rrahireg = ((rrahireg[:,bounds[i][0]:bounds[i][1], 1:].copy()).sum(axis=1)).sum(axis=1)

                        rqahereg = ((rqahereg[:,bounds[i][0]:bounds[i][1], 1:].copy()).sum(axis=1)).sum(axis=1)

                        rcxhireg = ((rcxhireg[:,bounds[i][0]:bounds[i][1], 1:].copy()).sum(axis=1)).sum(axis=1)
                        


                        fh_norm = numpy.max([numpy.max(numpy.abs(fhx_dat[:,:])),numpy.max(numpy.abs(fhy_dat[:,:]))])
                        
                        
                        f_1 = fhx + fhy
                        
                        f_2 = -rsahireg - rrahireg + rqahereg -rcxhireg
                        
                        f_3 = -b2divu - b2exb - b2visa - b2joule - b2fraa - b2stbr
                        
                        f_4 = - b2wrong1 -b2wrong2 -b2wrong3
                        
                        f_5 = rsap - rrap + b2bnrp
                        
                        
                        EB_formula = f_1 + f_2 + f_3 + f_4 + f_5
                        
                        
                        print('total', fh_norm, (EB_formula).mean(), ((EB_formula)/fh_norm).mean())
                        plt.plot(times[:],((EB_formula)/fh_norm),   
                        label= '{}'.format(A_dic[ad]), ls = '-', color = cl_dic[ad])
                    
            
                axs.set_xlabel("time (s)")
                axs.set_ylabel('normalised energy error')
                axs.legend(loc="best")
            

