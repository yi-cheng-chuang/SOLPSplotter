# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:50:05 2024

@author: user
"""

"""
The following functions (load_transcoefile_method,  comes from B2TransportParser code written by 
Richard Reksoatmodjo and Jameson Crouse
"""
import re
import numpy as np
import matplotlib.pyplot as plt
import load_coord_method as lcm

def load_transcoefile_method(file ='b2.transport.inputfile', plot=False):
    
    Coefficients = {'1':'Particle density-driven diffusivity',
                    '2': 'Particle pressure-driven diffusivity',
                    '3': 'Ion thermal anomalous diffusivity',
                    '4': 'Electron thermal anomalous diffusivity',
                    '5': 'Poloidal-component of the anomalous ”pinch” velocity',
                    '6': 'Radial-component of the anomalous ”pinch” velocity',
                    '7': 'Anomalous viscosity',
                    '8': 'Anomalous radial electrical conductivity',
                    '9': 'Anomalous radial thermo-electric coefficient'}    
    
    with open(file) as f:
        dataList=f.readlines()
    
    Points={}
    ii=1

    while ii<len(dataList)-2: 
        ndata = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", dataList[ii])
        CoeffID = ndata[1]
        PtNo = int(ndata[3])
        XList = []
        YList = []
        for mm in range(PtNo):
            XList.append(float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",dataList[mm+ii+1])[4]))
            YList.append(float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",dataList[mm+ii+1])[9]))
        Points[CoeffID] = np.array([XList,YList])
        ii=ii+ PtNo + 1
        
    if plot:
        dd=len(Points.keys())
        fig1,ax1=plt.subplots(nrows= dd,ncols=1,sharex=True)
        for ii, jj in enumerate(Points.keys()):
            ax1[ii].plot(Points[jj][0],Points[jj][1])
            ax1[ii].set_ylabel(r'{} $[m^2/s]$'.format(Coefficients[jj]))
        
        ax1[0].set_title(file)    
        ax1[-1].set_xlabel(r'$R-R_{sep}$')
        
    return Points

def Generate_transcoefile_method(trans_pts, CoeffID = 1, SpeciesID = 1, M=[1]):
    '''
    Function that is used to turn the radial points into a readable
    b2.transport.inputfile

    Parameters
    ----------
    trans_pts : nx2 array, x coordinates being r-r_sep and
    y coordinates the coefficient value at that point

    CoeffID : int, integer specifier of transport coefficient 
    type according to SOLPS manual (1 by default)
    
    SpeciesID : int, integer index of transport species (1 by default)
    
    M : float or list, factor to multiply all transport coefficient values by.
    If list, creates a separate multiplied string block for each listed factor

    Returns a formatted string block for use in the b2.transport.inputfile
    -------
    '''
    if type(M) is not list:
        M=[M]
        
    n = len(trans_pts)
    m = 0
    i = CoeffID
    j = SpeciesID
    r = trans_pts
    inputfile={}
    
    for MM in M:
        inputfile[MM] = ' ndata(1, {0}, {1})= {2},\n'.format(i,j,n)
        for m in range(n):
            inputfile[MM] = inputfile[MM] + ' tdata(1, {0}, {1}, {2})= {3}, tdata(2, {0}, {1}, {2})= {4},\n'.format(m+1,i,j,np.round(r[m][0],5),np.round(r[m][1]*MM,5))
            
    return inputfile

def Write_transcoefile_method(file='b2.transport.inputfile', points={},M_1 = True, M=[1]):
    inputfile={}
    if points:
        for k in points.keys():
            inputfile[k] = Generate_transcoefile_method(points[k].T,CoeffID=int(k),M=M)
    else:
        points = load_transcoefile_method(file)
        for k in points.keys():
            inputfile[k] = Generate_transcoefile_method(points[k].T,CoeffID=int(k),M=M)
    if M_1 == True:
        for MM in M:
            with open('{}'.format(file),'w') as f:
                f.write(' &TRANSPORT\n')
                for k in inputfile.keys():
                    f.writelines(inputfile[k][MM])
                f.write(' no_pflux=.false.\n /\n')
    else:
        for MM in M:
            with open('{}.f{}'.format(file,MM),'w') as f:
                f.write(' &TRANSPORT\n')
                for k in inputfile.keys():
                    f.writelines(inputfile[k][MM])
                f.write(' no_pflux=.false.\n /\n')



"""
The following functions (read_transport_files)  comes from SOLPSutils code written by 
Robert Wilcox and Jeremy Lore
"""


def read_transport_files(fileloc, dsa=None, geo=None, state=None, force_read_inpufile=False):
    # Attempts to get transport coefficient data, including pinch term.
    # Fills using b2.transport.parameters, then b2.transport.inputfile if
    # b2mn.dat indicates inputfile is active. 
    #
    # Uses f90nml
    #
    # Inputs:
    # fileloc : location of b2mn.dat, b2.transport.*
    # dsa : dsa from read_dsa (or whatever dx_sep to interpolate onto)
    # geo : geometry dict from read_b2fgmtry 
    # state : state dict from read_b2fstate
    #
    # Note, geo and state only used to get dimensions, only really need
    # geo['ny'], state['ns']
    
    if dsa is None:
        dsa = lcm.read_dsa('dsa')
    
    if geo is None:
        geo = lcm.read_b2fgmtry('../baserun/b2fgmtry')

    if state is None:
        state = lcm.read_b2fstate('b2fstate')
        
    read_inputfile = False
    if force_read_inpufile:
        read_inputfile = True

    # Read b2mn.dat unless overridden by input argument
    if not read_inputfile:
        b2mn = lcm.scrape_b2mn("b2mn.dat")
        if ("b2tqna_inputfile" in b2mn.keys()):
            if b2mn["b2tqna_inputfile"] == 1:
                read_inputfile = True

    # Initialize variables we want defined. The final arrays will have
    # a dimension [ny+2] but b2.transport.inputfile may not!
    dn = np.zeros((geo['ny']+2,state['ns']))
    dp = np.zeros((geo['ny']+2,state['ns']))
    chii = np.zeros((geo['ny']+2,state['ns']))
    chie = np.zeros(geo['ny']+2)
    vlax = np.zeros((geo['ny']+2,state['ns']))
    vlay = np.zeros((geo['ny']+2,state['ns']))
    vsa = np.zeros((geo['ny']+2,state['ns']))
    sig = np.zeros(geo['ny']+2)
    alf = np.zeros(geo['ny']+2)    

    try:
        import f90nml
        parser = f90nml.Parser()
        parser.global_start_index=1
    except:
        print('f90nml required to read transport input files!')
        return None
    
    try:
        nml = parser.read('b2.transport.parameters')
        this = nml['transport']['parm_dna']
        for ispec in range(len(this)):
            dn[:,ispec] = this[ispec]
    except:
        pass
    
    if read_inputfile:
        try:
            nml = parser.read('b2.transport.inputfile')
            tdata = nml['transport']['tdata']
            nS = len(tdata)
            for ispec in range(nS):
                nKinds = len(tdata[ispec])
                for jkind in range(nKinds):
                    this = tdata[ispec][jkind]

                    # Check if this kind was filled with none by f90nml (not defined)
                    test = [i[1] for i in this]
                    if all(val is None for val in test):
                        continue
                    
                    # Read and interpolate back on dsa
                    xRaw = [i[0] for i in this] 
                    yRaw = [i[1] for i in this]
                    yInterp = np.interp(dsa,xRaw,yRaw)
                    
                    if jkind+1 == 1:
                        dn[:,ispec] = yInterp
                    elif jkind+1 == 2:
                        dp[:,ispec] = yInterp
                    elif jkind+1 == 3:
                        chii[:,ispec] = yInterp
                    elif jkind+1 == 4:
                        chie[:] = yInterp
                    elif jkind+1 == 5:
                        vlax[:,ispec] = yInterp                        
                    elif jkind+1 == 6:
                        vlay[:,ispec] = yInterp                        
                    elif jkind+1 == 7:
                        vsa[:,ispec] = yInterp                        
                    elif jkind+1 == 8:
                        sig[:] = yInterp                        
                    elif jkind+1 == 9:
                        alf[:] = yInterp                        
                    
        except:
            pass
        
        
    return dict(dn=dn, dp=dp, chii=chii, chie=chie, vlax=vlax, vlay=vlay, vsa=vsa, sig=sig, alf=alf)




