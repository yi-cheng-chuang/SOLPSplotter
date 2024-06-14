# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 18:37:40 2023

@author: user
"""

from os import path, environ
import numpy as np



"""

The following functions (read_transport_files)  comes from SOLPSutils code written by 
Robert Wilcox and Jeremy Lore

"""


def scrape_b2mn(fname = 'b2mn.dat'):
    b2mn = {}
    if not path.exists(fname):
        print('ERROR: b2mn.dat file not found:',fname)
        return b2mn
    
    with open(fname, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        sline = line.strip()
        if len(line) <= 1:
            continue
        if sline[0] == "#" or sline[0] == "*":
            continue
        else:
            count_quotes = 0
            quote_pos = []
            for i, c in enumerate(sline):
                if c == "'":
                    count_quotes = count_quotes + 1
                    quote_pos.append(i)
            # count_quotes = sline.count("'")
            if count_quotes == 2:
                # For cases where variable is enclosed in single quotes but value is not: 'b2mwti_jxa'
                thisvar = sline[quote_pos[0]+1:quote_pos[1]]
                this = sline[quote_pos[1]+1:-1]
            else:
                # Typically this is for 4 single quotes
                # For cases where both variable and value are enclosed in single quotes: 'b2mwti_jxa'   '36'
                #
                # This can occur when some comment after the value has quotes in it
                # e.g., "'b2sicf_phm0'  '0.0'  Old value '1.0'"
                thisvar = sline[quote_pos[0]+1:quote_pos[1]]
                this = sline[quote_pos[2]+1:quote_pos[3]]


            if thisvar == "b2mwti_jxa":
                b2mn['jxa'] = int(this)
            if thisvar == "b2tqna_inputfile":
                b2mn['b2tqna_inputfile'] = int(this)

    return b2mn



def read_b2fgmtry(fname):
    if not path.exists(fname):
        print('ERROR: b2fgmtry file not found:',fname)
        return None

    data = []
    with open(fname, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):

        # Special handling for first few lines
        if i == 0:
            version = line.split()[0][7:-1]
            geo = {'version':version}
            continue
        elif i == 1:
            continue
        elif i == 2:
            # Assume starts with nx,ny after version
            geo['nx'] = int(line.split()[0])
            geo['ny'] = int(line.split()[1])
            numcells = (geo['nx']+2)*(geo['ny']+2)
            continue

        if line.split()[0] == '*cf:':
            vartype = line.split()[1]
            varsize = int(line.split()[2])
            varname = line.split()[3]
            # Some variables have no entries depending on config
            if varsize == 0:
                geo[varname] = None
            data = []
        else:
            # Parse by type
            if vartype == "char":
                geo[varname] = line.strip()
            else:
                splitline = line.split()
                for value in splitline:
                    if vartype == "int":
                        data.append(int(value))
                    else:
                        data.append(float(value))

                if len(data) == varsize:
                    if varsize%numcells == 0:
                        geo[varname] = np.array(data).reshape([geo['nx']+2,geo['ny']+2,int(varsize/numcells)], order = 'F')
                    else:
                        geo[varname] = np.array(data)

    return geo


 
def loadg(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()

    # read first line for case string and grid size
    line = lines[0]
    words = line.split()

    nw = int(words[-2])
    nh = int(words[-1])
    psi = np.linspace(0, 1, nw)

    # read in scalar parameters
    #   note: word size of 16 characters each is assumed for all of the data to be read

    # line 1
    line = lines[1]
    rdim = float(line[0:16])
    zdim = float(line[16:32])
    rcentr = float(line[32:48])
    rleft = float(line[48:64])
    zmid = float(line[64:80])

    # line 2
    line = lines[2]
    rmaxis = float(line[0:16])
    zmaxis = float(line[16:32])
    simag = float(line[32:48])
    sibry = float(line[48:64])
    bcentr = float(line[64:80])

    # line 3
    line = lines[3]
    current = float(line[0:16])

    # read in profiles
    #   start by reading entire file into single list then split into individual profiles
    #   first block has 5 profiles of length nw and one array of length nh*nw

    temp = []
    count = 0
    lnum = 5
    terms = 5 * nw + nw * nh
    while count < terms:
        line = lines[lnum]
        numchar = len(line)
        nwords = numchar // 16
        count1 = 0
        while count1 < nwords:
            i1 = count1 * 16
            i2 = i1 + 16
            temp.append(float(line[i1:i2]))
            count1 += 1
            count += 1
        lnum += 1

    fpol = temp[0:nw]
    pres = temp[nw:2 * nw]
    ffprime = temp[2 * nw:3 * nw]
    pprime = temp[3 * nw:4 * nw]
    psirz_temp = temp[4 * nw:(4 + nh) * nw]
    qpsi = temp[(4 + nh) * nw:]
    
    
    
    
    # split psirz up into 2D matrix
    count = 0
    psirz = []
    while count < nh:
        ind1 = count * nw
        ind2 = ind1 + nw
        psirz.append(psirz_temp[ind1:ind2])
        count += 1

    # scalars for length of boundary and limiter arrays
    line = lines[lnum]
    words = line.split()
    nbbbs = int(words[0])
    limitr = int(words[1])

    # read boundary and limiter points into temp array
    
    
    
    temp = []
    count = 0
    terms = 2 * (nbbbs + limitr)
    lnum += 1
    
    
    print('the line index number for boundary is: {}'.format(str(lnum)))
    
    while count < terms:
        line = lines[lnum]
        numchar = len(line)
        nwords = numchar // 16
        count1 = 0
        while count1 < nwords:
            i1 = count1 * 16
            i2 = i1 + 16
            temp.append(float(line[i1:i2]))
            count1 += 1
            count += 1
        lnum += 1
    bdry_temp = temp[0:(2 * nbbbs)]
    limit_temp = temp[(2 * nbbbs):]

    # split boundary into (R,Z) pairs
    count = 0
    rbdry = []
    zbdry = []
    while count < len(bdry_temp) - 1:
        rbdry.append(bdry_temp[count])
        zbdry.append(bdry_temp[count + 1])
        count += 2

    # split limiter into (R,Z) pairs
    count = 0
    rlim = []
    zlim = []
    while count < len(limit_temp) - 1:
        rlim.append(limit_temp[count])
        zlim.append(limit_temp[count + 1])
        count += 2

    g = dict(nw=nw, nh=nh, rdim=rdim, zdim=zdim, rcentr=rcentr, rleft=rleft, zmid=zmid,
             rmaxis=rmaxis, zmaxis=zmaxis, simag=simag, sibry=sibry, current=current,
             fpol=np.array(fpol),
             ffprime=np.array(ffprime), pprime=np.array(pprime), psirz=np.array(psirz),
             qpsi=np.array(qpsi), nbbbs=nbbbs, bcentr=bcentr,
             pres=np.array(pres), limitr=limitr, rbdry=np.array(rbdry),
             zbdry=np.array(zbdry), rlim=np.array(rlim), zlim=np.array(zlim))

    return g


def read_dsa(fname):
    if not path.exists(fname):
        print('ERROR: b2fgmtry file not found:',fname)
        return None

    dsa = []
    with open(fname, 'r') as f:
        lines = f.readlines()

    for i,line in enumerate(lines):
        dsa.append(float(line.split()[0]))

    return dsa
