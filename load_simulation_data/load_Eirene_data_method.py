# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 20:08:23 2024

@author: user
"""

import numpy as np
import os.path


def is_number(s):
    ''' checks to see if s is a number, useful for parsing outputfiles'''
    try:
        float(s)
        return True
    except ValueError:
        return False



def read_ft44_field(fid,ver,fieldname,dims,intField=False):
    '''Auxiliary routine to read fields from fort.44 file
    fid is the file object i.e fid = open(fileLocation)
    Verion 20160829: field label and size are specified in fort.44
    fieldname is the name of the variable to find
    dims is the dimension of that variable the array will be shaped into
    intField says whether to return the values as an integer or float'''
#    Do consistency check on the data
    if (ver >= 20160829):
        # Search the file until identifier 'fieldname' is found
        line = fid.readline().rstrip()
        while fieldname not in line:
            line = fid.readline().rstrip()
            if len(line) == 0: print('read_ft44_field: EOF reached without finding '+str(fieldname))
        # Consistency check: number of elements specified in the file should equal
        # prod(dims)
        for i in range(len(line.split())):
            if is_number(line.split()[i]): numin = int(line.split()[i])
	
        if (numin != np.prod(dims) and 'wld' not in fieldname):
            print('issue with field '+fieldname)
            print("numin="+str(numin))
            print("np.prod(dims)="+str(np.prod(dims)))
            print('read_ft44_rfield: inconsistent number of input elements.')
            print('if this is a wall paramter, could be fine, check it.')
            print('number of walls is hardcoded in, need to fix')
        elif (numin!= np.prod(dims) and ('wldnek' in fieldname or 'wldnep' in fieldname)):
            dims = [numin]
        elif (numin!= np.prod(dims) and ('wldna' in fieldname or 'ewlda' in fieldname or 'wldnm' in fieldname or 'ewldm' in fieldname)):
            dims[1] = int((numin)/dims[0])
        print(numin)
            

    # Read the data
    fieldVal=[]
    # collect field values
    while (len(fieldVal) != numin):
        line = fid.readline().rstrip()
        if ('wld' in fieldname) and len(fieldVal)>=numin-1: break
        for i in range(len(line.split())):
            if ('wlpump' in fieldname):
                if not is_number(line.split()[i]): continue
            if (intField): fieldVal.append(int(line.split()[i]))
            else: fieldVal.append(float(line.split()[i]))
    fieldVal=np.array(fieldVal)
    
    if (np.size(dims) > 1 and 'wld' not in fieldname): fieldVal = fieldVal.reshape(dims,order='F').copy()
    print(dims)
    if (np.size(dims) > 1 and 'wld' in fieldname): fieldVal = fieldVal.reshape(dims).copy() 

    return fieldVal



def read_ft44(fileName):
    '''Read fort.44 file
    For now
    - only fort.44 version 20081111 recognized
    - assuming nfla = 1 until a better fix
    - assuming nlwrmsh = 1 until a better fix

    note March 2023: I don't remember what nfla and nlwrmsh do, but I think it may be irrelevant at this point'''

    print('read_ft44: assuming nlwrmsh = 1, nfla = 1.')
    nlwrmsh = 1
    nfla = 1

    fid = open(fileName)
    if (fid == -1): print("read_ft44: can't open file")

# Read dimensions

# nx, ny, version
    dims = fid.readline().rstrip().split()
    nx = int(dims[0])
    ny = int(dims[1])
    ver = int(dims[2])

    if (ver != 20081111 and ver != 20160829 and ver != 20170328):
        print('read_ft44: unknown format of fort.44 file (this is usually fine)')

# go to new line (skip reading a possible git-hash)
#     fid.readline().rstrip()

# natm, nmol, nion
    dims = fid.readline().rstrip().split()
    natm = int(dims[0])
    nmol = int(dims[1])
    nion = int(dims[2])
# for now, ignore reading species labels

    for i in range(natm): line = fid.readline().rstrip()
    for i in range(nmol): line = fid.readline().rstrip()
    for i in range(nion): line = fid.readline().rstrip()

# Read basic data, there is more, I might grab it if I find out I need it
    class ft44Results:
        def __init__(self):
            self.dab2     = read_ft44_field(fid,ver,'dab2',[nx,ny,natm]);
            self.tab2     = read_ft44_field(fid,ver,'tab2',[nx,ny,natm]);
            self.dmb2     = read_ft44_field(fid,ver,'dmb2',[nx,ny,nmol]);
            self.tmb2     = read_ft44_field(fid,ver,'tmb2',[nx,ny,nmol]);
            self.dib2     = read_ft44_field(fid,ver,'dib2',[nx,ny,nion]);
            self.tib2     = read_ft44_field(fid,ver,'tib2',[nx,ny,nion]);
            self.rfluxa   = read_ft44_field(fid,ver,'rfluxa',[nx,ny,natm]);
            self.rfluxm   = read_ft44_field(fid,ver,'rfluxm',[nx,ny,nmol]);
            self.pfluxa   = read_ft44_field(fid,ver,'pfluxa',[nx,ny,natm]);
            self.pfluxm   = read_ft44_field(fid,ver,'pfluxm',[nx,ny,nmol]);
            self.refluxa  = read_ft44_field(fid,ver,'refluxa',[nx,ny,natm]);
            self.refluxm  = read_ft44_field(fid,ver,'refluxm',[nx,ny,nmol]);
            self.pefluxa  = read_ft44_field(fid,ver,'pefluxa',[nx,ny,natm]);
            self.pefluxm  = read_ft44_field(fid,ver,'pefluxm',[nx,ny,nmol]);
            self.emiss    = read_ft44_field(fid,ver,'emiss',[nx,ny,1]);
            self.emissmol = read_ft44_field(fid,ver,'emissmol',[nx,ny,1]);
            self.srcml    = read_ft44_field(fid,ver,'srcml',[nx,ny,nmol]);
            self.edissml  = read_ft44_field(fid,ver,'edissml',[nx,ny,nmol]);
            [wallN,someNumber,istra] = fid.readline().split()
            wldna = []
            ewlda = []
            wldnm = []
            ewldm = []
            wldnek = []
            wldnep = []
            wldra = []
            wldrm = []
            wldnek.append(read_ft44_field(fid,ver,'wldnek(0)',[-1]))
            wldnep.append(read_ft44_field(fid,ver,'wldnep(0)',[-1]))
            wldna.append(read_ft44_field(fid,ver,'wldna(0)',[natm,-1]))
            ewlda.append(read_ft44_field(fid,ver,'ewlda(0)',[natm,-1]))
            wldnm.append(read_ft44_field(fid,ver,'wldnm(0)',[nmol,-1]))
            ewldm.append(read_ft44_field(fid,ver,'ewldm(0)',[nmol,-1]))
            self.wall_geometry = read_ft44_field(fid,ver,'wall_geometry',[620])
            wldra.append(read_ft44_field(fid,ver,'wldra(0)',[natm,-1]))
            wldrm.append(read_ft44_field(fid,ver,'wldrm(0)',[nmol,-1]))
            for i in range(1,int(istra)+1):
                if i<10:
                    wldnek.append(read_ft44_field(fid,ver,'wldnek(  '+str(i)+')',[-1]))
                    wldnep.append(read_ft44_field(fid,ver,'wldnep(  '+str(i)+')',[-1]))
                    wldna.append(read_ft44_field(fid,ver,'wldna(  '+str(i)+')',[natm,-1]))
                    ewlda.append(read_ft44_field(fid,ver,'ewlda(  '+str(i)+')',[natm,-1]))
                    wldnm.append(read_ft44_field(fid,ver,'wldnm(  '+str(i)+')',[nmol,-1]))
                    ewldm.append(read_ft44_field(fid,ver,'ewldm(  '+str(i)+')',[nmol,-1]))
                    wldra.append(read_ft44_field(fid,ver,'wldra(  '+str(i)+')',[natm,-1]))
                    wldrm.append(read_ft44_field(fid,ver,'wldrm(  '+str(i)+')',[nmol,-1]))
                else:
                    wldnek.append(read_ft44_field(fid,ver,'wldnek( '+str(i)+')',[-1]))
                    wldnep.append(read_ft44_field(fid,ver,'wldnep( '+str(i)+')',[-1]))
                    wldna.append(read_ft44_field(fid,ver,'wldna( '+str(i)+')',[natm,-1]))
                    ewlda.append(read_ft44_field(fid,ver,'ewlda( '+str(i)+')',[natm,-1]))
                    wldnm.append(read_ft44_field(fid,ver,'wldnm( '+str(i)+')',[nmol,-1]))
                    ewldm.append(read_ft44_field(fid,ver,'ewldm( '+str(i)+')',[nmol,-1]))
                    wldra.append(read_ft44_field(fid,ver,'wldra( '+str(i)+')',[natm,-1]))
                    wldrm.append(read_ft44_field(fid,ver,'wldrm( '+str(i)+')',[nmol,-1]))
            self.wldna    = np.array(wldna)
            self.ewlda    = np.array(ewlda)
            self.wldnm    = np.array(wldnm)
            self.ewldm    = np.array(ewldm)
            self.wldnek   = np.array(wldnek)
            self.wldnep   = np.array(wldnep)
            self.wldra    = np.array(wldra)
            self.wldrm    = np.array(wldrm)
            wldpeb = []
            wldpeb.append(read_ft44_field(fid,ver,'wldpeb(0)',[-1]))
            for i in range(1,int(istra)+1):
                if i<10:
                    wldpeb.append(read_ft44_field(fid,ver,'wldpeb(  '+str(i)+')',[-1]))
                else:
                    wldpeb.append(read_ft44_field(fid,ver,'wldpeb( '+str(i)+')',[-1]))
            self.wldpeb   = np.array(wldpeb)
            self.wlarea   = read_ft44_field(fid,ver,'wlarea',[168])
#            self.wldspt0   = read_ft44_field(fid,ver,'wldspt(0)',[natm,-1])
#            self.wldspt1   = read_ft44_field(fid,ver,'wldspt(  1)',[natm,-1])
#            self.wldspt2   = read_ft44_field(fid,ver,'wldspt(  2)',[natm,-1])
#            self.wldspt3   = read_ft44_field(fid,ver,'wldspt(  3)',[natm,-1])
#            self.wldspt4   = read_ft44_field(fid,ver,'wldspt(  4)',[natm,-1])
#            self.wldspt5   = read_ft44_field(fid,ver,'wldspt(  5)',[natm,-1])
#            self.wldspt6   = read_ft44_field(fid,ver,'wldspt(  6)',[natm,-1])
 #           self.wlpumpA   = read_ft44_field(fid,ver,'wlpump(A)',[natm,88])
 #           self.wlpumpM   = read_ft44_field(fid,ver,'wlpump(M)',[nmol,88])
            self.eneutrad  = read_ft44_field(fid,ver,'eneutrad',[nx,ny,natm])
            self.emolrad  = read_ft44_field(fid,ver,'emolrad',[nx,ny,nmol])
            self.eionrad  = read_ft44_field(fid,ver,'eionrad',[nx,ny,nion]) 
    ft44 = ft44Results()
    fid.close()
    print('done reading ft44 file')
    return ft44



def read_ft46(fileName):
    '''reads fort.46 file and tries to convert to SI units,
    For now, only fort.46 version 20160513 recognized, though later versions probably work'''

    fid = open(fileName)
    if (fid == -1): print("read_ft44: can't open file")

    # Read dimensions

    # ntri, version, avoid reading git-hash
    line = fid.readline().rstrip().split()
    ntri = int(line[0])
    ver  = int(line[1])

    if ver != 20160513 and ver != 20160829 and ver != 20170930:
        print('read_ft46: unknown format of fort.46 file')

    # natm, nmol, nion
    dims = fid.readline().rstrip().split()
    natm = int(dims[0])
    nmol = int(dims[1])
    nion = int(dims[2])

    # for now, ignore reading species labels
    for i in range(natm): fid.readline().rstrip()
    for i in range(nmol): fid.readline().rstrip()
    for i in range(nion): fid.readline().rstrip()

    eV   = 1.6021765650000000E-019
    # Read data
    class ft46Results:
        def __init__(self):
            self.pdena  = read_ft44_field(fid,ver,'pdena',[ntri,natm])*1e6# m^{-3}
            self.pdenm  = read_ft44_field(fid,ver,'pdenm',[ntri,nmol])*1e6
            self.pdeni  = read_ft44_field(fid,ver,'pdeni',[ntri,nion])*1e6
            
            self.edena  = read_ft44_field(fid,ver,'edena',[ntri,natm])*1e6*eV# J m^{-3}
            self.edenm  = read_ft44_field(fid,ver,'edenm',[ntri,nmol])*1e6*eV
            self.edeni  = read_ft44_field(fid,ver,'edeni',[ntri,nion])*1e6*eV

            self.vxdena = read_ft44_field(fid,ver,'vxdena',[ntri,natm])*1e1# kg s^{-1} m^{-2}
            self.vxdenm = read_ft44_field(fid,ver,'vxdenm',[ntri,nmol])*1e1
            self.vxdeni = read_ft44_field(fid,ver,'vxdeni',[ntri,nion])*1e1

            self.vydena = read_ft44_field(fid,ver,'vydena',[ntri,natm])*1e1# kg s^{-1} m^{-2}
            self.vydenm = read_ft44_field(fid,ver,'vydenm',[ntri,nmol])*1e1
            self.vydeni = read_ft44_field(fid,ver,'vydeni',[ntri,nion])*1e1

            self.vzdena = read_ft44_field(fid,ver,'vzdena',[ntri,natm])*1e1# kg s^{-1} m^{-2}
            self.vzdenm = read_ft44_field(fid,ver,'vzdenm',[ntri,nmol])*1e1
            self.vzdeni = read_ft44_field(fid,ver,'vzdeni',[ntri,nion])*1e1
            self.vol =    read_ft44_field(fid,ver,'volumes',[ntri])
            self.pux    = read_ft44_field(fid,ver,'pux',[ntri])
            self.puy    = read_ft44_field(fid,ver,'puy',[ntri])
            
    ft46 = ft46Results()
    # Close file
    fid.close()
    print('done reading ft46 file')
    return ft46





"""

Save for old version:




def read_ft44_field(fid,ver,fieldname,dims,intField=False):
    '''Auxiliary routine to read fields from fort.44 file
    fid is the file object i.e fid = open(fileLocation)
    Verion 20160829: field label and size are specified in fort.44
    fieldname is the name of the variable to find
    dims is the dimension of that variable the array will be shaped into
    intField says whether to return the values as an integer or float'''
#    Do consistency check on the data
    if (ver >= 20160829):
        # Search the file until identifier 'fieldname' is found
        line = fid.readline().rstrip()
        while fieldname not in line:
            line = fid.readline().rstrip()
            if len(line) == 0: print('read_ft44_field: EOF reached without finding '+str(fieldname))
        # Consistency check: number of elements specified in the file should equal
        # prod(dims)
        for i in range(len(line.split())):
            if is_number(line.split()[i]): numin = int(line.split()[i])
	
        if (numin != np.prod(dims) and 'wld' not in fieldname):
            print('issue with field '+fieldname)
            print("numin="+str(numin))
            print("np.prod(dims)="+str(np.prod(dims)))
            print('read_ft44_rfield: inconsistent number of input elements.')
            print('if this is a wall paramter, could be fine, check it.')
            print('number of walls is hardcoded in, need to fix')
        elif (numin!= np.prod(dims) and ('wldnek' in fieldname or 'wldnep' in fieldname)):
            dims = [numin]
        elif (numin!= np.prod(dims) and ('wldna' in fieldname or 'ewlda' in fieldname or 'wldnm' in fieldname or 'ewldm' in fieldname)):
            dims[1] = int(numin/dims[0])
            

    # Read the data
    fieldVal=[]
    # collect field values
    while (len(fieldVal) != numin):
        line = fid.readline().rstrip()
        if ('wld' in fieldname) and len(fieldVal)>=numin-1: break
        for i in range(len(line.split())):
            if ('wlpump' in fieldname):
                if not is_number(line.split()[i]): continue
            if (intField): fieldVal.append(int(line.split()[i]))
            else: fieldVal.append(float(line.split()[i]))
    fieldVal=np.array(fieldVal)
    if (np.size(dims) > 1 and 'wld' not in fieldname): fieldVal = fieldVal.reshape(dims,order='F').copy()
    if (np.size(dims) > 1 and 'wld' in fieldname): fieldVal = fieldVal.reshape(dims).copy() 

    return fieldVal


"""




