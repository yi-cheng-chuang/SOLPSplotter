# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 16:40:01 2024

@author: user
"""

"""
The following functions (is_number, read_field, read_b2fstate, read_b2fplasmf)
comes from tutorialFunctionDoc code written by Eric Emdee

"""



import numpy as np
# import os.path
import os
# from os import path, environ


def is_number(s):
    ''' checks to see if s is a number, useful for parsing outputfiles'''
    try:
        float(s)
        return True
    except ValueError:
        return False


def read_field(f,fieldname,dims,intField=False):
    '''reads a single variable from an outputfile, used in various read_b2f* functions.
    f is the file object i.e f = open(fileLocation)
    fieldname is the variable name to be searched, i.e crx
    dims is the dimensions of the variable i.e [96,36]
    intField says if the variable should be returned as int or float'''

    line = f.readline().rstrip()
    # find the right field
    while fieldname not in line:
        line = f.readline().rstrip()
        if len(line) == 0:
            print('read_field: EOF reached without finding '+str(fieldname))
            print('The first variable not found is probably not in your file')
            print('Take out the search for that variable in the function doc or update your SOLPS so that the output is produced')
            return 0

    # Consistency check: number of elements specified
    # in the file should equal prod(dims)

    for i in range(len(line.split())):
        if is_number(line.split()[i]): numin = int(line.split()[i])
            
    if (numin != np.prod(dims)): 
        print(line)
        print('read_field: inconsistent number of input elements.');

    fieldVal=[]
    # collect field values
    while (len(fieldVal) != numin):
        line = f.readline().rstrip()
        for i in range(len(line.split())):
            if (intField): fieldVal.append(int(line.split()[i]))
            else:
                fieldVal.append(float(line.split()[i]))
    fieldVal=np.array(fieldVal)

    if (np.size(dims) > 1): fieldVal = fieldVal.reshape(dims,order='F').copy()
    return fieldVal


def read_b2fstate(b2fstateLoc):
    '''reads b2fstate and returns a class of the data
    b2fstateLoc is the file path to b2fgmtry i.e "/path/to/b2fstate"'''
    fieldname = 'nx,ny'
    fid = open(b2fstateLoc)
    line = fid.readline().rstrip()#read the first line
    version = line[7:17]
    print('read_b2fstate -- file version '+version)#check the version
    dim = read_field(fid,'nx,ny,ns',3,True)#check the grid size
    nx  = dim[0]
    ny  = dim[1]
    ns  = dim[2]
    fluxdim  = [nx+2,ny+2,2]
    fluxdimp = [nx+2,ny+2]
    fluxdims = [nx+2,ny+2,2,ns]
    if (version >= '03.001.000'):
        fluxdim  = [nx+2,ny+2,2,2]
        fluxdimp = fluxdim
        fluxdims = [nx+2,ny+2,2,2,ns]
    #initialize class that will hold all the data of b2fstate
    #note that it is read in order that it appears in the file, adding a variable won't work unless it is put in order that it appears in the file
    class stateResults:
        def __init__(self):
            # Read charges etc.

            self.zamin = read_field(fid,'zamin',ns)
            self.zamax = read_field(fid,'zamax',ns)
            self.zn    = read_field(fid,'zn',ns)
            self.am    = read_field(fid,'am',ns)

            # Read state variables

            self.na     = read_field(fid,'na',[nx+2,ny+2,ns])
            self.ne     = read_field(fid,'ne',[nx+2,ny+2])
            self.ua     = read_field(fid,'ua',[nx+2,ny+2,ns])
            self.uadia  = read_field(fid,'uadia',[nx+2,ny+2,2,ns])
            self.te     = read_field(fid,'te',[nx+2,ny+2])
            self.ti     = read_field(fid,'ti',[nx+2,ny+2])
            self.po     = read_field(fid,'po',[nx+2,ny+2])

            # Read fluxes

            self.fna    = read_field(fid,'fna',fluxdims)
            self.fhe    = read_field(fid,'fhe',fluxdim)
            self.fhi    = read_field(fid,'fhi',fluxdim)
            self.fch    = read_field(fid,'fch',fluxdim)
            self.fch_32 = read_field(fid,'fch_32',fluxdim)
            self.fch_52 = read_field(fid,'fch_52',fluxdim)
            self.kinrgy = read_field(fid,'kinrgy',[nx+2,ny+2,ns])
            self.time   = read_field(fid,'time',1)
            self.fch_p  = read_field(fid,'fch_p',fluxdimp)
            
    state = stateResults()#instantiate class
    # Close file
    fid.close()
    print('done reading state file')
    return state, dim



def read_b2fplasmf(fileName,nx,ny,ns):
    '''Read formatted b2fplasmf file created by B2.5 (b2run b2uf)
    returns class of SOME of the data (add what you want if it's not here)
    fileName is "/path/to/b2fplasmf" 
    nx is the poloidal cell number, including guard cells
    ny is the radial cell number, including guard cells
    ns is the number of species'''
    if not (os.path.isfile(fileName)): 
        print("b2fplasmf: Cannot find the filename")
        return 0
    fid = open(fileName)
    if (fid == -1): print("read_b2fplasmf: can't open file")


    # Get version of the b2fstate file

    line    = fid.readline().rstrip()
    version = line[7:17]

    print('read_b2fplasmf -- file version '+version)
    # Expected array sizes, gmtry
    qcdim = [nx+2,ny+2]
    if version >= '03.001.000': qcdim  = [nx+2,ny+2,2]


    # Expected array sizes, state
    fluxdim  = [nx+2,ny+2,2]
    fluxdims = [nx+2,ny+2,2,ns]
# Read basic data, there is more in b2fplasmf, I might grab it if I find out I need it 
# New variables also get added, this might be out of date, check the file if a var is missing
    class plasmfResults:
        def __init__(self):
            # Read gmtry variables

            self.crx  = read_field(fid,'crx' ,[nx+2,ny+2,4])
            self.cry  = read_field(fid,'cry' ,[nx+2,ny+2,4])

            # Read state variables

            self.fch    = read_field(fid,'fch'   ,fluxdim)
            self.fch0   = read_field(fid,'fch0'  ,fluxdim)
            self.fchp   = read_field(fid,'fchp'  ,fluxdim)
            self.fhe    = read_field(fid,'fhe'   ,fluxdim)
            self.fhe0   = read_field(fid,'fhe0'  ,fluxdim)
            self.fhep   = read_field(fid,'fhep'  ,fluxdim)
            self.fhet   = read_field(fid,'fhet'  ,fluxdim)
            self.fhi    = read_field(fid,'fhi'   ,fluxdim)
            self.fhi0   = read_field(fid,'fhi0'  ,fluxdim)
            self.fhip   = read_field(fid,'fhip'  ,fluxdim)
            self.fhit   = read_field(fid,'fhit'  ,fluxdim)
            self.fna    = read_field(fid,'fna'   ,fluxdims)
            self.fna0   = read_field(fid,'fna0'  ,fluxdims)
            self.fnap   = read_field(fid,'fnap'  ,fluxdims)
            self.fne    = read_field(fid,'fne'   ,fluxdim)
            self.fni    = read_field(fid,'fni'   ,fluxdim)
            self.na     = read_field(fid,'na'    ,[nx+2,ny+2,ns])
            self.na0    = read_field(fid,'na0'   ,[nx+2,ny+2,ns])
            self.nap    = read_field(fid,'nap'   ,[nx+2,ny+2,ns])
            self.ne     = read_field(fid,'ne'    ,[nx+2,ny+2])
            self.ne0    = read_field(fid,'ne0'   ,[nx+2,ny+2])
            self.ne2    = read_field(fid,'ne2'   ,[nx+2,ny+2])
            self.nep    = read_field(fid,'nep'   ,[nx+2,ny+2])
            self.ni     = read_field(fid,'ni'    ,[nx+2,ny+2,2])
            self.ni0    = read_field(fid,'ni0'   ,[nx+2,ny+2,2])
            self.pb     = read_field(fid,'pb'    ,[nx+2,ny+2])
            self.po     = read_field(fid,'po'    ,[nx+2,ny+2])
            self.po0    = read_field(fid,'po0'   ,[nx+2,ny+2])
            self.pop    = read_field(fid,'pop'   ,[nx+2,ny+2])
            self.te     = read_field(fid,'te'    ,[nx+2,ny+2])
            self.te0    = read_field(fid,'te0'   ,[nx+2,ny+2])
            self.tep    = read_field(fid,'tep'   ,[nx+2,ny+2])
            self.ti     = read_field(fid,'ti'    ,[nx+2,ny+2])
            self.ti0    = read_field(fid,'ti0'   ,[nx+2,ny+2])
            self.tip    = read_field(fid,'tip'   ,[nx+2,ny+2])
            self.ua     = read_field(fid,'ua'    ,[nx+2,ny+2,ns])
            self.ua0    = read_field(fid,'ua0'   ,[nx+2,ny+2,ns])
            self.uap    = read_field(fid,'uap'   ,[nx+2,ny+2,ns])
            self.uadia  = read_field(fid,'uadia' ,[nx+2,ny+2,2,ns])
            self.fchdia = read_field(fid,'fchdia',fluxdim)
            self.fmo    = read_field(fid,'fmo'   ,fluxdims)
            self.fna_32 = read_field(fid,'fna_32',fluxdims)
            self.fna_52 = read_field(fid,'fna_52',fluxdims)
            self.fni_32 = read_field(fid,'fni_32',fluxdim)
            self.fni_52 = read_field(fid,'fni_52',fluxdim)
            self.fne_32 = read_field(fid,'fne_32',fluxdim)
            self.fne_52 = read_field(fid,'fne_52',fluxdim)
            self.wadia  = read_field(fid,'wadia' ,[nx+2,ny+2,2,ns])
            self.vaecrb = read_field(fid,'vaecrb',[nx+2,ny+2,2,ns])
            self.facdrift     = read_field(fid,'facdrift'    ,[nx+2,ny+2])
            self.fac_ExB      = read_field(fid,'fac_ExB'     ,[nx+2,ny+2])
            self.fchvispar    = read_field(fid,'fchvispar'   ,fluxdim)
            self.fchvisper    = read_field(fid,'fchvisper'   ,fluxdim)
            self.fchin        = read_field(fid,'fchin'       ,fluxdim)
            self.fna_nodrift  = read_field(fid,'fna_nodrift' ,fluxdims)
            self.fac_vis      = read_field(fid,'fac_vis'     ,[nx+2,ny+2])
            self.fna_mdf      = read_field(fid,'fna_mdf'     ,fluxdims)
            self.fhe_mdf      = read_field(fid,'fhe_mdf'     ,fluxdim)
            self.fhi_mdf      = read_field(fid,'fhi_mdf'     ,fluxdim)
            self.fnaPSch      = read_field(fid,'fnaPSch'     ,fluxdims)
            self.fhePSch      = read_field(fid,'fhePSch'     ,fluxdim)
            self.fhiPSch      = read_field(fid,'fhiPSch'     ,fluxdim)
            self.fna_fcor     = read_field(fid,'fna_fcor'    ,fluxdims)
            self.fna_he       = read_field(fid,'fna_he'      ,fluxdims)
            self.fchvisq      = read_field(fid,'fchvisq'     ,fluxdim)
            self.fchinert     = read_field(fid,'fchinert'    ,fluxdim)
            if version>'3.000.006' or (version[0]=='0' and version[1:]>'3.000.006'):
              # these are only present in 3.000.007 and beyond
              self.fht          = read_field(fid,'fht'         ,fluxdim)
              self.fhj          = read_field(fid,'fhj'         ,fluxdim)
              self.fhm          = read_field(fid,'fhm'         ,fluxdims)
              self.fhp          = read_field(fid,'fhp'         ,fluxdims)
            self.resco        = read_field(fid,'resco'       ,[nx+2,ny+2,ns])
            self.reshe        = read_field(fid,'reshe'       ,[nx+2,ny+2])
            self.reshi        = read_field(fid,'reshi'       ,[nx+2,ny+2])
            self.resmo        = read_field(fid,'resmo'       ,[nx+2,ny+2,ns])
            self.resmt        = read_field(fid,'resmt'       ,[nx+2,ny+2])
            self.respo        = read_field(fid,'respo'       ,[nx+2,ny+2])
            self.sch          = read_field(fid,'sch'         ,[nx+2,ny+2,4])
            self.she          = read_field(fid,'she'         ,[nx+2,ny+2,4])
            self.shi          = read_field(fid,'shi'         ,[nx+2,ny+2,4])
            self.smo          = read_field(fid,'smo'         ,[nx+2,ny+2,4,ns])
            self.smq          = read_field(fid,'smq'         ,[nx+2,ny+2,4,ns])
            self.sna          = read_field(fid,'sna'         ,[nx+2,ny+2,2,ns])
            self.sne          = read_field(fid,'sne'         ,[nx+2,ny+2,2])
            self.rsana        = read_field(fid,'rsana'       ,[nx+2,ny+2,ns])
            self.rsahi        = read_field(fid,'rsahi'       ,[nx+2,ny+2,ns])
            self.rsamo        = read_field(fid,'rsamo'       ,[nx+2,ny+2,ns])
            self.rrana        = read_field(fid,'rrana'       ,[nx+2,ny+2,ns])
            self.rrahi        = read_field(fid,'rrahi'       ,[nx+2,ny+2,ns])
            self.rramo        = read_field(fid,'rramo'       ,[nx+2,ny+2,ns])
            self.rqahe        = read_field(fid,'rqahe'       ,[nx+2,ny+2,ns])
            self.rqrad        = read_field(fid,'rqrad'       ,[nx+2,ny+2,ns])
            self.rqbrm        = read_field(fid,'rqbrm'       ,[nx+2,ny+2,ns])
            self.rcxna        = read_field(fid,'rcxna'       ,[nx+2,ny+2,ns])
            self.rcxhi        = read_field(fid,'rcxhi'       ,[nx+2,ny+2,ns])
            self.rcxmo        = read_field(fid,'rcxmo'       ,[nx+2,ny+2,ns])
            self.b2stbr_sna   = read_field(fid,'b2stbr_sna'  ,[nx+2,ny+2,ns])
            self.b2stbr_smo   = read_field(fid,'b2stbr_smo'  ,[nx+2,ny+2,ns])
            self.b2stbr_she   = read_field(fid,'b2stbr_she'  ,[nx+2,ny+2])
            self.b2stbr_shi   = read_field(fid,'b2stbr_shi'  ,[nx+2,ny+2])
            self.b2stbr_sch   = read_field(fid,'b2stbr_sch'  ,[nx+2,ny+2])
            self.b2stbr_sne   = read_field(fid,'b2stbr_sne'  ,[nx+2,ny+2])
            self.b2stbc_sna   = read_field(fid,'b2stbc_sna'  ,[nx+2,ny+2,ns])
            self.b2stbc_smo   = read_field(fid,'b2stbc_smo'  ,[nx+2,ny+2,ns])
            self.b2stbc_she   = read_field(fid,'b2stbc_she'  ,[nx+2,ny+2])
            self.b2stbc_shi   = read_field(fid,'b2stbc_shi'  ,[nx+2,ny+2])
            self.b2stbc_sch   = read_field(fid,'b2stbc_sch'  ,[nx+2,ny+2])
            self.b2stbc_sne   = read_field(fid,'b2stbc_sne'  ,[nx+2,ny+2])
            self.b2stbm_sna   = read_field(fid,'b2stbm_sna'  ,[nx+2,ny+2,ns])
            self.b2stbm_smo   = read_field(fid,'b2stbm_smo'  ,[nx+2,ny+2,ns])
            self.b2stbm_she   = read_field(fid,'b2stbm_she'  ,[nx+2,ny+2])
            self.b2stbm_shi   = read_field(fid,'b2stbm_shi'  ,[nx+2,ny+2])
            self.b2stbm_sch   = read_field(fid,'b2stbm_sch'  ,[nx+2,ny+2])
            self.b2stbm_sne   = read_field(fid,'b2stbm_sne'  ,[nx+2,ny+2])
            self.b2sihs_divue = read_field(fid,'b2sihs_divue',[nx+2,ny+2])
            self.b2sihs_divua = read_field(fid,'b2sihs_divua',[nx+2,ny+2])
            self.b2sihs_exbe  = read_field(fid,'b2sihs_exbe' ,[nx+2,ny+2])
            self.b2sihs_exba  = read_field(fid,'b2sihs_exba' ,[nx+2,ny+2])
            self.b2sihs_visa  = read_field(fid,'b2sihs_visa' ,[nx+2,ny+2])
            self.b2sihs_joule = read_field(fid,'b2sihs_joule',[nx+2,ny+2])
            self.b2sihs_fraa  = read_field(fid,'b2sihs_fraa' ,[nx+2,ny+2])
            self.b2sihs_str   = read_field(fid,'b2sihs_str' ,[nx+2,ny+2])
            self.b2npmo_smaf  = read_field(fid,'b2npmo_smaf' ,[nx+2,ny+2,4,ns])
            self.b2npmo_smag  = read_field(fid,'b2npmo_smag' ,[nx+2,ny+2,4,ns])
            self.b2npmo_smav  = read_field(fid,'b2npmo_smav' ,[nx+2,ny+2,4,ns])
            self.smpr         = read_field(fid,'smpr'        ,[nx+2,ny+2,ns])
            self.smpt         = read_field(fid,'smpt'        ,[nx+2,ny+2,ns])
            self.smfr         = read_field(fid,'smfr'        ,[nx+2,ny+2,ns])
            self.smcf         = read_field(fid,'smcf'        ,[nx+2,ny+2,ns])
            self.ext_sna      = read_field(fid,'ext_sna'     ,[nx+2,ny+2,ns])
            self.ext_smo      = read_field(fid,'ext_smo'     ,[nx+2,ny+2,ns])
            self.ext_she      = read_field(fid,'ext_she'     ,[nx+2,ny+2])
            self.ext_shi      = read_field(fid,'ext_shi'     ,[nx+2,ny+2])
            self.ext_sch      = read_field(fid,'ext_sch'     ,[nx+2,ny+2])
            self.ext_sne      = read_field(fid,'ext_sne'     ,[nx+2,ny+2])
            self.calf         = read_field(fid,'calf'        ,fluxdim)
            self.cdna         = read_field(fid,'cdna'        ,fluxdims)
            self.cdpa         = read_field(fid,'cdpa'        ,fluxdims)
            self.ceqp         = read_field(fid,'ceqp'        ,[nx+2,ny+2])
            self.chce         = read_field(fid,'chce'        ,fluxdim)
            self.chci         = read_field(fid,'chci'        ,fluxdim)
            self.chve         = read_field(fid,'chve'        ,fluxdim)
            self.chvemx       = read_field(fid,'chvemx'      ,[nx+2,ny+2])
            self.chvi         = read_field(fid,'chvi'        ,fluxdim)
            self.chvimx       = read_field(fid,'chvimx'      ,[nx+2,ny+2])
            self.csig         = read_field(fid,'csig'        ,fluxdim)
            self.cvla         = read_field(fid,'cvla'        ,fluxdims)
            self.cvsa         = read_field(fid,'cvsa'        ,fluxdims)
            self.cthe         = read_field(fid,'cthe'        ,[nx+2,ny+2,ns])
            self.cthi         = read_field(fid,'cthi'        ,[nx+2,ny+2,ns])
            self.csigin       = read_field(fid,'csigin'      ,[fluxdims[0],fluxdims[1],fluxdims[2],fluxdims[3],ns])
            self.cvsa_cl      = read_field(fid,'cvsa_cl'     ,fluxdims)
            self.fllime       = read_field(fid,'fllime'      ,[nx+2,ny+2])
            self.fllimi       = read_field(fid,'fllimi'      ,[nx+2,ny+2])
            self.fllim0fna    = read_field(fid,'fllim0fna'   ,fluxdims)
            self.fllim0fhi    = read_field(fid,'fllim0fhi'   ,fluxdims)
            self.fllimvisc    = read_field(fid,'fllimvisc'   ,[nx+2,ny+2,ns])
            self.sig0         = read_field(fid,'sig0'        ,[nx+2,ny+2])
            self.hce0         = read_field(fid,'hce0'        ,[nx+2,ny+2])
            self.alf0         = read_field(fid,'alf0'        ,[nx+2,ny+2])
            self.hci0         = read_field(fid,'hci0'        ,[nx+2,ny+2])
            self.hcib         = read_field(fid,'hcib'        ,[nx+2,ny+2,ns])
            self.dpa0         = read_field(fid,'dpa0'        ,[nx+2,ny+2,ns])
            self.dna0         = read_field(fid,'dna0'        ,[nx+2,ny+2,ns])
            self.vsa0         = read_field(fid,'vsa0'        ,[nx+2,ny+2,ns])
            self.vla0         = read_field(fid,'vla0'        ,[nx+2,ny+2,2,ns])
            self.csig_an      = read_field(fid,'csig_an'     ,fluxdim)
            self.calf_an      = read_field(fid,'calf_an'     ,fluxdim)
            nstra              = read_field(fid,'nstra'       ,[1],True)
            if nstra!=0:
              self.sclstra      = read_field(fid,'sclstra'     ,[ns+1,nstra[0]])
              self.sclrtio      = read_field(fid,'sclrtio'     ,[ns+1,nstra[0]])
    plasmf = plasmfResults()
    fid.close()
    print('done reading b2fplasmf')
    return plasmf



"""

The following functions (read_transport_files)  comes from SOLPSutils code written by 
Robert Wilcox and Jeremy Lore

"""


def read_b2fstate_Bob(fname):
    if not os.path.exists(fname):
        print('ERROR: b2fstate file not found: ',fname)
        return None

    DEBUG = False

    data = []
    with open(fname, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):

        # Special handling for first few lines
        if i == 0:
            version = line.split()[0][7:-1]
            state = {'version':version}
            continue
        elif i == 1:
            continue
        elif i == 2:
            # Assume starts with nx,ny,ns after version
            state['nx'] = int(line.split()[0])
            state['ny'] = int(line.split()[1])
            state['ns'] = int(line.split()[2])
            numcells = (state['nx']+2)*(state['ny']+2)
            continue

        if line.split()[0] == '*cf:':
            vartype = line.split()[1]
            varsize = int(line.split()[2])
            varname = line.split()[3]
            if DEBUG:
                print(varname,vartype,varsize,state['nx'],state['ny'],state['ns'],numcells)
            # Some variables have no entries depending on config
            if varsize == 0:
                state[varname] = None
            data = []
        else:
            # Parse by type
            if vartype == "char":
                state[varname] = line.strip()
            else:
                splitline = line.split()
                for value in splitline:
                    if vartype == "int":
                        data.append(int(value))
                    else:
                        data.append(float(value))

                if len(data) == varsize:
                    if varsize == numcells:
                        # This is a scalar quantity
                        state[varname] = np.array(data).reshape([state['nx']+2,state['ny']+2], order = 'F')
                    elif varsize == 2*numcells:
                        # This is a flux quantity
                        state[varname] = np.array(data).reshape([state['nx']+2,state['ny']+2,2], order = 'F')
                    elif varsize == numcells*state['ns']:
                        # This is a scalar quantity by species
                        state[varname] = np.array(data).reshape([state['nx']+2,state['ny']+2,state['ns']], order = 'F')
                    elif varsize == 2*numcells*state['ns']:
                        # This is a flux quantity by species
                        state[varname] = np.array(data).reshape([state['nx']+2,state['ny']+2,2,state['ns']], order = 'F')
                    elif varsize == 4*numcells:
                        # This is a flux quantity in 3.1 format
                        state[varname] = np.array(data).reshape([state['nx']+2,state['ny']+2,2,2], order = 'F')
                    elif varsize == 4*numcells*state['ns']:
                        # This is a flux quantity by species in 3.1 format
                        state[varname] = np.array(data).reshape([state['nx']+2,state['ny']+2,2,2,state['ns']], order = 'F')
                    elif varsize%numcells == 0:
                            print("Warning, must have missed some dimension checks for variable:",varname)
                    else:
                        # For other dimensions assign as is (e.g., zamin)
                        state[varname] = np.array(data)
    return state



def read_iout_method(fdir, fname, nx, ny):
    
    if not os.path.exists(fdir):
        print('ERROR: {} file not found: '.format(fname), fdir)
        return None

    DEBUG = False

    data = np.zeros([ny+2, nx+2])
    with open(fdir, 'r') as f:
        lines = f.readlines()
    
    
    for i, line in enumerate(reversed(lines)):
        splitline = line.split()
        value_list = []
        
        for k, value in enumerate(splitline):
            if k == 0:
                pass
            else:
                value_list.append(float(value))
                
            # if value.type() is int:
            #     data.append(int(value))
            # else:
                
        # print(len(value_list))
        if i < ny+2:
            data[i, :] = value_list
        
        iout = data[1:37, 1:97]
            
        
        # iout = np.array(data).reshape([nx+2, ny+2], order = 'F')
    
    return iout
        

        
def read_b2wdat_field(filename):
# reads a .out file produced by setting b2wdat_iout='4' in b2mn.dat
    f = open(filename)
    line = f.readline().rstrip().split()
    fieldVal = []
    while (line!=[]):
        line = f.readline().rstrip().split()
        if line==[]: break
        fieldVal.append([float(i) for i in line][1:])
    return np.array(fieldVal[::-1]).T


def is_neutral(a):
# checks if the species is a neutral 
# DOESNT WORK WITH MORE THAN TWO TYPES OF IONS
    if a>=6: print('WARNING: bigger species index than is_neutral was made for, proceed with caution')
    if a==0 or a==2:
        return True
    else:
        return False




def read_b2wdat(b2wdatLoc,nSpec):
# reads .out files produced by setting b2wdat_iout='4' in b2mn.dat and returns a class with the data
# currently only grabs what I have needed there are literally hundreds more
# this isn't a very robust function, might not work if it isn't a D-only or D+Li case. Be careful
# adjusting is_neutral to be more robust might be all it needs but not sure
    nas = []
    uas = []
    ues = []
    b2srdt_smodts = []
    b2npmo_fmoxs = []
    b2npmo_fmoys = []
    b2sigp_smogpis = []
    b2sigp_smogpos = []
    b2npmo_smbs = []
    b2stcx_smqs = []
    b2npmo_smocfs = []
    b2stel_smq_ions = []
    b2stel_smq_recs = []
    b2npmo_smotfias = []
    b2npmo_smotfeas = []
    b2npmo_smofreas = []
    b2npmo_smofrias = []
    b2npmo_smoans = []
    b2stbc_phys_smos = []
    b2npmo_smovvs = []
    b2npmo_smovhs = []
    b2stbr_smos = []
    b2trcl_luciani_fllim_cvsahzxs = []
    b2trcl_luciani_cvsahzxs = []
    b2npmo_resmos = []
    b2stbr_sna_eirs = []
    b2stel_sna_ions = []
    b2stel_sna_recs = []
    b2npc_snas = []
    b2npc_fnaxs = []
    b2npc_fnays = []
    b2tfnb_dPat_mdf_gradnax = []
    b2trno_cdnax = []
    crxs = []
    crys = []
    for i in range(4): crxs.append(read_b2wdat_field(b2wdatLoc+'output/'+'crx'+str(i)+'.dat'))
    for i in range(4): crys.append(read_b2wdat_field(b2wdatLoc+'output/'+'cry'+str(i)+'.dat'))
    for spIdx in range(nSpec):
        if is_neutral(spIdx): continue
        if spIdx==1:nas.append(read_b2wdat_field(b2wdatLoc+'output/'+'b2npc11_na00'+str(spIdx)+'.dat'))
        if spIdx>1: nas.append(read_b2wdat_field(b2wdatLoc+'output/'+'b2npc11_na00'+str(spIdx)+'.dat'))
        uas.append(np.array(read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_ua00'+str(spIdx)+'.dat')))
        b2srdt_smodts.append( -read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_madnavadt00'+str(spIdx)+'.dat'))
        b2npmo_fmoxs.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_fmox00'+str(spIdx)+'.dat'))
        b2npmo_fmoys.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_fmoy00'+str(spIdx)+'.dat'))
        b2sigp_smogpis.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2sigp_smogpi00'+str(spIdx)+'.dat'))
        b2sigp_smogpos.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2sigp_smogpo00'+str(spIdx)+'.dat'))
        b2npmo_smbs.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_smb00'+str(spIdx)+'.dat'))
        b2stcx_smqs.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2stcx_smq00'+str(spIdx)+'.dat'))
        b2npmo_smocfs.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_smocf00'+str(spIdx)+'.dat'))
        b2stel_smq_ions.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2stel_smq_ion00'+str(spIdx)+'.dat'))
        b2stel_smq_recs.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2stel_smq_rec00'+str(spIdx)+'.dat'))
        b2npmo_smotfias.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_smotfia00'+str(spIdx)+'.dat'))
        b2npmo_smotfeas.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_smotfea00'+str(spIdx)+'.dat'))
        b2npmo_smofreas.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_smofrea00'+str(spIdx)+'.dat'))
        b2npmo_smofrias.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_smofria00'+str(spIdx)+'.dat'))
        b2npmo_smoans.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_smoan00'+str(spIdx)+'.dat'))
        b2stbc_phys_smos.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2stbc_phys_smo00'+str(spIdx)+'.dat'))
        b2npmo_smovvs.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_smovv00'+str(spIdx)+'.dat'))
        b2npmo_smovhs.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_smovh00'+str(spIdx)+'.dat'))
        b2stbr_smos.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2stbr_smo_eir00'+str(spIdx)+'.dat'))
        b2trcl_luciani_fllim_cvsahzxs.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2trcl_luciani_fllim_cvsahzx00'+str(spIdx)+'.dat'))
        b2trcl_luciani_cvsahzxs.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2trcl_luciani_cvsahzx00'+str(spIdx)+'.dat'))
        b2npmo_resmos.append( read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_resmo00'+str(spIdx)+'.dat'))
        b2stbr_sna_eirs.append(read_b2wdat_field(b2wdatLoc+'output/'+'b2stbr_sna_eir00'+str(spIdx)+'.dat'))
        b2stel_sna_ions.append(read_b2wdat_field(b2wdatLoc+'output/'+'b2stel_sna_ion00'+str(spIdx)+'.dat'))
        b2stel_sna_recs.append(read_b2wdat_field(b2wdatLoc+'output/'+'b2stel_sna_rec00'+str(spIdx)+'.dat'))
        b2tfnb_dPat_mdf_gradnax.append(read_b2wdat_field(b2wdatLoc+'output/'+'b2tfnb_dPat_mdf_gradnax00'+str(spIdx)+'.dat'))
        b2trno_cdnax.append(read_b2wdat_field(b2wdatLoc+'output/'+'b2trno_cdnax00'+str(spIdx)+'.dat'))
        if spIdx==0 or spIdx==2: b2npc_fnaxs.append(read_b2wdat_field(b2wdatLoc+'output/'+'b2npco_fnax00'+str(spIdx)+'.dat'))
        else: b2npc_fnaxs.append(read_b2wdat_field(b2wdatLoc+'output/'+'b2npc11_fnax00'+str(spIdx)+'.dat'))

        if spIdx==0 or spIdx==2: b2npc_fnays.append(read_b2wdat_field(b2wdatLoc+'output/'+'b2npco_fnaxy00'+str(spIdx)+'.dat'))
        else: b2npc_fnays.append(read_b2wdat_field(b2wdatLoc+'output/'+'b2npc11_fnay00'+str(spIdx)+'.dat'))

        if spIdx==0 or spIdx==2: b2npc_snas.append(read_b2wdat_field(b2wdatLoc+'output/'+'b2npc_sna00'+str(spIdx)+'.dat'))
        else: b2npc_snas.append(read_b2wdat_field(b2wdatLoc+'output/'+'b2npc11_sna00'+str(spIdx)+'.dat'))
    #initialize class that will hold all the data of bwdat output
    class b2wdatResults:
        def __init__(self):
            #LHS of the momentum eqn
            self.b2srdt_smodt = b2srdt_smodts
            self.b2npmo_fmox = b2npmo_fmoxs
            self.b2npmo_fmoy = b2npmo_fmoys
            self.b2sigp_smogpi = b2sigp_smogpis
            self.b2sigp_smogpo = b2sigp_smogpos

            #RHS of the momentum equation
            self.b2npmo_smb = b2npmo_smbs
            self.b2stcx_smq = b2stcx_smqs
            self.b2npmo_smocf = b2npmo_smocfs
            self.b2stel_smq_ion = b2stel_smq_ions
            self.b2stel_smq_rec = b2stel_smq_recs
            self.b2npmo_smotfia = b2npmo_smotfias
            self.b2npmo_smotfea = b2npmo_smotfeas
            self.b2npmo_smofrea = b2npmo_smofreas
            self.b2npmo_smofria = b2npmo_smofrias
            self.b2npmo_smoan = b2npmo_smoans
            self.b2stbc_phys_smo = b2stbc_phys_smos
            self.b2npmo_smovv = b2npmo_smovvs
            self.b2npmo_smovh = b2npmo_smovhs
            self.b2stbr_smo = b2stbr_smos
            self.b2trcl_luciani_fllim_cvsahzx = b2trcl_luciani_fllim_cvsahzxs
            self.b2trcl_luciani_cvsahzx = b2trcl_luciani_cvsahzxs
            self.b2npmo_resmo = b2npmo_resmos

            #energy balance terms
            self.b2stbr_shi_eir = read_b2wdat_field(b2wdatLoc+'output/'+'b2stbr_shi_eir.dat')
            self.b2stbr_she_eir = read_b2wdat_field(b2wdatLoc+'output/'+'b2stbr_she_eir.dat')
            self.b2stel_she_rad = read_b2wdat_field(b2wdatLoc+'output/'+'b2stel_she_rad.dat')

            #particle sources
            self.b2stbr_sna_eir = b2stbr_sna_eirs
            self.b2stel_sna_ion = b2stel_sna_ions
            self.b2stel_sna_rec = b2stel_sna_recs
            self.b2npc_sna = b2npc_snas
            self.b2npc_fnaxs = b2npc_fnaxs
            self.b2npc_fnays = b2npc_fnays

            #geo info
            self.hx  = read_b2wdat_field(b2wdatLoc+'output/'+'hx.dat')
            self.hy  = read_b2wdat_field(b2wdatLoc+'output/'+'hy.dat')
            self.hz  = read_b2wdat_field(b2wdatLoc+'output/'+'hz.dat')
            self.vol = read_b2wdat_field(b2wdatLoc+'output/'+'vol.dat')
            self.bbx = read_b2wdat_field(b2wdatLoc+'output/'+'bbx.dat')
            self.bb  = read_b2wdat_field(b2wdatLoc+'output/'+'bb.dat')
            self.bx  = read_b2wdat_field(b2wdatLoc+'output/'+'bbx.dat')/read_b2wdat_field(b2wdatLoc+'output/'+'bb.dat')
            self.crx = crxs
            self.cry = crys

            #plasma parameters
            self.na = nas
            self.ua = uas
            self.ue = read_b2wdat_field(b2wdatLoc+'output/'+'b2npmo_ue.dat')
            self.phi = read_b2wdat_field(b2wdatLoc+'output/'+'b2nppo_po.dat')
            self.ti = read_b2wdat_field(b2wdatLoc+'output/'+'ti_eV.dat')
            #misc
            self.fchanmly = read_b2wdat_field(b2wdatLoc+'output/'+'b2tfch__fchanmly.dat')
            self.fchanmlx = read_b2wdat_field(b2wdatLoc+'output/'+'b2tfch__fchanmlx.dat')
            self.fchdiay  = read_b2wdat_field(b2wdatLoc+'output/'+'b2tfch__fchdiay.dat')
            self.fchdiax  = read_b2wdat_field(b2wdatLoc+'output/'+'b2tfch__fchdiax.dat')

    b2wdat = b2wdatResults()#instantiate class
    print('done reading b2wdat files')
    return b2wdat

    
    
    








