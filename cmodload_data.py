# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 18:46:59 2023

@author: user
"""

elif DEV=='cmod':
    
    if ROOTSHOT == '':
        
        BASEDRT = '{}cmod/{}home'.format(BASEDRT,Shot)
        GFILE = glob.glob('{}gfileProcessing/cmod_files/g{}*'.format(TOPDRT,Shot))
        print(GFILE)
        GF = eq.equilibrium(gfile=GFILE[-1])
        if EXP:
            ExpFile = Shot
            ExpData = loadmat('{}gfileProcessing/cmod_files/{}.mat'.format(TOPDRT, ExpFile))
    
    else: 
        
        BASEDRT = '{}cmod/0{}home'.format(BASEDRT, Shot[-2:])
        
        #GFILE = '{}gfileProcessing/cmod_files/g11607180{}.01209_974'.format(TOPDRT, Shot[-2:])
        GFILE = glob.glob('{}gfileProcessing/cmod_files/g{}*{}*'.format(TOPDRT,ROOTSHOT,Shot[-2:]))
        print(GFILE)
        GF = eq.equilibrium(gfile=GFILE[0])
        
        #ExpFile = '11607180{}'.format(Shot[-2:])
        ExpFile = glob.glob('{}gfileProcessing/cmod_files/{}*{}.mat'.format(TOPDRT,ROOTSHOT,Shot[-2:]))
        #print(ExpFile)
        ExpData = loadmat(ExpFile[0])
        
    if EXP:
        TR=[ii for ii, tt in enumerate(ExpData['time'][0,:]) if tt > TimeRange[0] and tt < TimeRange[1]]
        
        Rmid = ExpData['rmid'][:,TR]
        RmidAvg = np.mean(Rmid, axis=1)

        try:
            Psin = ExpData['psin'][:,TR]
            PsinAvg = np.mean(Psin, axis=1)
        except:
            print('Psin coordinates not found. Attempting to approximate experimental psin from gfile')
            PsinAvg = GF.psiN(RmidAvg,np.zeros(RmidAvg.shape))[0]
        
        #Robust Statistics -> Use MEDIAN, not MEAN, of TS Data -> Need to filter out zeros before taking Median!
        
        Nemid = ExpData['ne'][:,TR]
        NemidAvg=np.zeros(Nemid.shape[0])
        for ii in range(Nemid.shape[0]):
            if np.sum(Nemid[ii,:]) == 0:
                continue
            else:
                NemidAvg[ii]=np.median(Nemid[ii][np.nonzero(Nemid[ii])])
                
        ErrNe=np.median(ExpData['nerr'][:,TR], axis=1)
        '''
        Nemid[Nemid == 0] = np.nan
        NemidAvg = np.nanmean(Nemid, axis=1)
        ErrNe = np.nanmean(ExpData['nerr'][:,ti:tf], axis=1)
        '''
        NeThresh = (ErrNe*2)/NemidAvg

        for NT in range(len(NeThresh)):
            if np.abs(NeThresh[NT]) > 5.0:
                NemidAvg[NT] = np.nan
                ErrNe[NT] = np.nan
        
        Temid = ExpData['te'][:,TR]
        TemidAvg=np.zeros(Temid.shape[0])
        for ii in range(Temid.shape[0]):
            if np.sum(Temid[ii,:]) == 0:
                continue
            else:
                TemidAvg[ii]=np.median(Temid[ii][np.nonzero(Temid[ii])])
        ErrTe=np.median(ExpData['terr'][:,TR], axis=1)
        '''
        Temid[Temid == 0] = np.nan
        TemidAvg = np.nanmean(Temid, axis=1)
        ErrTe = np.nanmean(ExpData['terr'][:,ti:tf], axis=1)
        '''
        TeThresh = (ErrTe*2)/TemidAvg
        
        for TT in range(len(TeThresh)):
            if np.abs(TeThresh[TT]) > 5.0:
                TemidAvg[TT] = np.nan
                ErrTe[TT] = np.nan
        
        
        self.ExpDict['NemidAvg'] = NemidAvg
        self.ExpDict['ErrNe'] = ErrNe
        self.ExpDict['TemidAvg'] = TemidAvg
        self.ExpDict['ErrTe'] = ErrTe