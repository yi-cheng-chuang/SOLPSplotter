# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 22:34:31 2023

@author: user
"""

import xarray as xr
import numpy as np
import load_mast_expdata_method as lmem

a = np.zeros([2, 2, 1])
b = np.ones([2, 2])
a[:, :, 0] = b

pts = np.array([[2.1, 6.2, 8.3],
                [3.3, 5.2, 7.1]])


k = lmem.mast_series_dir()



# P0=np.array([RadLoc.loc[1,PolPos,Attempt[-1]].values,VertLoc.loc[1,PolPos,Attempt[-1]].values])
# P1 = np.array([RadLoc.loc[SEP,PolPos,Attempt[-1]].values,VertLoc.loc[SEP,PolPos,Attempt[-1]].values])
# RR_SEP = np.sqrt((ChordX-P1[0])**2 + 
# (ChordY-P1[1])**2)*np.sign(np.sqrt((ChordX-P0[0])**2 + (ChordY-P0[1])**2)-R_SEP)

# RR = binned_statistic(RR_SEP, RR_SEP,statistic='mean',bins=36)[0]
# RR_SEP_avg[int(pol_loc)] = RR[~np.isnan(RR)]

# RR_SEP_avg = {}
# Psin_avg = {}
# Ne_SOLPS_med = {}
# NeuDen_SOLPS_med = {}

# efold={}
# efold_adj={}
# efold_adj_err={}
# yparam={}
# eparam={}
# e_err={}
# efold_err={}
# x0 = {}
# xi = {}
# fluxpsn = {}
# fluxpsnparam={}
# fluxpsn_err={}

# zmask=np.where(NeuDen_SOLPS_med[int(pol_loc)][RR_i] != 0)

# exfit = curve_fit(fm.expfit, RR_exp[zmask],NeuDen_SOLPS_med[int(pol_loc)][RR_i][zmask],e0)
# eparam[pol_loc] = exfit[0]
# e_err[int(pol_loc)] = np.sqrt(np.diag(exfit[1]))
# efold[pol_loc] = 1000/eparam[int(pol_loc)][1]
# efold_err[int(pol_loc)] = efold[int(pol_loc)]*(e_err[int(pol_loc)][1]/eparam[int(pol_loc)][1])




# Mask=Bounds.contains_points(np.vstack((PR,PZ)).T)
# ChordX=np.ma.array(PR,mask=~Mask).compressed()
# ChordY=np.ma.array(PZ,mask=~Mask).compressed()

# RR_SEP = np.sqrt((ChordX-P1[0])**2 + (ChordY-P1[1])**2)*np.sign(np.sqrt((ChordX-P0[0])**2 + (ChordY-P0[1])**2)-R_SEP)

# RR = binned_statistic(RR_SEP, RR_SEP,statistic='mean',bins=36)[0]
# RR_SEP_avg[int(pol_loc)] = RR[~np.isnan(RR)]

# Psin_ALL=Psin.loc[:,:,Attempt[-1]].values.flatten()
# Psin_SOLPS = np.ma.array(Psin_ALL,mask=~Mask).compressed()
# PS = binned_statistic(RR_SEP,Psin_SOLPS,statistic='mean',bins=36)[0] 
# Psin_avg[int(pol_loc)] = PS[~np.isnan(PS)]

# RR_i = np.where((RR_SEP_avg[int(pol_loc)]>(xi[int(pol_loc)])) & 
#                 (RR_SEP_avg[int(pol_loc)]<(x0[int(pol_loc)])))[0]
# RR_i = np.arange((RR_i[0]-1),(RR_i[-1]+2),1)
# while np.count_nonzero(NeuDen_SOLPS_med[int(pol_loc)][RR_i]) < 3:
#     RR_i = np.arange((RR_i[0]),(RR_i[-1]+2),1)
#     print(RR_i)
# RR_exp= RR_SEP_avg[int(pol_loc)][RR_i] 

# fluxpsnparam[int(pol_loc)] = np.polyfit(RR_exp, Psin_avg[int(pol_loc)][RR_i], 1, cov=True)
# fluxpsn[PolPos] = fluxpsnparam[PolPos][0][0]/fluxpsnparam[JXA][0][0]
# efold_adj[PolPos] = fluxpsn[PolPos]*efold[PolPos]
