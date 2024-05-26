# -*- coding: utf-8 -*-
"""
Created on Sun May 26 14:30:55 2024

@author: ychuang
"""

import numpy as np
import SOLPS_set as sps
from SOLPS_load_directory import load_directory 
import load_mast_expdata_method as lmem
import load_B2_data_method as lbdm
import load_coord_method as lcm
import fitting_method as fm 
from scipy import interpolate



simu_loc = 'C:/Users/ychuang/Documents/SOLPS_data/simulation_data/mast/027205'
simu_case = 'AD3D'
file_dir = '{}/{}/baserun'.format(simu_loc, simu_case)
geo_dir = '{}/b2fgmtry'.format(file_dir)


simu_loc = 'C:/Users/ychuang/Documents/SOLPS_data/simulation_data/mast/027205'
simu_case_2 = 'org_new_series'
file_dir_2 = '{}/{}/baserun'.format(simu_loc, simu_case)
geo_dir_2 = '{}/b2fgmtry'.format(file_dir)



def load_b2fgmtry(file_name):
          
    try:
        geo = lcm.read_b2fgmtry(file_name)
        # print(type(geo))
    except:
        print('can not generate geo')
    
    return geo
    
    

geo = load_b2fgmtry(file_name= geo_dir)
geo_org = load_b2fgmtry(file_name= geo_dir_2)






