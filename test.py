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


# k = lmem.mast_series_dir()





c = [3, 50, 7]
c_total = [3, 4, 7, 9, 10, 19, 20, 21, 49, 50]

print(all([item in c_total for item in c]))





