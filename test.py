# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 22:34:31 2023

@author: user
"""

import xarray as xr
import numpy as np

a = np.zeros([2, 2, 1])
b = np.ones([2, 2])
a[:, :, 0] = b
