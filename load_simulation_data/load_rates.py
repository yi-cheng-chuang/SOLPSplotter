# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 16:00:40 2026

@author: ychuang
"""

import pickle

pkl_path = "/mnt/data/380afcf2-26e6-4aee-b994-e1cf2f32a8bc.pkl"

with open(pkl_path, "rb") as f:
    data = pickle.load(f)


print(type(data))
