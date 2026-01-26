# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 23:47:06 2026

@author: ychuang
"""

import pickle
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np


class load_Matteo_pkl:

    def __init__(self, DF, data):

        self.DF = DF
        self.data = data

    def load_ft44pkl(self):

        data_dic = {}

        for aa in self.data['dircomp']['Attempt'].keys():

            simu_dir = self.data['dirdata']['simudir'][aa]

            npz_path = "{}/ft44.npz".format(simu_dir)

            data = dict(np.load(npz_path, allow_pickle=True))

            data_dic[aa] = data

        self.data['Matteo_ft44'] = data_dic
