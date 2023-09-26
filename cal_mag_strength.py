# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 12:04:08 2023

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt



r_0 = 0.2
r_I = 0.85
a = 1
R_left = np.linspace(0.2, 0.8, 53)
mag_magnitude_left = np.sqrt((1/(r_0 + R_left))**2 + (1/(abs(R_left - r_I)))**2)
R_right = np.linspace(0.9, 1.5, 53)
mag_magnitude_right = np.sqrt((1/(r_0 + R_right))**2 + (1/(abs(R_right - r_I)))**2)


plt.plot(R_left, mag_magnitude_left, 'o-', color = 'g')
plt.plot(R_right, mag_magnitude_right, 'o-', color = 'g')
plt.show()




