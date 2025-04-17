# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 12:08:46 2025

@author: ychuang
"""



import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


test_list = np.linspace(250, -100, 8)
print(test_list)




# Scatter plot with triangle markers
plt.scatter([1, 2, 3], [4, 5, 6], marker='^', label='Data', s=100)

# Create a custom legend handle
triangle_marker = mlines.Line2D([], [], color='blue', marker='^',
                                 linestyle='None', markersize=10, label='Triangle Marker')

plt.legend(handles=[triangle_marker])
plt.show()
