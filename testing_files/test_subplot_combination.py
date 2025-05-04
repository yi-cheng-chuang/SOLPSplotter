# -*- coding: utf-8 -*-
"""
Created on Fri May  2 13:07:22 2025

@author: ychuang
"""



import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 10, 100)
datasets = [np.sin(x), np.cos(x), np.sin(2*x)]

# Set up the figure and GridSpec
fig = plt.figure()
gs = fig.add_gridspec(3, 1)

# First two subplots with shared x-axis
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)

# Third subplot independent
ax2 = fig.add_subplot(gs[2, 0])  # No sharex

# Plotting
ax0.plot(x, datasets[0])
ax1.plot(x, datasets[1])
ax2.plot(x, datasets[2])

ax0.set_title("Plot 0 (shared x)")
ax1.set_title("Plot 1 (shared x with Plot 0)")
ax2.set_title("Plot 2 (independent x)")

# Hide x tick labels for the first two (optional, for cleaner look)
plt.setp(ax0.get_xticklabels(), visible=False)
plt.setp(ax1.get_xticklabels(), visible=False)

plt.tight_layout()
plt.show()







x = np.linspace(0, 10, 100)
y = np.sin(x)

plt.plot(x, y, color='blue', alpha=0.3, linewidth=2)  # transparent line
plt.title("Transparent Line Example")
plt.show()
