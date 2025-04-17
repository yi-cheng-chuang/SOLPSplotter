# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 04:53:01 2025

@author: ychuang
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Data
A = [[1, 4], [7, 10]]
B = [[2, 6], [8, 9]]
C = [[11, 13], [16, 17]]
D = [[1, 2], [3, 4]]

# Plot settings
datasets = {'A': (A, 'red'), 'B': (B, 'blue'), 'C': (C, 'green')}
markers = ['^', 'o']  # triangle, circle

# Plotting
fig, ax = plt.subplots()

for label, (data, color) in datasets.items():
    for i in range(2):  # Two sublists
        ax.plot(D[i], data[i], marker=markers[i], color=color, linestyle='', 
                label=f'{label} - {"Group 1" if i == 0 else "Group 2"}')

# Remove duplicate legends
handles, labels = ax.get_legend_handles_labels()
unique = dict(zip(labels, handles))
ax.legend(unique.values(), unique.keys())

# Labels and title
ax.set_xlabel('D values')
ax.set_ylabel('Data values')
ax.set_title('Plot of A, B, C vs D')

plt.grid(True)
plt.show()


# Data
A = [[1, 4], [7, 10]]
B = [[2, 6], [8, 9]]
C = [[11, 13], [16, 17]]
D = [[1, 2], [3, 4]]

# Datasets and styles
datasets = {'A': (A, 'red'), 'B': (B, 'blue'), 'C': (C, 'green')}
markers = ['^', 'o']  # triangle for Group 1, circle for Group 2

fig, ax = plt.subplots()

# Plot the actual data
for label, (data, color) in datasets.items():
    for i in range(2):  # two groups
        ax.plot(D[i], data[i], marker=markers[i], color=color, linestyle='')

# Create dummy handles for color legend (dataset)
color_handles = [
    Line2D([0], [0], color='red', marker='o', linestyle='', label='A'),
    Line2D([0], [0], color='blue', marker='o', linestyle='', label='B'),
    Line2D([0], [0], color='green', marker='o', linestyle='', label='C'),
]

# Create dummy handles for marker legend (group)
marker_handles = [
    Line2D([0], [0], color='black', marker='^', linestyle='', label='Group 1'),
    Line2D([0], [0], color='black', marker='o', linestyle='', label='Group 2'),
]

# Add the two legends separately
legend1 = ax.legend(handles=color_handles, title='Dataset (Color)', loc='upper left')
legend2 = ax.legend(handles=marker_handles, title='Group (Marker)', loc='upper right')

# Add both legends to the plot
ax.add_artist(legend1)

# Labeling
ax.set_xlabel('D values')
ax.set_ylabel('Data values')
ax.set_title('Split Legend: Color (Dataset) & Marker (Group)')

plt.grid(True)
plt.show()

