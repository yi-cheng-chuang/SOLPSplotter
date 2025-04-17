# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 22:34:31 2023

@author: user
"""

import xarray as xr
import numpy as np
import math
import matplotlib.pyplot as plt


startnum = 5.015
space_num = 4


k = np.linspace(startnum, space_num + startnum, space_num + 1)
print(k)



number = 3.7e-13
mantissa, exponent = math.frexp(number)
# Adjust the exponent to match the scientific notation
# exponent = exponent - 1

print(f"Mantissa: {mantissa}, Exponent: {exponent}")



# Data
categories = ['A', 'B', 'C', 'D']  # Categories
group1 = [10, 15, 7, 12]           # Values for group 1
group2 = [8, 12, 6, 10]            # Values for group 2
group3 = [5, 18, 9, 14]            # Values for group 3

# Number of categories
x = np.arange(len(categories))  # Positions for the categories

# Width of each bar
bar_width = 0.25

# Create the bars
plt.bar(x - bar_width, group1, width=bar_width, label='Group 1', color='blue')
plt.bar(x, group2, width=bar_width, label='Group 2', color='orange')
plt.bar(x + bar_width, group3, width=bar_width, label='Group 3', color='green')

# Add labels, title, and legend
plt.xlabel('Categories')
plt.ylabel('Values')
plt.title('Grouped Bar Chart')
plt.xticks(x, categories)  # Set the category labels at the center of the group
plt.legend()

# Show the plot
plt.tight_layout()
plt.show()




# Data
categories = ['A', 'B', 'C', 'D']  # Categories
values = [
    [10, 15, 7, 12],  # Group 1
    [8, 12, 6, 10],   # Group 2
    [5, 18, 9, 14]    # Group 3
]
group_labels = ['Group 1', 'Group 2', 'Group 3']
colors = ['blue', 'orange', 'green']

# Number of categories and groups
n_categories = len(categories)
n_groups = len(values)

# X positions for the categories
x = np.arange(n_categories)

# Bar width
bar_width = 0.2

# Create the plot
plt.figure(figsize=(10, 6))

# Loop through each group to plot
for i, group_values in enumerate(values):
    plt.bar(
        x + i * bar_width,  # Shift each group by bar_width
        group_values,       # Heights of the bars
        width=bar_width,    # Width of the bars
        label=group_labels[i],  # Label for the group
        color=colors[i]     # Color for the group
    )

# Add labels, title, and legend
plt.xlabel('Categories')
plt.ylabel('Values')
plt.title('Grouped Bar Chart')
plt.xticks(x + bar_width * (n_groups - 1) / 2, categories)  # Center the category labels
plt.legend()

# Show the plot
plt.tight_layout()
plt.show()




"""

Te_J = 1.66094e-16
ev = 1.6021766339999999 * pow(10, -19)
Te_data = Te_J / ev

print('Te = {}'.format(Te_data))


string = '80_cnf6.0_leakbsol_nts5_a'

if 'leakbsol_nts5_a' in string:
    
    print('we have tail recognize function')

else:
    
    print('we need to look for something else!')







a = np.zeros([2, 2, 1])
b = np.ones([2, 2])
a[:, :, 0] = b

pts = np.array([[2.1, 6.2, 8.3],
                [3.3, 5.2, 7.1]])


# k = lmem.mast_series_dir()


print(0.5*np.log(2+ np.sqrt(3)))
print(0.5*np.log(2- np.sqrt(3)))


c = [3, 50, 7]
c_total = [3, 4, 7, 9, 10, 19, 20, 21, 49, 50]

print(all([item in c_total for item in c]))


mod = [0, 0.3, 0.5, 0.7, 1]

for i in mod:
    A = (0.7 + i)/0.5
    print(A)

test_zero = np.zeros(10)
test_ones = np.ones(10)*0.5
print('try ones: {}'.format(test_ones))
test_id = np.identity(4, dtype=int)

test_value = np.all(test_id == 0)

print('the test result is: {}'.format(test_value))


test_dict = {"geeks": 7, "for": 1, "geek": 2}
dic_list = test_dict.keys()


for keys in test_dict:
    
    print(keys)



print('The data type of dic_list is:')

print(type(dic_list))


test_array_a = np.array([[1, 3, 5, 7], [2, 4, 6, 8]])

print(test_array_a)

test_array_b = np.array([[2, 4, 2, 8], [2.5, 4.5, 6.5, 8.5]])

print(test_array_b)

ratio_array = np.divide(test_array_a, test_array_b)

print(ratio_array)

"""

