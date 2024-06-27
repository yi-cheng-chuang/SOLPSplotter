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



