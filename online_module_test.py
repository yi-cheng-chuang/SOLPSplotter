# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 20:45:25 2024

@author: user
"""

import os




aa = os.environ['OS']
print('The environment is {}'.format(aa))

# bb = os.environ['USERNAME']

# print('The username is {}'.format(bb))

bb = os.environ['USER']

print(f'The user is {bb}')


"""


for key, value in os.environ.items():
    print(f'{key}: {value}')
    

"""