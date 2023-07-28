# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 16:43:03 2023

@author: user
"""

if 'gas' in Shot:
    BASEDRT = '{}gaspuff/'.format(BASEDRT)


if AVG is True and 'AVG' not in Attempts:
    Attempts.append('AVG')
    print('Attempts {} Requested...'.format(Attempts))

        