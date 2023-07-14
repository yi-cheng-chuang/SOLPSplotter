# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 12:45:19 2023

@author: user
"""
import B2plotter_set as b2s
import B2plotter_class as b2c

d = b2s.Setting_dic()


# xd = b2c.B2plotter(DEV= d['DEV'], Shot= d['Shot'], shift= d['shift'], 
#                    series= d['series'], Attempts= d['Attempts'], 
#                    Parameters= d['Parameters'], 
#                    DefaultSettings= d['DefaultSettings'])


xp = b2c.mastdata(DEV= d['DEV'], Shot= d['Shot'], shift= d['shift'], 
                    series= d['series'], Attempts= d['Attempts'], 
                    Parameters= d['Parameters'], loadDS= d['loadDS'], 
                    DefaultSettings= d['DefaultSettings'])
# print(type(xd.Parameters))
xp.loadmastdata(a_shift='org')
xp.mastcalcpsi(a_shift='org', geo=None, b2mn=None, dsa=None, shift= 0)
xpp = xp.data['dirdic']['fitloc']
print(xpp)