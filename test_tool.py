# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 20:19:42 2023

@author: user
"""

import os


def SET_WDIR(): #Function to set correct Working Directory Path depending on which machine is in use
    if os.environ['OS'] == 'Windows_NT':
        if os.environ['USERNAME'] == 'rmreksoatmodjo':
            BASEDRT = r"C:/Users/rmreksoatmodjo/Desktop/WMGDrive/College of William and Mary/Research/SOLPS Stuff/SOLPS_2D_prof/"
            TOPDRT = r"C:/Users/rmreksoatmodjo/Desktop/WMGDrive/College of William and Mary/Research/SOLPS Stuff/"
        elif os.environ['USERNAME'] == '18313':
            BASEDRT = r"G:/My Drive/College of William and Mary/Research/SOLPS Stuff/SOLPS_2D_prof/"
            TOPDRT = r"G:/My Drive/College of William and Mary/Research/SOLPS Stuff/"
        elif os.environ['USERNAME'] == 'Richard':
            BASEDRT = r"C:/Users/Richard/WMGDrive/College of William and Mary/Research/SOLPS Stuff/SOLPS_2D_prof/"
            TOPDRT = r"C:/Users/Richard/WMGDrive/College of William and Mary/Research/SOLPS Stuff/"
        elif os.environ['USERNAME'] == 'Yi-Cheng':
            BASEDRT = r"C:/Users/Yi-Cheng/Documents/SOLPS_Data/Simulation_Data"
            TOPDRT = r"C:/Users/Yi-Cheng/Documents/SOLPS_Data/Experimental_Data"
        elif os.environ['USERNAME'] == 'user':
            BASEDRT = r"C:/Users/user/Documents/SOLPS data/simulation data"
            TOPDRT = r"C:/Users/user/Documents/SOLPS data/experiment data"
    else:
        print("no user is listed")
    
    dirdic = {}
    dirdic['simulation directory'] = str(BASEDRT)
    dirdic['experiment directory'] = str(TOPDRT)
    
    return dirdic

if __name__ == '__main__':
    A = SET_WDIR()
    d1 = A['simulation directory']
    print(A)
    