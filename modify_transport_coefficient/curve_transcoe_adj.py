# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 17:10:01 2025

@author: ychuang
"""

from matplotlib.offsetbox import AnchoredText
from modify_transport_coefficient.transport_coefficient_adjust_method import transcoe_method
from load_directory.grab_attempt_number import grab_aptn_method
from SOLPS_input.input_setting import set_figdir
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np



class curve_trancoe_adjustment:
    
    def __init__(self, DF, data, gam: grab_aptn_method, tm: transcoe_method):
        
        self.DF = DF
        self.data = data
        self.gam = gam
        self.tm = tm
    
    
    
    def fit(self):

        
        # Given points (x, y)
        x_points = np.array([1, 2, 3, 4, 5])
        y_points = np.array([2, 4, 9, 16, 25])
        
        # Fit a 2nd degree polynomial (parabola)
        coefficients = np.polyfit(x_points, y_points, 2)
        
        # coefficients will return [a, b, c] for ax^2 + bx + c
        a, b, c = coefficients
        
        print(f"The parabola equation is: y = {a}x^2 + {b}x + {c}")
        
        # Plotting the points and the fitted parabola
        x_fit = np.linspace(min(x_points), max(x_points), 100)
        y_fit = a * x_fit**2 + b * x_fit + c
        
        plt.scatter(x_points, y_points, color='red', label='Data points')
        plt.plot(x_fit, y_fit, label=f'Fitted Parabola: $y = {a:.2f}x^2 + {b:.2f}x + {c:.2f}$')
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Fitting a Parabola to Data Points')
        plt.show()

    
    
    
    
    def mod_transco_method(self,file_loc, withmod, de_SOL, ki_SOL, ke_SOL, log_flag):
        
        
        DEV = self.DF.DEV
        
        trans_list = self.tm.load_transcoefile_method(file_loc, plot= False)
        cod = trans_list['1'].T
        coki = trans_list['3'].T
        coke = trans_list['4'].T
        x= cod[:,0]  #the coordinate here is R-R_sep
        yd= cod[:,1]
        yki = coki[:,1]
        yke = coke[:,1]

        m = len(yd)
        if withmod:
            mod_y = np.zeros(m)
            for j in range(m):
                if j<= de_SOL:
                    mod_y[j] = cod[j,1]
                else:
                    mod_y[j] = 20.0
            cod[:,1] = mod_y

            mod_yki = np.zeros(m)
            for j in range(m):
                if j<= ki_SOL:
                    mod_yki[j] = coki[j,1]  
                else:
                    mod_yki[j] = 20.0
            coki[:,1] = mod_yki

            mod_yke = np.zeros(m)
            for j in range(m):
                if j<= ke_SOL:
                    mod_yke[j] = coke[j,1]  
                else:
                    mod_yke[j] = 20.0
            coke[:,1] = mod_yke
        else:
            pass


        self.tm.Generate_transcoefile_method(cod, CoeffID=1, SpeciesID=2, M=[1])
        
        shift = 'org'
        
        if DEV == 'mast':
            n = str(self.gam.s_number(file_loc)[0])
        
        elif DEV == 'mastu':
            
            n = str(self.gam.mastu_atp_number(file_loc, usage = 'transcoe')[0])
            

        simu_dir = file_loc.rsplit("/",2)[0]
        # k = str(ss.s_number(file_loc, series_flag= None))
        print(simu_dir)
        
        self.tm.Write_transcoefile_method(file = '{}/b2.transport.inputfile_mod_{}{}'.format(simu_dir, shift, n), points= trans_list ,M_1 = True, M=[1])


        specieslist = ['1','3','4']
        transcoe_unit = self.tm.transport_coe_unit()

        for k in specieslist:
            
            plt.figure()
            
            if log_flag:
                plt.yscale('log')
            else:
                pass
            plt.plot(trans_list[k][0,:], trans_list[k][1,:], 'o-', color = 'orange')
            plt.xlabel('Radial coordinate: $R- R_{sep}$')
            # plt.ylabel(transcoe_unit[k][1])
            plt.title(transcoe_unit[k][0])

        plt.show()