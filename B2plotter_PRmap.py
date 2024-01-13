# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 17:27:01 2023

@author: user
"""

from SOLPSplotter_load_data import load_outputgen_data
import matplotlib.pyplot as plt
import load_coord_method as lcm
import fitting_method as fm 
from scipy.interpolate import interp1d
import numpy as np


class RP_mapping(load_outputgen_data):
    def __init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters):
        load_outputgen_data.__init__(self, DEV, withshift, withseries, DefaultSettings, loadDS, Parameters)
    
  
#-----------------------------------------------------------------------------

#calculate separatrix index, R-Rsep, R-Rsep psiN mapping  
  
    
    def calc_sep_index(self, psi, rad_range):
        index_low = []
        index_high = []
        index = []
        # print(psi)
        if psi[0] > psi[-1]:
            order = 'reverse_order_large2small'
        elif psi[0] < psi[-1]:
            order = 'order_small2large'
        
        if order == 'order_small2large':
            for y in range(len(psi)):
                if psi[y] <= 1:
                    index_low.append(y)
                elif psi[y] > 1:
                    index_high.append(y)
                else:
                    pass
        elif order == 'reverse_order_large2small':
            print('The input array is in reverse order')
        
        else:
            print('calc_sep_index function has a bug')
        
        index.append(index_low[-1])
        index.append(index_high[0])
        
        # print(type(index_high)) 
        # print('the following is index_low')
        # print(index_low)
        
        index_dic = {'index_low': index_low, 'index_high': index_high, 
                     'index': index}
        return index_dic
    
    
    
    def calc_RRsep_method(self, itername, plotRR):
        
        if itername == None:
        
            pol_range = int(self.data['b2fgeo']['nx'] + 2)
            rad_range = int(self.data['b2fgeo']['ny'] + 2)
            
            RadLoc = self.data['grid']['RadLoc']
            VertLoc = self.data['grid']['VertLoc']
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            
        elif itername != None:
            geo = self.data['b2fgeo'][itername]
            pol_range = int(geo['nx'] + 2)
            rad_range = int(geo['ny'] + 2)
            
            RadLoc = self.data['grid']['RadLoc'][itername]
            VertLoc = self.data['grid']['VertLoc'][itername]
            
        if self.withshift == False and self.withseries == False:
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            
        elif self.withshift == True and self.withseries == False:
            psiNinterp_RBS = self.data['gfile']['gcomp'][itername]['interp_dic']['RBS']
        
        elif self.withshift == False and self.withseries == True:
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            
        
        
        
        
        mag_axis_z = self.data['gfile']['g']['zmaxis']
        # print(mag_axis_z)
        
        
        "Calculate midplane R for Z=0"
        
        "Calculate weight"
        
        crup = RadLoc[:, 58]
        crlow = RadLoc[:, 60]
        czup = VertLoc[:, 58]
        czlow = VertLoc[:, 60]
        

        weight_mid = np.zeros(rad_range)
        for x in range(rad_range):
            weight_mid[x] = (mag_axis_z - czlow[x])/ (czup[x] - czlow[x])
        
        mid_choice = np.zeros((rad_range, 4))
        mid_choice[:, 0] = czup
        mid_choice[:, 1] = czlow
        mid_choice[:, 2] = crup
        mid_choice[:, 3] = crlow
        
        
        
        mid_R = np.zeros(rad_range)
        for xa in range(rad_range):
            mid_R[xa] = weight_mid[xa]*crup[xa] + (1 - weight_mid[xa])*crlow[xa]
        
        mid_Z = np.zeros(rad_range)
        for xb in range(rad_range):
            mid_Z[xb] = weight_mid[xb]*czup[xb] + (1 - weight_mid[xb])*czlow[xb]
            
            
        psi_solps_mid = np.zeros(rad_range)
        
        # psi_solps_cp = psiNinterp_RBS(crloc, czloc)
        for i in range(rad_range):
            psi_mid = psiNinterp_RBS(mid_R[i], mid_Z[i])
            psi_solps_mid[i] = psi_mid
        
        
        sep_index_dic = self.calc_sep_index(psi = psi_solps_mid, rad_range = rad_range)
        sep_index_high = int(sep_index_dic['index'][1])
        
        
        weight_psi = (psi_solps_mid[sep_index_high] - 1)/(psi_solps_mid[sep_index_high] -psi_solps_mid[sep_index_high -1])
        
        R_sep = (1 - weight_psi)*mid_R[sep_index_high] + weight_psi*mid_R[sep_index_high -1]
        R_Rsep = mid_R - R_sep
        
        midplane_dic = {'weight': weight_mid, 'mid_choice': mid_choice, 
                        'mid_R': mid_R, 'mid_Z': mid_Z, 
                        'psi_solps_mid': psi_solps_mid, 
                        'weight_psi': weight_psi, 'R_Rsep': R_Rsep}
        
        psi_dsa_dic = fm.dsa_psi_fit(dsa= R_Rsep, psi= psi_solps_mid)
        
        psi_dsa_ratio = psi_dsa_dic['dsa_psi_fitcoe'][0]
        
        
        pol_list = [57, 58, 59, 60]
        if plotRR:
            plt.figure(figsize=(7,7))
            for in_pol in pol_list:
                crloc = RadLoc[:, int(in_pol)]
                czloc = VertLoc[:, int(in_pol)]
                plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
            plt.plot(mid_R, mid_Z, color = 'r', label= 'R_Rsep')
            plt.xlabel('R: [m]')
            plt.ylabel('Z: [m]')
    
    
        return midplane_dic, psi_dsa_ratio, sep_index_high
    
    
    
    def calc_RRsep(self, plotRR, plot_psi_dsa_align):
        
        if self.withshift == False and self.withseries == False:
            midplane_dic, psi_dsa_ratio, sep_index_high = self.calc_RRsep_method(itername = None, plotRR = plotRR)
            
            self.data['DefaultSettings']['psi_dsa'] = psi_dsa_ratio
            self.data['midplane_calc'] = midplane_dic
            self.data['DefaultSettings']['SEP'] = sep_index_high
        
        
        elif self.withshift == True and self.withseries == False:        
            
            sep_index_dic = {}
            midplane_series_dic = {}
            psi_dsa_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                midplane_dic, psi_dsa_ratio, sep_index_high = self.calc_RRsep_method(itername = aa, plotRR = plotRR)
                
                midplane_series_dic[aa] = midplane_dic
                psi_dsa_dic[aa] = psi_dsa_ratio
                sep_index_dic[aa] = sep_index_high
            
            self.data['DefaultSettings']['psi_dsa'] = psi_dsa_dic
            self.data['midplane_calc'] = midplane_series_dic
            self.data['DefaultSettings']['SEP'] = sep_index_dic
                
                
            if plot_psi_dsa_align:
                color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                             'dot7': 'blue', 'one': 'purple'}
                plt.figure(figsize=(7,7))
                for in_sh in self.data['dircomp']['multi_shift']:
                    # plt.plot(midplane_dic[in_sh]['psi_solps_mid'], midplane_dic[in_sh]['R_Rsep'], 'o-', 
                    #   label= 'R-Rsep:[m] with modify {} m'.format(self.data['dircomp']['shift_dic'][in_sh]))
                    plt.plot(midplane_dic[in_sh]['psi_solps_mid'], midplane_dic[in_sh]['R_Rsep'], 'o-', 
                      color= color_dic[in_sh])
                plt.xlabel('psiN')
                plt.title('psiN R-Rsep plot')
                # plt.legend()
                    
           
        elif self.withshift == False and self.withseries == True:
            den_rep = list(self.data['dircomp']['Attempt'].keys())[0]
            
            midplane_dic, psi_dsa_ratio, sep_index_high = self.calc_RRsep_method(self, itername = den_rep, plotRR = plotRR)
            
            self.data['DefaultSettings']['psi_dsa'] = psi_dsa_ratio
            self.data['midplane_calc'] = midplane_dic
            self.data['DefaultSettings']['SEP'] = sep_index_high
        
        elif self.withshift == True and self.withseries == True:
            print('calc_RRsep function is not there yet!')
        
        else:
            print('calc_RRsep function has a bug')
            
            
#-----------------------------------------------------------------------------

#calculate the poloidal angle function    
        
    def calc_pol_angle_method(self, itername, pol_list, plot_angle):
        
        if itername == None:
        
            pol_range = int(self.data['b2fgeo']['nx'] + 2)
            rad_range = int(self.data['b2fgeo']['ny'] + 2)
            
            RadLoc = self.data['grid']['RadLoc']
            VertLoc = self.data['grid']['VertLoc']
            
        
        elif itername != None:
            
            geo = self.data['b2fgeo'][itername] 
            pol_range = int(geo['nx'] + 2)
            rad_range = int(geo['ny'] + 2)
                
            RadLoc = self.data['grid']['RadLoc'][itername]
            VertLoc = self.data['grid']['VertLoc'][itername]
        
               
        gfile = self.data['gfile']['g']
        mag_axis_z = gfile['zmaxis']
        # print(mag_axis_z)
        
        
        if self.withshift == False and self.withseries == False:
            mag_axis_r = self.data['gfile']['g']['rmaxis'] + self.data['dircomp']['shift_value']
            # print(mag_axis_r)
        
        elif self.withshift == True and self.withseries == False:
            mag_axis_r = self.data['gfile']['g']['rmaxis'] + self.data['dircomp']['shift_dic'][itername]
            # print(mag_axis_r)
            
        elif self.withshift == False and self.withseries == True:
            mag_axis_r = self.data['gfile']['g']['rmaxis'] + self.data['dircomp']['shift_value']
            # print(mag_axis_r)
            
        elif self.withshift == True and self.withseries == True:
            print('calc_pol_angle_method function is not there yet!')
        
        else:
            print('calc_pol_angle_method function has a bug')
            
               
        xprloc = RadLoc[:, 72][19]
        xpzloc = VertLoc[:, 72][19]
        
        xpoint_R = xprloc - mag_axis_r
        xpoint_Z = xpzloc - mag_axis_z
        xpoint_angle = np.arctan2(xpoint_Z, xpoint_R)*180 / np.pi
        if self.withshift == False:
            print('xpoint angle is {}'.format(str(xpoint_angle)))
        elif self.withshift:
            print('xpoint angle is {:.2f} for {} case'.format(xpoint_angle, itername))
        else:
            print('xpoint angle has a bug')
            
            
        angle_list = []
        
        for pol_index in pol_list:
            rloc = RadLoc[:, int(pol_index)][19]
            zloc = VertLoc[:, int(pol_index)][19]
            
            x = rloc - mag_axis_r
            y = zloc - mag_axis_z
            
            
            angle = np.arctan2(y, x)*180 / np.pi
            if angle <= xpoint_angle:
                angle = angle + 360
            angle_list.append(angle)
            
        pol_loc = np.asarray(int(pol_list))
        
        if plot_angle: 
            plt.figure(figsize=(7,7))
            plt.plot(pol_loc, angle_list, 'o', color = 'r', label= 'poloidal angle')
            plt.xlabel('poloidal index')
            plt.title('poloidal angle verses poloidal index')
            plt.legend()
                
        return angle_list, xpoint_angle
    
    
    
    def calc_pol_angle(self, pol_list, plot_angle):
        
        if self.withshift == False and self.withseries == False:
            angle_list, xpoint_angle = self.calc_pol_angle_method(itername = None, 
                                        plot_angle = plot_angle, pol_list = pol_list)
            
            angle_dic = {'angle_list': angle_list, 'xpoint_angle': xpoint_angle}
            self.data['angle'] = angle_dic
            
        elif self.withshift == True and self.withseries == False:        
            angle_dic = {}
            angle_list_dic = {}
            xpoint_angle_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                angle_list, xpoint_angle = self.calc_pol_angle_method(itername = aa, 
                                            plot_angle = plot_angle, pol_list = pol_list)
                
                angle_list_dic[aa] = angle_list
                xpoint_angle_dic[aa] = xpoint_angle
                
    
            angle_dic = {'angle_list': angle_list_dic, 'xpoint_angle': xpoint_angle_dic}
            self.data['angle'] = angle_dic
            
           
        elif self.withshift == False and self.withseries == True:
            
            den_list = list(self.data['dircomp']['Attempt'].keys())
            aa = den_list[0]
            
            angle_list, xpoint_angle = self.calc_pol_angle_method(itername = aa, 
                                        plot_angle = plot_angle, pol_list = pol_list)
            
            angle_dic = {'angle_list': angle_list, 'xpoint_angle': xpoint_angle}
            self.data['angle'] = angle_dic
        
        elif self.withshift == True and self.withseries == True:
            print('calc_pol_angle function is not there yet!')
        
        else:
            print('calc_pol_angle function has a bug')
            
            
            
#-----------------------------------------------------------------------------

#calculate the arc length of each poloidal index curve
    
    """
    from stackoverflow
    How to interpolate a 2D curve in Python
    https://stackoverflow.com/questions/52014197/how-to-interpolate-a-2d-curve-in-python


    """

    def RR_sep_calculator(self, cr, cz):
        # Define some points:
        p = len(cr)
        points = np.zeros((p, 2))
        points[:, 0] = cr
        points[:, 1] = cz
        # points = np.array([[0, 1, 8, 2, 2],
        #                    [1, 0, 6, 7, 2]]).T  # a (nbre_points x nbre_dim) array
        
        "Linear length along the line:"
        
        # crzdiff = np.diff(points, axis=0)
        # crzsum = np.sum( np.diff(points, axis=0)**2, axis=1 )
        distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )))
        distance = np.insert(distance, 0, 0)/distance[-1]
        
        arclength = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )))
        arclength = np.insert(arclength, 0, 0)
        
        
        "Interpolation for different methods:"
        
        interpolations_methods = ['slinear', 'quadratic', 'cubic']
        alpha = np.linspace(0, 1, 38)
        
        interpolated_points = {}
        for method in interpolations_methods:
            interpolator =  interp1d(distance, points, kind=method, axis=0)
            interpolated_points[method] = interpolator(alpha)
        
        "Graph:"
        # plt.figure(figsize=(7,7))
        # for method_name, curve in interpolated_points.items():
        #     plt.plot(*curve.T, '-', label=method_name);
        
        # plt.plot(*points.T, 'ok', label='original points');
        # plt.axis('equal'); plt.legend(); plt.xlabel('r'); plt.ylabel('z');

        return arclength, interpolated_points       

    

    def calc_dsa_method(self, itername, pol_loc):
        
        if itername == None:
            RadLoc = self.data['grid']['RadLoc']
            VertLoc = self.data['grid']['VertLoc']
            
        elif itername != None:
            RadLoc = self.data['grid']['RadLoc'][itername]
            VertLoc = self.data['grid']['VertLoc'][itername]
        else:
            print('calc_dsa_method function has a bug')
            
            
        if self.withshift == False:
            SEP = self.data['DefaultSettings']['SEP']
            
        elif self.withshift == True and self.withseries == False:
            SEP = self.data['DefaultSettings']['SEP'][itername]
            
        elif self.withshift == True and self.withseries == True:
            print('calc_dsa_method function is not there yet!')
        
        else:
            print('calc_dsa_method has a bug')
        
        
        pol_index = int(pol_loc)
        crloc = RadLoc[:, pol_index]
        czloc = VertLoc[:, pol_index]
        # print(crloc)
        
        arclength, interpfunc_dic = self.RR_sep_calculator(cr = crloc, cz = czloc)
        
        sep_dist = np.mean([arclength[int(SEP)+ 1], arclength[int(SEP)]])
        # print(sep_dist)
                
        RR_sep = arclength - sep_dist
        
        return arclength, interpfunc_dic, RR_sep
        
    
    
    
   
    def calc_dsa(self, pol_loc):
        dsa = lcm.read_dsa(self.data['dirdata']['simudir'] + '/dsa')
        if self.withshift == False and self.withseries == False:
            
            arclength, interpfunc_dic, RR_sep = self.calc_dsa_method(itername = None, pol_loc = pol_loc)
            
            dsa_dic = {'arclength': arclength, 'interpfunc': interpfunc_dic,
                        'dsa_{}_val'.format(pol_loc): RR_sep}
            
            self.data['dsa']['dsa_{}'.format(pol_loc)] = dsa_dic
            
        
        elif self.withshift == True and self.withseries == False:
            dsa_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                
            
                dsa_dic[aa] = {'arclength': arclength, 'interpfunc': interpfunc_dic,
                            'dsa_{}_val'.format(pol_index): RR_sep}
            
            self.data['dsa']['dsa_{}'.format(pol_index)] = dsa_dic
    
    
    
        elif self.withshift == False and self.withseries == True:
            
            den_list = list(self.data['dircomp']['Attempt'].keys())
            aa = den_list[0] 
            
            # pol_range = int(self.data['b2fgeo']['nx'] + 2)
            rad_range = int(self.data['b2fgeo']['ny'] + 2)
            
            
                       
            Attempt = self.data['dircomp']['Attempt'][aa]
            DRT = self.data['dirdata']['outputdir'][aa]['Output']
            XDIM = self.data['b2fgeo']['nx'] + 2
            YDIM = self.data['b2fgeo']['ny'] + 2
            # dsa = lcm.read_dsa(self.data['dirdata']['simudir'] + '/dsa')
            
            RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                        usecols = (3)).reshape((YDIM, XDIM))
            VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                          usecols = (3)).reshape((YDIM,XDIM))
            
            pol_index = int(pol_loc)
            crloc = RadLoc[:, pol_index]
            czloc = VertLoc[:, pol_index]
            # print(crloc)
            
            arclength, interpfunc_dic = self.RR_sep_calculator(cr = crloc, cz = czloc)
            
            SEP = self.data['DefaultSettings']['SEP']
            
            sep_dist = np.mean([arclength[int(SEP)+ 1], arclength[int(SEP)]])
            # print(sep_dist)
            
            
            RR_sep = arclength - sep_dist
            
            dsa_dic = {'arclength': arclength, 'interpfunc': interpfunc_dic,
                        'dsa_{}_val'.format(pol_index): RR_sep}
            
            self.data['dsa']['dsa_{}'.format(pol_index)] = dsa_dic


#-----------------------------------------------------------------------------
      
# Additional dsa check function 1: left_test           
        
    def left_dsa_test_method(self, itername, pol_index, sep_loc):
        "Test for left only result"
        
        
        if self.withshift == False:
            geo = self.data['b2fgeo']
            rad_range = int(geo['ny'] + 2)
            
        elif self.withshift == True and self.withseries == False:
            geo = self.data['b2fgeo'][itername]
            rad_range = int(geo['ny'] + 2)
        
        elif self.withshift == True and self.withseries == True:
            print('left_dsa_test function is not there yet!')
            
        else:
            print('left_dsa_test function has a bug')
            
            
        crLowerLeft = geo['crx'][pol_index,:,0]
        crLowerRight = geo['crx'][pol_index,:,1]
        crUpperLeft = geo['crx'][pol_index,:,2]
        crUpperRight = geo['crx'][pol_index,:,3]
        czLowerLeft = geo['cry'][pol_index,:,0]
        czLowerRight = geo['cry'][pol_index,:,1]
        czUpperLeft = geo['cry'][pol_index,:,2]
        czUpperRight = geo['cry'][pol_index,:,3]
        
        SEP = sep_loc
        
        LLsep = np.mean([crLowerLeft[int(SEP)-1], crLowerLeft[int(SEP)-2]])
        LRsep = np.mean([crLowerRight[int(SEP)-1], crLowerRight[int(SEP)-2]])
        ULsep = np.mean([crUpperLeft[int(SEP)-1], crUpperLeft[int(SEP)-2]])
        URsep = np.mean([crUpperRight[int(SEP)-1], crUpperRight[int(SEP)-2]])
        
        
        LLdsa = crLowerLeft - LLsep
        LRdsa = crLowerRight - LRsep
        ULdsa = crUpperLeft - ULsep
        URdsa = crUpperRight - URsep
        avag_rad = np.zeros(rad_range)
        for j in range(rad_range):
            avag_rad[j] = np.mean([LLdsa[j], ULdsa[j], URdsa[j], LRdsa[j]])
      
        crzsum, points = fm.RR_sep_calculator(cr = crLowerLeft, cz = czLowerLeft)
        
               
        if self.withshift == False and self.withseries == False:
            Attempt = self.data['dircomp']['Attempt']
            DRT2 = self.data['dirdata']['outputdir']['Output2']
                  
        elif self.withshift == True and self.withseries == False:
            Attempt = self.data['dircomp']['Attempt'][itername]
            DRT2 = self.data['dirdata']['infolderdir'][itername]['outputdir']['Output2']
        
        elif self.withshift == False and self.withseries == True:
            Attempt = self.data['dircomp']['Attempt'][itername]
            DRT2 = self.data['dirdata']['outputdir'][itername]['Output2']
        
        elif self.withshift == True and self.withseries == True:
            print('left_dsa_test function is not there yet!')
            
        else:
            print('left_dsa_test function has a bug')
            
            
        XDIM = geo['nx'] + 2
        YDIM = geo['ny'] + 2   
        
        Rad0Cor = np.loadtxt('{}/Rad0Cor{}'.format(DRT2, str(Attempt)),
                    usecols = (3)).reshape((YDIM, XDIM))
        Vert0Cor = np.loadtxt('{}/Vert0Cor{}'.format(DRT2, str(Attempt)), 
                      usecols = (3)).reshape((YDIM,XDIM))
        
        compare_left = np.zeros((int(rad_range), 2))
        compare_left[:, 0] = Rad0Cor[:, pol_index]
        compare_left[:, 1] = crLowerLeft
        
        return compare_left
    
        
    def left_test(self, pol_loc):
        pol_pos = int(pol_loc)
        
        if self.withshift == False and self.withseries == False:
            SEP = self.data['DefaultSettings']['SEP']
            compare_left = self.left_dsa_test(itername = None, pol_index = pol_pos, 
                               sep_loc = SEP)
            self.data['dsa_compare_left'] = compare_left
        
        elif self.withshift == True and self.withseries == False:
            compare_left_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                SEP = self.data['DefaultSettings']['SEP'][aa]
                compare_left = self.left_dsa_test(itername = aa, pol_index = pol_pos, 
                                   sep_loc = SEP)
                
                compare_left_dic[aa] = compare_left
                
                self.data['dsa_compare_left'] = compare_left_dic
        
        elif self.withshift == False and self.withseries == True:
            den_list = list(self.data['dircomp']['Attempt'].keys())
            aa = den_list[0]
            
            SEP = self.data['DefaultSettings']['SEP']
            compare_left = self.left_dsa_test(itername = aa, pol_index = pol_pos, 
                               sep_loc = SEP)
            
            self.data['dsa_compare_left'] = compare_left
            
        
        elif self.withshift == True and self.withseries == True:
            print('left_test function is not there yet!')
        
        else:
            print('left_test function has a bug')
            
        
# Additional dsa check function 2: dsa_and_arclength_check

    def dsa_and_arclength_check(self):
        dsa = lcm.read_dsa(self.data['dirdata']['simudir'] + '/dsa')






"""
backup

# sep = []
# for in_sep in psi_solps_mid:
#     if in_sep >= 1:
#         sep.append(list(psi_solps_mid).index(in_sep))

        
# sep_index = sep[0]
# # print(sep)

pol_range = int(self.data['b2fgeo']['nx'] + 2)
rad_range = int(self.data['b2fgeo']['ny'] + 2)
mag_axis_z = self.data['gfile']['g']['zmaxis']
# print(mag_axis_z)

# Attempt = self.data['dircomp']['Attempt']
# DRT = self.data['dirdata']['outputdir']['Output']
# DRT2 = self.data['dirdata']['outputdir']['Output2']
# XDIM = self.data['b2fgeo']['nx'] + 2
# YDIM = self.data['b2fgeo']['ny'] + 2


RadLoc = self.data['grid']['RadLoc']
VertLoc = self.data['grid']['VertLoc']


"Calculate midplane R for Z=0"

"Calculate weight"

crup = RadLoc[:, 58]
crlow = RadLoc[:, 60]
czup = VertLoc[:, 58]
czlow = VertLoc[:, 60]

pol_list = [52, 53, 54, 55, 56, 57]


weight_mid = np.zeros(rad_range)
for x in range(rad_range):
    weight_mid[x] = (mag_axis_z - czlow[x])/ (czup[x] - czlow[x])

mid_choice = np.zeros((rad_range, 4))
mid_choice[:, 0] = czup
mid_choice[:, 1] = czlow
mid_choice[:, 2] = crup
mid_choice[:, 3] = crlow



mid_R = np.zeros(rad_range)
for xa in range(rad_range):
    mid_R[xa] = weight_mid[xa]*crup[xa] + (1 - weight_mid[xa])*crlow[xa]

mid_Z = np.zeros(rad_range)
for xb in range(rad_range):
    mid_Z[xb] = weight_mid[xb]*czup[xb] + (1 - weight_mid[xb])*czlow[xb]
    
pol_list = [57, 58, 59, 60]



psi_solps_mid = np.zeros(rad_range)
psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
# psi_solps_cp = psiNinterp_RBS(crloc, czloc)
for i in range(rad_range):
    psi_mid = psiNinterp_RBS(mid_R[i], mid_Z[i])
    psi_solps_mid[i] = psi_mid

sep = []
for in_sep in psi_solps_mid:
    if in_sep >= 1:
        sep.append(list(psi_solps_mid).index(in_sep))

        
sep_index = sep[0]
# print(sep)

weight_psi = (psi_solps_mid[sep_index] - 1)/(psi_solps_mid[sep_index] -psi_solps_mid[sep_index -1])

R_sep = (1 - weight_psi)*mid_R[sep_index] + weight_psi*mid_R[sep_index -1]
R_Rsep = mid_R - R_sep

midplane_dic = {'weight': weight_mid, 'mid_choice': mid_choice, 
                'mid_R': mid_R, 'mid_Z': mid_Z, 
                'psi_solps_mid': psi_solps_mid, 
                'weight_psi': weight_psi, 'R_Rsep': R_Rsep}

psi_dsa_dic = fm.dsa_psi_fit(dsa= R_Rsep, psi= psi_solps_mid)

self.data['DefaultSettings']['psi_dsa'] = psi_dsa_dic['dsa_psi_fitcoe'][0]
            
self.data['midplane_calc'] = midplane_dic

if plotRR:
    plt.figure(figsize=(7,7))
    for in_pol in pol_list:
        crloc = RadLoc[:, int(in_pol)]
        czloc = VertLoc[:, int(in_pol)]
        plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
    plt.plot(mid_R, mid_Z, color = 'r', label= 'R_Rsep')
    plt.xlabel('R: [m]')
    plt.ylabel('Z: [m]')
    

pol_range = int(self.data['b2fgeo'][aa]['nx'] + 2)
# print('xdim is {}'.format(str(pol_range)))
rad_range = int(self.data['b2fgeo'][aa]['ny'] + 2)
# print('ydim is {}'.format(str(rad_range)))
pol_range_dic[aa] = pol_range
rad_range_dic[aa] = rad_range

mag_axis_z = self.data['gfile']['g']['zmaxis']
# print(mag_axis_z)

Attempt = self.data['dircomp']['Attempt'][aa]
DRT = self.data['dirdata']['infolderdir'][aa]['outputdir']['Output']
DRT2 = self.data['dirdata']['infolderdir'][aa]['outputdir']['Output2']


RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
            usecols = (3)).reshape((rad_range, pol_range))
VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
              usecols = (3)).reshape((rad_range, pol_range))


"Calculate midplane R for Z=0"

"Calculate weight"

if aa == 'one':
    crup = RadLoc[:, 57]
    crlow = RadLoc[:, 59]
    czup = VertLoc[:, 57]
    czlow = VertLoc[:, 59]
else:
    crup = RadLoc[:, 57]
    crlow = RadLoc[:, 59]
    czup = VertLoc[:, 57]
    czlow = VertLoc[:, 59]


weight_mid = np.zeros(rad_range)
for x in range(rad_range):
    weight_mid[x] = (czup[x] - mag_axis_z)/ (czup[x] - czlow[x])

mid_choice = np.zeros((rad_range, 4))
mid_choice[:, 0] = czup
mid_choice[:, 1] = czlow
mid_choice[:, 2] = crup
mid_choice[:, 3] = crlow



mid_R = np.zeros(rad_range)
for xa in range(rad_range):
    mid_R[xa] = (1 - weight_mid[xa])*crup[xa] + weight_mid[xa]*crlow[xa]

mid_Z = np.zeros(rad_range)
for xb in range(rad_range):
    mid_Z[xb] = (1 - weight_mid[xb])*czup[xb] + weight_mid[xb]*czlow[xb]
    

pol_list = [57, 58, 59, 60]


psi_solps_mid = np.zeros(rad_range)
psiNinterp_RBS = self.data['gfile']['gcomp'][aa]['interp_dic']['RBS']
# psi_solps_cp = psiNinterp_RBS(crloc, czloc)
for i in range(rad_range):
    psi_mid = psiNinterp_RBS(mid_R[i], mid_Z[i])
    psi_solps_mid[i] = psi_mid

sep = []
for in_sep in psi_solps_mid:
    if in_sep >= 1:
        sep.append(list(psi_solps_mid).index(in_sep))

        
sep_index = sep[0]
SEP_dic[aa] = sep_index

weight_psi = (psi_solps_mid[sep_index] - 1)/(psi_solps_mid[sep_index] -psi_solps_mid[sep_index -1])

R_sep = (1 - weight_psi)*mid_R[sep_index] + weight_psi*mid_R[sep_index -1]
R_Rsep = mid_R - R_sep

midplane_dic[aa] = {'weight': weight_mid, 'mid_choice': mid_choice, 
                'mid_R': mid_R, 'mid_Z': mid_Z, 
                'psi_solps_mid': psi_solps_mid, 
                'weight_psi': weight_psi, 'R_Rsep': R_Rsep}

pd_dic = fm.dsa_psi_fit(dsa= R_Rsep, psi= psi_solps_mid)

psi_dsa_dic[aa] = pd_dic['dsa_psi_fitcoe'][0]



if plotRR:
    plt.figure(figsize=(7,7))
    for in_pol in pol_list:
        crloc = RadLoc[:, int(in_pol)]
        czloc = VertLoc[:, int(in_pol)]
        plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
    plt.plot(mid_R, mid_Z, color = 'r', label= 'R_Rsep')
    plt.xlabel('R: [m]')
    plt.ylabel('Z: [m]')
    plt.title('{} R-Z plot'.format(aa))
    
            self.data['DefaultSettings']['psi_dsa'] = psi_dsa_dic
            self.data['midplane_calc'] = midplane_dic
            self.data['DefaultSettings']['XDIM'] = pol_range_dic
            self.data['DefaultSettings']['YDIM'] = rad_range_dic
            self.data['DefaultSettings']['SEP'] = SEP_dic
            
            
     pol_range = self.data['b2fgeo']['nx'] + 2
     rad_range = self.data['b2fgeo']['ny'] + 2
     
     
     solps_dsa_dic = {}
     for aa in self.data['dircomp']['Attempt'].keys():
         solps_dsa_dic[aa] = lcm.read_dsa(self.data['dirdata']['simudir'][aa] + '/dsa')
     
     den_list = list(self.data['dircomp']['Attempt'].keys())
     aa = den_list[0]
     
     mag_axis_z = self.data['gfile']['g']['zmaxis']
     print(mag_axis_z)
     
     Attempt = self.data['dircomp']['Attempt'][aa]
     DRT = self.data['dirdata']['outputdir'][aa]['Output']
     # DRT2 = self.data['dirdata']['outputdir']['Output2']
     # XDIM = self.data['b2fgeo']['nx'] + 2
     # YDIM = self.data['b2fgeo']['ny'] + 2
     
     
     RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                 usecols = (3)).reshape((rad_range, pol_range))
     VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                   usecols = (3)).reshape((rad_range, pol_range))
     
     
     "Calculate midplane R for Z=0"
     
     "Calculate weight"
     crup = RadLoc[:, 58]
     crlow = RadLoc[:, 60]
     czup = VertLoc[:, 58]
     czlow = VertLoc[:, 60]
     
     weight_mid = np.zeros(rad_range)
     for x in range(rad_range):
         weight_mid[x] = (czup[x] - mag_axis_z)/ (czup[x] - czlow[x])
     
     mid_choice = np.zeros((rad_range, 4))
     mid_choice[:, 0] = czup
     mid_choice[:, 1] = czlow
     mid_choice[:, 2] = crup
     mid_choice[:, 3] = crlow
     
     
     
     mid_R = np.zeros(rad_range)
     for xa in range(rad_range):
         mid_R[xa] = (1 - weight_mid[xa])*crup[xa] + weight_mid[xa]*crlow[xa]
     
     mid_Z = np.zeros(rad_range)
     for xb in range(rad_range):
         mid_Z[xb] = (1 - weight_mid[xb])*czup[xb] + weight_mid[xb]*czlow[xb]
         
     
     pol_list = [57, 58, 59, 60]
     
     
     psi_solps_mid = np.zeros(rad_range)
     psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
     # psi_solps_cp = psiNinterp_RBS(crloc, czloc)
     for i in range(rad_range):
         psi_mid = psiNinterp_RBS(mid_R[i], mid_Z[i])
         psi_solps_mid[i] = psi_mid
     
     sep = []
     for in_sep in psi_solps_mid:
         if in_sep >= 1:
             sep.append(list(psi_solps_mid).index(in_sep))
     
             
     sep_index = sep[0]
     # print(sep)
     
     weight_psi = (psi_solps_mid[sep_index] - 1)/(psi_solps_mid[sep_index] -psi_solps_mid[sep_index -1])
     
     R_sep = (1 - weight_psi)*mid_R[sep_index] + weight_psi*mid_R[sep_index -1]
     R_Rsep = mid_R - R_sep
     
     psi_dsa_dic = fm.dsa_psi_fit(dsa= R_Rsep, psi= psi_solps_mid)           
     self.data['DefaultSettings']['psi_dsa'] = psi_dsa_dic['dsa_psi_fitcoe'][0]
     
     midplane_dic = {'weight': weight_mid, 'mid_choice': mid_choice, 
                     'mid_R': mid_R, 'mid_Z': mid_Z, 
                     'psi_solps_mid': psi_solps_mid, 
                     'weight_psi': weight_psi, 'R_Rsep': R_Rsep}
         
     self.data['midplane_calc'] = midplane_dic
     
     if plotRR:
         plt.figure(figsize=(7,7))
         for in_pol in pol_list:
             crloc = RadLoc[:, int(in_pol)]
             czloc = VertLoc[:, int(in_pol)]
             plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
         plt.plot(mid_R, mid_Z, color = 'r', label= 'R_Rsep')
         plt.xlabel('R: [m]')
         plt.ylabel('Z: [m]')
         plt.title('{} R-Z plot'.format(aa))       


          
Attempt = self.data['dircomp']['Attempt']
DRT = self.data['dirdata']['outputdir']['Output']
XDIM = self.data['b2fgeo']['nx'] + 2
YDIM = self.data['b2fgeo']['ny'] + 2


RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
            usecols = (3)).reshape((YDIM, XDIM))
VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
              usecols = (3)).reshape((YDIM,XDIM))

mag_axis_z = self.data['gfile']['g']['zmaxis']
print(mag_axis_z)
mag_axis_r = self.data['gfile']['g']['rmaxis'] + self.data['dircomp']['shift_value']
print(mag_axis_r)

xprloc = RadLoc[:, 72][19]
xpzloc = VertLoc[:, 72][19]

xpoint_R = xprloc - mag_axis_r
xpoint_Z = xpzloc - mag_axis_z
xpoint_angle = np.arctan2(xpoint_Z, xpoint_R)*180 / np.pi
print('xpoint angle is {}'.format(str(xpoint_angle)))

         
angle_list = []

for pol_index in pol_list:
    rloc = RadLoc[:, int(pol_index)][19]
    zloc = VertLoc[:, int(pol_index)][19]
    
    x = rloc - mag_axis_r
    y = zloc - mag_axis_z
    
    
    angle = np.arctan2(y, x)*180 / np.pi
    if angle <= xpoint_angle:
        angle = angle + 360
    angle_list.append(angle)
    
ln = len(pol_list)
pol_loc = np.zeros(ln)
i = 0

for ii in pol_list:
    pol_loc[i] = int(ii)
    i = i + 1

if plot_angle: 
    plt.figure(figsize=(7,7))
    plt.plot(pol_loc, angle_list, 'o', color = 'r', label= 'poloidal angle')
    plt.xlabel('poloidal index')
    plt.title('poloidal angle verses poloidal index')
    plt.legend()

angle_dic = {'angle_list': angle_list, 'xpoint_angle': xpoint_angle}

pol_range = int(self.data['b2fgeo'][aa]['nx'] + 2)
    # print('xdim is {}'.format(str(pol_range)))
    rad_range = int(self.data['b2fgeo'][aa]['ny'] + 2)
    # print('ydim is {}'.format(str(rad_range)))
    
    mag_axis_z = self.data['gfile']['g']['zmaxis']
    mag_axis_r = self.data['gfile']['g']['rmaxis'] + self.data['dircomp']['shift_dic'][aa]
    # print(mag_axis_z)
    
    Attempt = self.data['dircomp']['Attempt'][aa]
    DRT = self.data['dirdata']['infolderdir'][aa]['outputdir']['Output']
    
              
    RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                usecols = (3)).reshape((rad_range, pol_range))
    VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                  usecols = (3)).reshape((rad_range, pol_range))
    
    xprloc = RadLoc[:, 72][19]
    xpzloc = VertLoc[:, 72][19]
    
    xpoint_R = xprloc - mag_axis_r
    xpoint_Z = xpzloc - mag_axis_z
    xpoint_angle = np.arctan2(xpoint_Z, xpoint_R)*180 / np.pi
    print('xpoint angle is {:.2f} for {} case'.format(xpoint_angle, aa))
    
    xpoint_angle_dic[aa] = xpoint_angle
    
    angle_list = []
    
    for pol_index in pol_list:
        rloc = RadLoc[:, int(pol_index)][19]
        zloc = VertLoc[:, int(pol_index)][19]
        
        x = rloc - mag_axis_r
        y = zloc - mag_axis_z
        
        
        angle = np.arctan2(y, x)*180 / np.pi
        if angle <= xpoint_angle:
            angle = angle + 360
        angle_list.append(angle)
        
        
                      
    angle_list_dic[aa] = angle_list
        
ln = len(pol_list)
pol_loc = np.zeros(ln)
i = 0

for ii in pol_list:
    pol_loc[i] = int(ii)
    i = i + 1

plt.figure(figsize=(7,7))
color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
             'dot7': 'blue', 'one': 'purple'}
for ab in self.data['dircomp']['multi_shift']: 
    plt.plot(pol_loc, angle_list_dic[ab], 'o', color = color_dic[ab])
    plt.xlabel('poloidal index')
    plt.title('poloidal angle verses poloidal index from {} to {}'.format(int(min(pol_loc)), int(max(pol_loc))))
    
angle_dic = {'angle_list': angle_list_dic, 
             'xpoint_angle': xpoint_angle_dic}


pol_range = self.data['b2fgeo']['nx'] + 2
rad_range = self.data['b2fgeo']['ny'] + 2


solps_dsa_dic = {}
for aa in self.data['dircomp']['Attempt'].keys():
    solps_dsa_dic[aa] = lcm.read_dsa(self.data['dirdata']['simudir'][aa] + '/dsa')

den_list = list(self.data['dircomp']['Attempt'].keys())
aa = den_list[0]

mag_axis_z = self.data['gfile']['g']['zmaxis']
print(mag_axis_z)

Attempt = self.data['dircomp']['Attempt'][aa]
DRT = self.data['dirdata']['outputdir'][aa]['Output']


RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
            usecols = (3)).reshape((rad_range, pol_range))
VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
              usecols = (3)).reshape((rad_range, pol_range))

mag_axis_z = self.data['gfile']['g']['zmaxis']
print(mag_axis_z)
mag_axis_r = self.data['gfile']['g']['rmaxis']
print(mag_axis_r)

angle_list = []

for pol_index in pol_list:
    rloc = RadLoc[:, int(pol_index)][0]
    zloc = VertLoc[:, int(pol_index)][0]
    
    x = rloc - mag_axis_r
    y = zloc - mag_axis_z
    
    
    angle = np.arctan2(y, x)*180 / np.pi
    if angle <= -90:
        angle = angle + 360
    angle_list.append(angle)
    
ln = len(pol_list)
pol_loc = np.zeros(ln)
i = 0

for ii in pol_list:
    pol_loc[i] = int(ii)
    i = i + 1


# plt.figure(figsize=(7,7))
# plt.plot(pol_loc, angle_list, color = 'r', label= 'angle')
# plt.xlabel('poloidal index')
# plt.ylabel('angle')
# plt.title('angle verses poloidal index')

self.data['angle'] = angle_list


 # pol_range = int(self.data['b2fgeo']['nx'] + 2)
 rad_range = int(self.data['b2fgeo']['ny'] + 2)
 
 
 Attempt = self.data['dircomp']['Attempt']
 DRT = self.data['dirdata']['outputdir']['Output']
 XDIM = self.data['b2fgeo']['nx'] + 2
 YDIM = self.data['b2fgeo']['ny'] + 2


        RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
                    usecols = (3)).reshape((YDIM, XDIM))
        VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
                      usecols = (3)).reshape((YDIM,XDIM))


dsa = lcm.read_dsa(self.data['dirdata']['simudir'] + '/dsa')


"Debug tool"
# comp_dsa = np.zeros((int(rad_range), 2))
# comp_dsa[:, 0] = dsa
# comp_dsa[:, 1] = RR_sep

# dsa_dic = {'arclength': arclength, 'interpfunc': interpfunc_dic,
#            'RadLoc': RadLoc, 'VertLoc': VertLoc,
#            'dsa_{}_val'.format(pol_index): RR_sep, 'comp_dsa': comp_dsa}


# pol_range = int(self.data['b2fgeo']['nx'] + 2)
rad_range = int(self.data['b2fgeo']['ny'] + 2)


Attempt = self.data['dircomp']['Attempt']
DRT = self.data['dirdata']['outputdir']['Output']
XDIM = self.data['b2fgeo']['nx'] + 2
YDIM = self.data['b2fgeo']['ny'] + 2
dsa = lcm.read_dsa(self.data['dirdata']['simudir'] + '/dsa')

RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
            usecols = (3)).reshape((YDIM, XDIM))
VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
              usecols = (3)).reshape((YDIM,XDIM))

pol_index = int(pol_loc)
crloc = RadLoc[:, pol_index]
czloc = VertLoc[:, pol_index]
# print(crloc)

arclength, interpfunc_dic = self.RR_sep_calculator(cr = crloc, cz = czloc)

SEP = self.data['DefaultSettings']['SEP']

sep_dist = np.mean([arclength[int(SEP)+ 1], arclength[int(SEP)]])
# print(sep_dist)


RR_sep = arclength - sep_dist

dsa_dic = {'arclength': arclength, 'interpfunc': interpfunc_dic,
            'dsa_{}_val'.format(pol_index): RR_sep}

self.data['dsa']['dsa_{}'.format(pol_index)] = dsa_dic


"Debug tool"
# comp_dsa = np.zeros((int(rad_range), 2))
# comp_dsa[:, 0] = dsa
# comp_dsa[:, 1] = RR_sep

# dsa_dic = {'arclength': arclength, 'interpfunc': interpfunc_dic,
#            'RadLoc': RadLoc, 'VertLoc': VertLoc,
#            'dsa_{}_val'.format(pol_index): RR_sep, 'comp_dsa': comp_dsa}


"Test for left only result"
# LLdsa = crLowerLeft - LLsep
# LRdsa = crLowerRight - LRsep
# ULdsa = crUpperLeft - ULsep
# URdsa = crUpperRight - URsep
# avag_rad = np.zeros(rad_range)
# for j in range(rad_range):
#     avag_rad[j] = np.mean([LLdsa[j], ULdsa[j], URdsa[j], LRdsa[j]])

# crzsum, points = fm.RR_sep_calculator(cr = crLowerLeft, cz = czLowerLeft)

# Rad0Cor = np.loadtxt('{}/Rad0Cor{}'.format(DRT2, str(Attempt)),
#             usecols = (3)).reshape((YDIM, XDIM))
# Vert0Cor = np.loadtxt('{}/Vert0Cor{}'.format(DRT2, str(Attempt)), 
#               usecols = (3)).reshape((YDIM,XDIM))

# comp = np.zeros((int(rad_range), 2))
# comp[:, 0] = Rad0Cor[:, pol_index]
# comp[:, 1] = crLowerLeft

pol_range = int(self.data['b2fgeo'][aa]['nx'] + 2)
# print('xdim is {}'.format(str(pol_range)))
rad_range = int(self.data['b2fgeo'][aa]['ny'] + 2)
# print('ydim is {}'.format(str(rad_range)))

Attempt = self.data['dircomp']['Attempt'][aa]
DRT = self.data['dirdata']['infolderdir'][aa]['outputdir']['Output']

RadLoc = np.loadtxt('{}/RadLoc{}'.format(DRT, str(Attempt)),
            usecols = (3)).reshape((rad_range, pol_range))
VertLoc = np.loadtxt('{}/VertLoc{}'.format(DRT, str(Attempt)), 
              usecols = (3)).reshape((rad_range, pol_range))

pol_index = int(pol_loc)
crloc = RadLoc[:, pol_index]
czloc = VertLoc[:, pol_index]
# print(crloc)

arclength, interpfunc_dic = self.RR_sep_calculator(cr = crloc, cz = czloc)

SEP = self.data['DefaultSettings']['SEP'][aa]

sep_dist = np.mean([arclength[int(SEP)+ 1], arclength[int(SEP)]])
# print(sep_dist)


RR_sep = arclength - sep_dist


"""


