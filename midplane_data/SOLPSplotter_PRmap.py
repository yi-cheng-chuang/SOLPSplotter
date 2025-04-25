# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 17:27:01 2023

@author: user
"""


from load_coordinate.load_coord_method import load_coordgeo_method
from fit_data.fitting_method import fit_method_collection
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np


class RP_mapping:
    
    
    def __init__(self, DF, data, lcm: load_coordgeo_method, fmc: fit_method_collection):
        
        self.DF = DF
        self.data = data
        self.fmc = fmc
        self.lcm = lcm
    
  
#-----------------------------------------------------------------------------

#calculate separatrix index, R-Rsep, R-Rsep psiN mapping  
  
    
    def calc_sep_index(self, psi, sep_value):
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
                if psi[y] <= sep_value:
                    index_low.append(y)
                elif psi[y] > sep_value:
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
    
    
    
    def calc_RRsep_method(self, geo, radgrid, vertgrid, gfile_data, psiNinterp_function, plotRR):
               
        pol_range = int(geo['nx'] + 2)
        rad_range = int(geo['ny'] + 2)
        
        RadLoc = radgrid
        VertLoc = vertgrid
        psiNinterp_RBS = psiNinterp_function
        
        
        # mag_axis_z = gfile_data['zmaxis']
        # mag_axis_z = 0.15
        # print(mag_axis_z)
        
        
        "Calculate midplane R for Z at the magnetic axis"
        
        "Calculate weight"
        
        DEV = self.DF.DEV
        
        if DEV == 'mast':
            mag_axis_z = gfile_data['zmaxis']
            
            crup = RadLoc[:, 58]
            crlow = RadLoc[:, 60]
            czup = VertLoc[:, 58]
            czlow = VertLoc[:, 60]
            
            average_pair = (58, 60)
        
        elif DEV == 'mastu':
            mag_axis_z = 0.15
            
            crup = RadLoc[:, 71]
            crlow = RadLoc[:, 73]
            czup = VertLoc[:, 71]
            czlow = VertLoc[:, 73]
            
            average_pair = (71, 73)
        
        else:
            print("please check the device input")
        

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
        
        sep_index_dic = self.calc_sep_index(psi = psi_solps_mid, sep_value = 1)
        sep_index_high = int(sep_index_dic['index'][1])
        
        
        weight_psi = (psi_solps_mid[sep_index_high] -1)/(psi_solps_mid[sep_index_high] - psi_solps_mid[sep_index_high -1])
        
        R_sep = (1 - weight_psi)*mid_R[sep_index_high] + weight_psi*mid_R[sep_index_high -1]
        
        print('rsep is {:.2f}'.format(R_sep))
        
        # psi_high = psi_solps_mid[sep_index_high]
        # psi_low = psi_solps_mid[sep_index_high -1]
        
        # print('psi_high is {:.3f}'.format(psi_high))
        # print('psi_low is {:.3f}'.format(psi_low))
        
        R_Rsep = mid_R - R_sep
        
        dsa_psi_func = interp1d(R_Rsep, psi_solps_mid, kind='quadratic', fill_value = 'extrapolate')
        
        
        midR_psi_func = interp1d(mid_R, psi_solps_mid, kind='quadratic', fill_value = 'extrapolate')
        
        
        rsep_fit = midR_psi_func(1.37)
        
        print('rsep_mid return psiN value: {:.2f}'.format(rsep_fit))
        
        
        midplane_dic = {'weight': weight_mid, 'mid_choice': mid_choice, 
                        'mid_R': mid_R, 'mid_Z': mid_Z, 
                        'psi_solps_mid': psi_solps_mid, 
                        'weight_psi': weight_psi, 'R_Rsep': R_Rsep,
                        'dsa_psi_func': dsa_psi_func, 
                        'midR_psi_func': midR_psi_func,
                        'average_pair': average_pair}
        
        psi_dsa_dic = self.fmc.dsa_psi_fit(dsa= R_Rsep, psi= psi_solps_mid)
        
        psi_dsa_ratio = psi_dsa_dic['dsa_psi_fitcoe'][0]
        
        
        
        if DEV == 'mast':
            pol_list = [57, 58, 59, 60, 61]
            
            
        
        elif DEV == 'mastu':
            pol_list = [70, 71, 72, 73, 74]
            
                   
        if plotRR:
            plt.figure()
            for in_pol in pol_list:
                crloc = RadLoc[:, int(in_pol)]
                czloc = VertLoc[:, int(in_pol)]
                plt.plot(crloc, czloc, color = 'g', label = 'R&Zlocation')
            plt.plot(mid_R, mid_Z, color = 'r', label= 'R_Rsep')
            plt.xlabel('R: [m]')
            plt.ylabel('Z: [m]')
    
    
        return midplane_dic, psi_dsa_ratio, sep_index_high
    
    
    
    def calc_RRsep(self, plotRR, plot_psi_dsa_align):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        
        if withshift == False and withseries == False:
            b2fgeo = self.data['b2fgeo']
            radloc = self.data['grid']['RadLoc']
            vertloc = self.data['grid']['VertLoc']
            gfile = self.data['gfile']['g']
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            
            
            midplane_dic, psi_dsa_ratio, sep_index_high = self.calc_RRsep_method(geo = b2fgeo, 
                        radgrid = radloc, vertgrid = vertloc, gfile_data = gfile, psiNinterp_function = psiNinterp_RBS, 
                        plotRR = plotRR)
            
            self.data['DefaultSettings']['psi_dsa'] = psi_dsa_ratio
            self.data['midplane_calc'] = midplane_dic
            
            self.data['DefaultSettings']['sep_index_RRsep'] = sep_index_high
        
        
        elif withshift == True and withseries == False:        
            
            sep_index_dic = {}
            midplane_series_dic = {}
            psi_dsa_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                
                b2fgeo = self.data['b2fgeo'][aa]
                radloc = self.data['grid']['RadLoc'][aa]
                vertloc = self.data['grid']['VertLoc'][aa]
                gfile = self.data['gfile']['g']
                psiNinterp_RBS = self.data['gfile']['gcomp'][aa]['interp_dic']['RBS']
                
                               
                midplane_dic, psi_dsa_ratio, sep_index_high = self.calc_RRsep_method(geo = b2fgeo, 
                            radgrid = radloc, vertgrid = vertloc, gfile_data = gfile, psiNinterp_function = psiNinterp_RBS, 
                            plotRR = plotRR)
                
                midplane_series_dic[aa] = midplane_dic
                psi_dsa_dic[aa] = psi_dsa_ratio
                sep_index_dic[aa] = sep_index_high
            
            self.data['DefaultSettings']['psi_dsa'] = psi_dsa_dic
            self.data['midplane_calc'] = midplane_series_dic
            self.data['DefaultSettings']['sep_index_RRsep'] = sep_index_dic
                
                
            if plot_psi_dsa_align:
                color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                             'dot7': 'blue', 'one': 'purple'}
                plt.figure()
                for in_sh in self.data['dircomp']['multi_shift']:
                    # plt.plot(midplane_dic[in_sh]['psi_solps_mid'], midplane_dic[in_sh]['R_Rsep'], 'o-', 
                    #   label= 'R-Rsep:[m] with modify {} m'.format(self.data['dircomp']['shift_dic'][in_sh]))
                    plt.plot(midplane_dic[in_sh]['psi_solps_mid'], midplane_dic[in_sh]['R_Rsep'], 'o-', 
                      color= color_dic[in_sh])
                plt.xlabel('psiN')
                plt.title('psiN R-Rsep plot')
                # plt.legend()
                    
           
        elif withshift == False and withseries == True:
            den_rep = list(self.data['dircomp']['Attempt'].keys())[0]
            
            b2fgeo = self.data['b2fgeo']
            radloc = self.data['grid']['RadLoc']
            vertloc = self.data['grid']['VertLoc']
            gfile = self.data['gfile']['g']
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            
            
            midplane_dic, psi_dsa_ratio, sep_index_high = self.calc_RRsep_method(geo = b2fgeo, 
                        radgrid = radloc, vertgrid = vertloc, gfile_data = gfile, psiNinterp_function = psiNinterp_RBS, 
                        plotRR = plotRR)
            
            self.data['DefaultSettings']['psi_dsa'] = psi_dsa_ratio
            self.data['midplane_calc'] = midplane_dic
            self.data['DefaultSettings']['sep_index_RRsep'] = sep_index_high
        
        elif self.withshift == True and self.withseries == True:
            print('calc_RRsep function is not there yet!')
        
        else:
            print('calc_RRsep function has a bug')
            
            

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

    

    """
    find the separatrix location using dsa and confirm that the dsa = arclength
    
    """

    
    def calc_sep_dsa_method(self, RadLoc, VertLoc, dsa, jxa):
        
                
        pol_index = jxa
        crloc = RadLoc[:, pol_index]
        czloc = VertLoc[:, pol_index]
        # print(crloc)
        
        arclength, interpfunc_dic = self.RR_sep_calculator(cr = crloc, cz = czloc)
        
        sep_dist = arclength - dsa
        # print(sep_dist)
        sep_val = np.mean(sep_dist)
        
        sep_index_dic = self.calc_sep_index(psi = arclength, sep_value = sep_val)
        sep_index_high = int(sep_index_dic['index'][1])
        
        
        return sep_dist, sep_index_high
        
        
    def calc_sep_dsa(self):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == False:
            
            rad_grid = self.data['grid']['RadLoc']
            vert_grid = self.data['grid']['VertLoc']     
            dsa = self.lcm.read_dsa(self.data['dirdata']['simudir'] + '/dsa')
            jxa = self.data['b2mn']['jxa']
            
            dist, index = self.calc_sep_dsa_method(RadLoc = rad_grid, VertLoc = vert_grid, 
                                            dsa = dsa, jxa = jxa)
            
            self.data['dist'] = dist
            self.data['DefaultSettings']['sep_index_dsa'] = index
        
        elif withshift == True and withseries == False:
            
            dist_dic = {}
            index_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                rad_grid = self.data['grid']['RadLoc'][aa]
                vert_grid = self.data['grid']['VertLoc'][aa]
                
                
                series_compare = self.DF.series_compare
                
                if series_compare == True:
                    
                    dsa = self.lcm.read_dsa(self.data['dirdata']['simudir'][aa]['fixed'] + '/dsa')
                
                else:
                    
                    dsa = self.lcm.read_dsa(self.data['dirdata']['simudir'][aa] + '/dsa')
                    
                

                jxa = self.data['b2mn'][aa]['jxa']
                
                dist, index = self.calc_sep_dsa_method(RadLoc = rad_grid, VertLoc = vert_grid, 
                                                dsa = dsa, jxa = jxa)
                
                dist_dic[aa] = dist
                index_dic[aa] = index
            
            
            self.data['dist'] = dist_dic
            self.data['DefaultSettings']['sep_index_dsa'] = index_dic
        
        elif withshift == False and withseries == True:
            
            den_list = list(self.data['dircomp']['Attempt'].keys())
            aa = den_list[0]
            
            rad_grid = self.data['grid']['RadLoc']
            vert_grid = self.data['grid']['VertLoc']
            
            
            series_flag = self.DF.series_flag
            
            
            if series_flag == 'twin_scan':
                
                nf = aa[0]
                tf = aa[1]
                
                dsa = self.lcm.read_dsa(self.data['dirdata']['simudir'][nf][tf] + '/dsa')
            
            else:
                
                dsa = self.lcm.read_dsa(self.data['dirdata']['simudir'][aa] + '/dsa')


           
            jxa = self.data['b2mn']['jxa']
            
            dist, index = self.calc_sep_dsa_method(RadLoc = rad_grid, VertLoc = vert_grid, 
                                            dsa = dsa, jxa = jxa)
            
            self.data['dist'] = dist
            self.data['DefaultSettings']['sep_index_dsa'] = index
            
        
        elif withshift == True and withseries == True:
            print('calc_sep_dsa function is not there yet!')
        
        else:
            print('calc_sep_dsa function has a bug')
                
            
    
    """
    Calculate dsa using arclength method

    """       
    
        
    def calc_dsa_method(self, RadLoc, VertLoc, SEP, pol_loc):
        
        
        pol_index = int(pol_loc)
        crloc = RadLoc[:, pol_index]
        czloc = VertLoc[:, pol_index]
        # print(crloc)
        
        arclength, interpfunc_dic = self.RR_sep_calculator(cr = crloc, cz = czloc)
        
        sep_dist = np.mean([arclength[int(SEP)], arclength[int(SEP)- 1]])
        # print(sep_dist)
                
        RR_sep = arclength - sep_dist
        
        return arclength, interpfunc_dic, RR_sep
        
       
    def calc_dsa(self, pol_loc):
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == False:
            
            rad_grid = self.data['grid']['RadLoc']
            vert_grid = self.data['grid']['VertLoc']
            sep_loc = self.data['DefaultSettings']['sep_index_dsa']
            
            
            arclength, interpfunc_dic, RR_sep = self.calc_dsa_method(RadLoc = rad_grid, 
                                 VertLoc = vert_grid, SEP = sep_loc,  pol_loc = pol_loc)
            
            dsa_dic = {'arclength': arclength, 'interpfunc': interpfunc_dic,
                        'dsa_{}_val'.format(pol_loc): RR_sep}
            
            self.data['dsa']['dsa_{}'.format(pol_loc)] = dsa_dic
            
        
        elif withshift == True and withseries == False:
            
            dsa_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                
                rad_grid = self.data['grid']['RadLoc'][aa]
                vert_grid = self.data['grid']['VertLoc'][aa]
                sep_loc = self.data['DefaultSettings']['sep_index_dsa'][aa]
                
                
                arclength, interpfunc_dic, RR_sep = self.calc_dsa_method(RadLoc = rad_grid, 
                                     VertLoc = vert_grid, SEP = sep_loc, pol_loc = pol_loc)
                
                dsa_dic[aa] = {'arclength': arclength, 'interpfunc': interpfunc_dic,
                            'dsa_{}_val'.format(pol_loc): RR_sep}
            
            self.data['dsa']['dsa_{}'.format(pol_loc)] = dsa_dic
    
    
    
        elif withshift == False and withseries == True:
            
            den_list = list(self.data['dircomp']['Attempt'].keys())
            aa = den_list[0]
            
            rad_grid = self.data['grid']['RadLoc']
            vert_grid = self.data['grid']['VertLoc']
            sep_loc = self.data['DefaultSettings']['sep_index_dsa']
            
            arclength, interpfunc_dic, RR_sep = self.calc_dsa_method(RadLoc = rad_grid, 
                                 VertLoc = vert_grid, SEP = sep_loc, pol_loc = pol_loc)
            
            dsa_dic = {'arclength': arclength, 'interpfunc': interpfunc_dic,
                        'dsa_{}_val'.format(pol_loc): RR_sep}
            
            self.data['dsa']['dsa_{}'.format(pol_loc)] = dsa_dic
        
        
        elif withshift == True and withseries == True:
            print('calc_dsa_method function is not there yet!')
        
        else:
            print('calc_dsa_method has a bug')


#---------------------------------------------------------------------------          

# flux expansion calculation

                
    def calc_flux_expansion_method(self, arcR, RR_sep, ped_index):
        
        
        arcR_inv = list(reversed(arcR))
        RRsep_inv = list(reversed(RR_sep))
                   
        arcR_cut = []
        RRsep_cut = []
        
        for p_in in ped_index:
            arcR_cut.append(arcR_inv[p_in])
            RRsep_cut.append(RRsep_inv[p_in])
            
                   
        flux_fit_dic = self.fmc.flux_expand_fit(RRsep = RRsep_cut, arclength = arcR_cut)
        
        flux_expand = flux_fit_dic['flux_fitcoe'][0]
        
        return flux_expand
    
    
    def calc_flux_expansion(self, pol_loc, ped_index, iter_index):
       
       withshift = self.DF.withshift
       withseries = self.DF.withseries
       
       if withshift == False and withseries == False:
           
           arcR = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]
           RR_sep = self.data['midplane_calc']['R_Rsep']
       
       elif withshift == True and withseries == False:
           
           arcR = self.data['dsa']['dsa_{}'.format(pol_loc)][iter_index]['dsa_{}_val'.format(pol_loc)]
           RR_sep = self.data['midplane_calc'][iter_index]['R_Rsep']
       
       elif withshift == False and withseries == True:
           
           arcR = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]
           RR_sep = self.data['midplane_calc']['R_Rsep']
       
       elif withshift == True and withseries == True:
           print('calc_flux_expansion is not there yet, to be continue...')
       
       flux_expand = self.calc_flux_expansion_method(arcR = arcR, RR_sep = RR_sep, 
                                                     ped_index = ped_index)
       
       return flux_expand
   
    
    def calc_flux_expansion_wait(self, pol_loc, ped_index):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        
        if withshift == False and withseries == False:
            
            arcR = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]
            RR_sep = self.data['midplane_calc']['R_Rsep']
            
            flux_expand = self.calc_flux_expansion_method(arcR = arcR, RR_sep = RR_sep, 
                                                          ped_index = ped_index)
            
            return flux_expand
        
        elif withshift == True and withseries == False:
            
            flux_expand_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
            
                arcR = self.data['dsa']['dsa_{}'.format(pol_loc)][aa]['dsa_{}_val'.format(pol_loc)]
                RR_sep = self.data['midplane_calc'][aa]['R_Rsep']
                
                flux_expand = self.calc_flux_expansion_method(arcR = arcR, RR_sep = RR_sep, 
                                                              ped_index = ped_index)
                
                flux_expand_dic[aa] = flux_expand
            
            return flux_expand_dic
                
                        
        elif withshift == False and withseries == True:
            
            arcR = self.data['dsa']['dsa_{}'.format(pol_loc)]['dsa_{}_val'.format(pol_loc)]
            RR_sep = self.data['midplane_calc']['R_Rsep']
            
            
            flux_expand = self.calc_flux_expansion_method(arcR = arcR, RR_sep = RR_sep, 
                                                          ped_index = ped_index)
            
            return flux_expand
            
            
        elif withshift == True and withseries == True:
            print('calc_flux_expansion is not there yet, to be continue...')
        
        else:
            print('calc_flux_expansion has a bug.')


#---------------------------------------------------------------------------------------           


      
# Additional dsa check function 1: left_test          
        
    def left_dsa_test_method(self, geo, pol_index, sep_loc, Attempt, DRT2):
        
        rad_range = int(geo['ny'] + 2)
        
        "Test for left only result"
        
        crLowerLeft = geo['crx'][pol_index,:,0]
        crLowerRight = geo['crx'][pol_index,:,1]
        crUpperLeft = geo['crx'][pol_index,:,2]
        crUpperRight = geo['crx'][pol_index,:,3]
        czLowerLeft = geo['cry'][pol_index,:,0]
        czLowerRight = geo['cry'][pol_index,:,1]
        czUpperLeft = geo['cry'][pol_index,:,2]
        czUpperRight = geo['cry'][pol_index,:,3]
        
        SEP = sep_loc
        
        LLsep = np.mean([crLowerLeft[int(SEP)], crLowerLeft[int(SEP)-1]])
        LRsep = np.mean([crLowerRight[int(SEP)], crLowerRight[int(SEP)-1]])
        ULsep = np.mean([crUpperLeft[int(SEP)], crUpperLeft[int(SEP)-1]])
        URsep = np.mean([crUpperRight[int(SEP)], crUpperRight[int(SEP)-1]])
        
        
        LLdsa = crLowerLeft - LLsep
        LRdsa = crLowerRight - LRsep
        ULdsa = crUpperLeft - ULsep
        URdsa = crUpperRight - URsep
        avag_rad = np.zeros(rad_range)
        for j in range(rad_range):
            avag_rad[j] = np.mean([LLdsa[j], ULdsa[j], URdsa[j], LRdsa[j]])
      
        crzsum, points = self.fmc.RR_sep_calculator(cr = crLowerLeft, cz = czLowerLeft)
        
                      
        XDIM = geo['nx'] + 2
        YDIM = geo['ny'] + 2   
        
        Rad0Cor = np.loadtxt('{}/Rad0Cor{}'.format(DRT2, str(Attempt)),
                    usecols = (3)).reshape((YDIM, XDIM))
        # Vert0Cor = np.loadtxt('{}/Vert0Cor{}'.format(DRT2, str(Attempt)), 
        #               usecols = (3)).reshape((YDIM,XDIM))
        
        compare_left = np.zeros((int(rad_range), 2))
        compare_left[:, 0] = Rad0Cor[:, pol_index]
        compare_left[:, 1] = crLowerLeft
        
        return compare_left
    
        
    def left_test(self, pol_loc):
        pol_pos = int(pol_loc)
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        if withshift == False and withseries == False:
            
            b2fgeo = self.data['b2fgeo']
            pol_pos = int(pol_loc)
            sep = self.data['DefaultSettings']['sep_index_dsa']
            attempt = self.data['dircomp']['Attempt']
            drt2 = self.data['dirdata']['outputdir']['Output2']
            
            
            compare_left = self.left_dsa_test_method(geo = b2fgeo, pol_index = pol_pos, 
                                                sep_loc = sep, Attempt = attempt, 
                                                DRT2 = drt2)
            
            
            self.data['dsa_compare_left'] = compare_left
        
        elif withshift == True and withseries == False:
            
            compare_left_dic = {}
            
            for aa in self.data['dircomp']['multi_shift']:
                
                
                SEP = self.data['DefaultSettings']['SEP'][aa]
                b2fgeo = self.data['b2fgeo'][aa]
                pol_pos = int(pol_loc)
                sep = self.data['DefaultSettings']['sep_index_dsa'][aa]
                attempt = self.data['dircomp']['Attempt'][aa]
                drt2 = self.data['dirdata']['infolderdir'][aa]['outputdir']['Output2']
                
                
                compare_left = self.left_dsa_test_method(geo = b2fgeo, pol_index = pol_pos, 
                                                    sep_loc = sep, Attempt = attempt, 
                                                    DRT2 = drt2)
                
                compare_left_dic[aa] = compare_left
                
                self.data['dsa_compare_left'] = compare_left_dic
        
        elif withshift == False and withseries == True:
            
            den_list = list(self.data['dircomp']['Attempt'].keys())
            aa = den_list[0]
            
            SEP = self.data['DefaultSettings']['SEP']
            b2fgeo = self.data['b2fgeo']
            pol_pos = int(pol_loc)
            sep = self.data['DefaultSettings']['sep_index_dsa'][aa]
            attempt = self.data['dircomp']['Attempt'][aa]
            drt2 = self.data['dirdata']['outputdir'][aa]['Output2']
            
            compare_left = self.left_dsa_test_method(geo = b2fgeo, pol_index = pol_pos, 
                                                sep_loc = sep, Attempt = attempt, 
                                                DRT2 = drt2)
            
            self.data['dsa_compare_left'] = compare_left
            
        
        elif withshift == True and withseries == True:
            print('left_test function is not there yet!')
        
        else:
            print('left_test function has a bug')


#-----------------------------------------------------------------------------






"""
backup


"""


