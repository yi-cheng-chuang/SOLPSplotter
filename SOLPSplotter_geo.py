# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 11:54:55 2023

@author: user
"""

import numpy as np
from SOLPS_load_directory import load_directory 
import load_coord_method as lcm
from scipy import interpolate


class load_geometry(load_directory):
    
    def __init__(self, DefaultSettings):
        load_directory.__init__(self, DefaultSettings)
                


   
    def check_b2mn(self, itername):
        if self.withshift == False and self.withseries == False:
            if self.data['b2mn']['jxa'] == None:
                b2mn = lcm.scrape_b2mn(self.data['dirdata']['simudir']
                                      + '/b2mn.dat')
                self.data['b2mn'] = b2mn
            else:
                pass
        elif self.withshift == True and self.withseries == False:
            if self.data['b2mn'][itername]['jxa'] == None:
                b2mn = lcm.scrape_b2mn(self.data['dirdata']['infolderdir'][itername]['simudir']
                                          + '/b2mn.dat')
            else:
                pass
        elif self.withshift == False and self.withseries == True:
            if self.data['b2mn']['jxa'] == None:
                b2mn = lcm.scrape_b2mn(self.data['dirdata']['simudir'][itername]
                                      + '/b2mn.dat')
            else:
                pass
                
#loadgeometry                  
    """
    loadb2fgmtry method for single case, changing aspect ratio cases and 
    changing density cases
    
    """

    def loadgeo_method(self, attempt_loc, simufile_loc, g_data, shift_value):
        
        try:
            b2mn = lcm.scrape_b2mn(attempt_loc + '/b2mn.dat')
            
        except:
            print('can not generate b2mn')
        
            
        try:
            geo = lcm.read_b2fgmtry(simufile_loc + '/baserun/b2fgmtry')
            # print(type(geo))
        except:
            print('can not generate geo')
        
        psiN = (g_data['psirz'] - g_data['simag']) / (g_data['sibry'] - g_data['simag'])

        dR = g_data['rdim'] / (g_data['nw'] - 1)
        dZ = g_data['zdim'] / (g_data['nh'] - 1)
        
        gZ = np.zeros(g_data['nh'])
        for i in range(g_data['nh']):
            gZ[i] = g_data['zmid'] - 0.5 * g_data['zdim'] + i * dZ
                

        gR = np.zeros(g_data['nw'])
        for i in range(g_data['nw']):
            gR[i] = g_data['rleft'] + i * dR + float(shift_value)
            

        psiNinterp_RBS = interpolate.RectBivariateSpline(gR, gZ, np.transpose(psiN))
        psiNinterp_2d = interpolate.interp2d(gR, gZ, psiN, kind = 'cubic')
        psiNinterp_RGI = interpolate.RegularGridInterpolator((gR, gZ), np.transpose(psiN))
        
        interp_dic = {}
        interp_dic['RBS'] = psiNinterp_RBS
        interp_dic['2d'] = psiNinterp_2d
        interp_dic['RGI'] = psiNinterp_RGI
        
        
        gfilesum = {'psiN': psiN, 'dR': dR, 'dZ': dZ, 'gR': gR, 'gZ': gZ,
                    'check': 'yeah! shift is {} and series is {}'.format(self.withshift, self.withseries),
                    'interp_dic': interp_dic}
            
        return b2mn, geo, gfilesum
            
    """
    load_solpsgeo function is to load the geometric file for solps and is
    for single case, changing aspect ratio cases and changing density cases
       
    """     
    def load_solpsgeo(self):
        
        if self.terminal == True:
            
            print('check the g file dir')
            print(self.data['dirdata']['gdir'])
            g_loc = self.data['dirdata']['gdir'][0]
            
            
        elif self.terminal == False:
            g_loc = self.data['dirdata']['gbase'] + '/MAST__RMP_results/g027205.00275_efitpp'
            

        gfile_data = lcm.loadg(g_loc)
        self.data['gfile']['g'] = gfile_data
        
        
        if self.withshift == False and self.withseries == False:
            simudir = self.data['dirdata']['simudir']
            simutop = self.data['dirdata']['simutop']
            shift = self.data['dircomp']['shift_value']
            
            b2mn, geo, gfilesum = self.loadgeo_method(attempt_loc = simudir, 
                            simufile_loc = simutop, g_data = gfile_data, shift_value = shift)
            
            self.data['b2mn'] = b2mn
            self.data['b2fgeo'] = geo
            self.data['gfile']['gcomp'] = gfilesum
            
        elif self.withshift == True and self.withseries == False:
            b2mn_dic = {}
            geo_dic = {}
            gfilesum_dic = {}
            for shiftname in self.data['dircomp']['multi_shift']:
                
                simudir = self.data['dirdata']['simudir'][shiftname]
                simutop = self.data['dirdata']['simutop'][shiftname]
                shift = self.data['dircomp']['shift_dic'][shiftname]
                b2mn, geo, gfilesum = self.loadgeo_method(attempt_loc = simudir, 
                        simufile_loc = simutop, g_data = gfile_data, shift_value = shift)
                b2mn_dic[shiftname] = b2mn
                geo_dic[shiftname] = geo
                gfilesum_dic[shiftname] = gfilesum
            
            self.data['b2mn'] = b2mn_dic
            self.data['b2fgeo'] = geo_dic      
            self.data['gfile']['gcomp'] = gfilesum_dic
                  
        elif self.withshift == False and self.withseries == True:
            
            seriesname = list(self.data['dircomp']['Attempt'].keys())[0]
            
            if self.series_flag == 'twin_scan':
                
                nf = seriesname[0]
                tf = seriesname[1]
                        
                simudir = self.data['dirdata']['simudir'][nf][tf]
            
            else:
                
                simudir = self.data['dirdata']['simudir'][seriesname]
                
            simutop = self.data['dirdata']['simutop']
            shift = self.data['dircomp']['shift_value']                          
            b2mn, geo, gfilesum = self.loadgeo_method(attempt_loc = simudir, 
                    simufile_loc = simutop, g_data = gfile_data, shift_value = shift)
            
            
            self.data['b2mn'] = b2mn
            self.data['b2fgeo'] = geo
            self.data['gfile']['gcomp'] = gfilesum
                            
        elif self.withshift == True and self.withseries == True:
            print('load_solpsgeo is not there yet, to be continue...')
        
        else:
            print('There is a bug')


#---------------------------------------------------------------------------

#psi function mapping

    """
    calcpsi method for single case, changing aspect ratio cases and 
    changing density cases
       
    """ 
    
    def calcpsi_method_avcr(self, geo, psi_RBS):
        
        pol_range = int(geo['nx'] + 2)
        rad_range = int(geo['ny'] + 2)
        psival = np.zeros((rad_range, pol_range))
        
        # Assuming you have four 98x38 matrices stored in variables

        # Replace these with your actual matrices or generate them as needed
        
        crLowerLeft = geo['crx'][:,:,0]
        crLowerRight = geo['crx'][:,:,1]
        crUpperLeft = geo['crx'][:,:,2]
        crUpperRight = geo['crx'][:,:,3]
        czLowerLeft = geo['cry'][:,:,0]
        czLowerRight = geo['cry'][:,:,1]
        czUpperLeft = geo['cry'][:,:,2]
        czUpperRight = geo['cry'][:,:,3]
        
        
        # Create a list of matrices
        crx_collect = [crLowerLeft, crLowerRight, crUpperLeft, crUpperRight]
        cry_collect = [czLowerLeft, czLowerRight, czUpperLeft, czUpperRight]
        
        # Convert the list of matrices to a NumPy array
        crx_array = np.array(crx_collect)
        cry_array = np.array(cry_collect)
        
        # Calculate the average along the first axis (axis=0) to get the desired matrix
        rad_mean = np.mean(crx_array, axis = 0)
        vert_mean = np.mean(cry_array, axis = 0)
        
        
        rad_loc = rad_mean.transpose()
        vert_loc = vert_mean.transpose()
        
        # Now, average_matrix contains the average values of the corresponding elements in the four matrices
        # self.data['grid_check'] = rad_loc
        for pol_loc in range(pol_range):
            for i in range(rad_range):
                # print(i)
                psival[i, pol_loc] = psi_RBS(rad_loc[i, pol_loc], 
                                                      vert_loc[i, pol_loc])
        
                
        return rad_loc, vert_loc, psival, pol_range, rad_range
    
    
    def calcpsi_avcr(self):
            
        if self.withshift == False and self.withseries == False:
            
            b2fgeo = self.data['b2fgeo']
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            
            
            RadLoc, VertLoc, psival, pol_range, rad_range = self.calcpsi_method_avcr(geo = b2fgeo, 
                                                                    psi_RBS = psiNinterp_RBS)
            
            
            coord_dic = {'RadLoc': RadLoc, 'VertLoc': VertLoc}
            self.data['grid'] = coord_dic
            self.data['psi']['psival'] = psival
            self.data['DefaultSettings']['XDIM'] = pol_range
            self.data['DefaultSettings']['YDIM'] = rad_range
        
        elif self.withshift == True and self.withseries == False:
            psival_dic = {}
            RadLoc_dic = {}
            VertLoc_dic = {}
            xdim_dic = {}
            ydim_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                
                b2fgeo = self.data['b2fgeo'][aa]
                psiNinterp_RBS = self.data['gfile']['gcomp'][aa]['interp_dic']['RBS']
                
                RadLoc, VertLoc, psival, pol_range, rad_range = self.calcpsi_method_avcr(geo = b2fgeo, 
                                                                    psi_RBS = psiNinterp_RBS)
                
                psival_dic[aa] = psival
                RadLoc_dic[aa] = RadLoc
                VertLoc_dic[aa] = VertLoc
                xdim_dic[aa] = pol_range
                ydim_dic[aa] = rad_range
            
            coord_dic = {'RadLoc': RadLoc_dic, 'VertLoc': VertLoc_dic}
            self.data['grid'] = coord_dic
            self.data['psi']['psival'] = psival_dic
            self.data['DefaultSettings']['XDIM'] = xdim_dic
            self.data['DefaultSettings']['YDIM'] = ydim_dic
                
            
        elif self.withshift == False and self.withseries == True:

            b2fgeo = self.data['b2fgeo']
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            
            
            RadLoc, VertLoc, psival, pol_range, rad_range = self.calcpsi_method_avcr(geo = b2fgeo, 
                                                                psi_RBS = psiNinterp_RBS)
            
            coord_dic = {'RadLoc': RadLoc, 'VertLoc': VertLoc}
            self.data['grid'] = coord_dic
            self.data['psi']['psival'] = psival       
            self.data['DefaultSettings']['XDIM'] = pol_range
            self.data['DefaultSettings']['YDIM'] = rad_range
        
        elif self.withshift == True and self.withseries == True:
            print('calcpsi_avcr function is not there yet!')
        
        else:
            print('calcpsi_avcr function has a bug')


    

#-----------------------------------------------------------------------------
    

    def calcpsi_method(self, rad_range, pol_range, psi_RBS, attempt_number, outputdir):
               
        psival = np.zeros((rad_range, pol_range))
       
        RadLoc = np.loadtxt('{}/RadLoc{}'.format(outputdir, str(attempt_number)),
                    usecols = (3)).reshape((rad_range, pol_range))
        VertLoc = np.loadtxt('{}/VertLoc{}'.format(outputdir, str(attempt_number)), 
                      usecols = (3)).reshape((rad_range, pol_range))
        
        # coord_dic = {'RadLoc': RadLoc, 'VertLoc': VertLoc}
            
        for pol_loc in range(pol_range):
            for i in range(rad_range):
                # print(i)
                psival[i, pol_loc] = psi_RBS(RadLoc[i, pol_loc], 
                                                      VertLoc[i, pol_loc])
        return RadLoc, VertLoc, psival, pol_range, rad_range


                
    def calcpsi(self):
            
        if self.withshift == False and self.withseries == False:
            n_pol = int(self.data['b2fgeo']['nx'] + 2)
            # print('xdim is {}'.format(pol_range))
            n_rad = int(self.data['b2fgeo']['ny'] + 2)
            # print('ydim is {}'.format(rad_range))
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            # psiNinterp_RGI = self.data['gfile']['gcomp']['interp_dic']['RGI'] 
            # psiNinterp_2d = self.data['gfile']['gcomp']['interp_dic']['2d']
            Attempt = self.data['dircomp']['Attempt']
            DRT = self.data['dirdata']['outputdir']['Output']
            

            RadLoc, VertLoc, psival, pol_range, rad_range = self.calcpsi_method(rad_range = n_rad, 
                                pol_range = n_pol, psi_RBS = psiNinterp_RBS, 
                                attempt_number = Attempt, outputdir = DRT)
            
            
            coord_dic = {'RadLoc': RadLoc, 'VertLoc': VertLoc}
            self.data['grid'] = coord_dic
            self.data['psi']['psival'] = psival
            self.data['DefaultSettings']['XDIM'] = pol_range
            self.data['DefaultSettings']['YDIM'] = rad_range
        
        elif self.withshift == True and self.withseries == False:
            psival_dic = {}
            RadLoc_dic = {}
            VertLoc_dic = {}
            xdim_dic = {}
            ydim_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                
                n_pol = int(self.data['b2fgeo'][aa]['nx'] + 2)
                # print('xdim is {}'.format(str(pol_range)))
                n_rad = int(self.data['b2fgeo'][aa]['ny'] + 2)
                # print('ydim is {}'.format(str(rad_range)))
                psiNinterp_RBS = self.data['gfile']['gcomp'][aa]['interp_dic']['RBS']
                Attempt = self.data['dircomp']['Attempt'][aa]
                DRT = self.data['dirdata']['outputdir'][aa]['Output']
                
                RadLoc, VertLoc, psival, pol_range, rad_range = self.calcpsi_method(rad_range = n_rad, 
                                    pol_range = n_pol, psi_RBS = psiNinterp_RBS, 
                                    attempt_number = Attempt, outputdir = DRT)
                
                psival_dic[aa] = psival
                RadLoc_dic[aa] = RadLoc
                VertLoc_dic[aa] = VertLoc
                xdim_dic[aa] = pol_range
                ydim_dic[aa] = rad_range
            
            coord_dic = {'RadLoc': RadLoc_dic, 'VertLoc': VertLoc_dic}
            self.data['grid'] = coord_dic
            self.data['psi']['psival'] = psival_dic
            self.data['DefaultSettings']['XDIM'] = xdim_dic
            self.data['DefaultSettings']['YDIM'] = ydim_dic
                
            
        elif self.withshift == False and self.withseries == True:
            # print('we are working on calcpsi for series case')
            
            series_rep = list(self.data['dircomp']['Attempt'].keys())[0]
            n_pol = int(self.data['b2fgeo']['nx'] + 2)
            # print('xdim is {}'.format(pol_range))
            n_rad = int(self.data['b2fgeo']['ny'] + 2)
            # print('ydim is {}'.format(rad_range))
            psiNinterp_RBS = self.data['gfile']['gcomp']['interp_dic']['RBS']
            Attempt = self.data['dircomp']['Attempt'][series_rep]
            DRT = self.data['dirdata']['outputdir'][series_rep]['Output']
            
            
            RadLoc, VertLoc, psival, pol_range, rad_range = self.calcpsi_method(rad_range = n_rad, 
                                pol_range = n_pol, psi_RBS = psiNinterp_RBS, 
                                attempt_number = Attempt, outputdir = DRT)
            
            coord_dic = {'RadLoc': RadLoc_dic, 'VertLoc': VertLoc_dic}
            self.data['grid'] = coord_dic
            self.data['psi']['psival'] = psival       
            self.data['DefaultSettings']['XDIM'] = pol_range
            self.data['DefaultSettings']['YDIM'] = rad_range
        
        elif self.withshift == True and self.withseries == True:
            print('calcpsi function is not there yet!')
        
        else:
            print('calcpsi function has a bug!')
        
    

#----------------------------------------------------------------------------


    def calcpsi_block_method(self, file_loc, shift):
        b2mn = lcm.scrape_b2mn(file_loc + '/b2mn.dat')
        print(file_loc)
        print(b2mn)
        jxa = b2mn['jxa']
        pol_index = int(jxa)
        
        simu_dir = file_loc.rsplit("/",1)[0]
        fname = file_loc.rsplit("/",1)[1]
        print('plasmf input: {}'.format(simu_dir))
        print('filename input: {}'.format(fname))
        
        geo = lcm.read_b2fgmtry(simu_dir + '/baserun/b2fgmtry')
        rad_range = int(geo['ny'] + 2)
        
        crLowerLeft = geo['crx'][pol_index,:,0]
        crLowerRight = geo['crx'][pol_index,:,1]
        crUpperLeft = geo['crx'][pol_index,:,2]
        crUpperRight = geo['crx'][pol_index,:,3]
        czLowerLeft = geo['cry'][pol_index,:,0]
        czLowerRight = geo['cry'][pol_index,:,1]
        czUpperLeft = geo['cry'][pol_index,:,2]
        czUpperRight = geo['cry'][pol_index,:,3]
        
        
        g = lcm.loadg(self.data['dirdata']['gbase'] 
                               + '/MAST__RMP_results/g027205.00275_efitpp')
        # g = lcm.loadg(self.data['dirdata']['gdir'][0])
        psiN = (g['psirz'] - g['simag']) / (g['sibry'] - g['simag'])

        dR = g['rdim'] / (g['nw'] - 1)
        dZ = g['zdim'] / (g['nh'] - 1)
        
        gZ = np.zeros(g['nh'])
        for i in range(g['nh']):
            gZ[i] = g['zmid'] - 0.5 * g['zdim'] + i * dZ
            
        
        gR = np.zeros(g['nw'])
        for i in range(g['nw']):
            gR[i] = g['rleft'] + i * dR + float(shift)
            

        psiNinterp_RBS = interpolate.RectBivariateSpline(gR, gZ, np.transpose(psiN))
        
        
        psi_solps_RBS = np.zeros(rad_range)
        for i in range(rad_range):
            psi_LL_RBS = psiNinterp_RBS(crLowerLeft[i], czLowerLeft[i])
            psi_UL_RBS = psiNinterp_RBS(crUpperLeft[i], czUpperLeft[i])
            psi_LR_RBS = psiNinterp_RBS(crLowerRight[i], czLowerRight[i])
            psi_UR_RBS = psiNinterp_RBS(crUpperRight[i], czUpperRight[i])
            psi_solps_RBS[i] = np.mean([psi_LL_RBS, psi_UL_RBS, 
                                        psi_LR_RBS, psi_UR_RBS])
        
        return psi_solps_RBS



#---------------------------------------------------------------------------
   
#1D psi function mapping  

    def calcpsi_1D_method(self, geo, psiNinterp_dic, pol_loc, no_coord_avg_check, grid_data):
        
        
        rad_range = int(geo['ny'] + 2)
        psiNinterp_RGI = psiNinterp_dic['RGI'] 
        psiNinterp_2d = psiNinterp_dic['2d']
        psiNinterp_RBS = psiNinterp_dic['RBS']
        
        
        pol_index = int(pol_loc)
        
        crLowerLeft = geo['crx'][pol_index,:,0]
        crLowerRight = geo['crx'][pol_index,:,1]
        crUpperLeft = geo['crx'][pol_index,:,2]
        crUpperRight = geo['crx'][pol_index,:,3]
        czLowerLeft = geo['cry'][pol_index,:,0]
        czLowerRight = geo['cry'][pol_index,:,1]
        czUpperLeft = geo['cry'][pol_index,:,2]
        czUpperRight = geo['cry'][pol_index,:,3]
            
        
        "dsa check"
        # kk = np.zeros((int(rad_range),3))
        # kk[:, 0] = dsa
        # kk[:, 1] = avag_rad
        
        # del_dsa = np.zeros(rad_range)
        # for ia in range(rad_range):
        #     del_dsa[ia] = avag_rad[ia] - dsa[ia]
        
        # kk[:, 2] = del_dsa
        # self.data['psi']['dsa_{}_val'.format(pol_loc)] = avag_rad
        
        
        crz_LL = np.stack([crLowerLeft.ravel(), czLowerLeft.ravel()], -1)  # shape (N, 2) in 2d
        crz_UL = np.stack([crUpperLeft.ravel(), czUpperLeft.ravel()], -1)
        crz_LR = np.stack([crLowerRight.ravel(), czLowerRight.ravel()], -1)  # shape (N, 2) in 2d
        crz_UR = np.stack([crUpperRight.ravel(), czUpperRight.ravel()], -1)
        psi_solps_LL = psiNinterp_RGI(crz_LL)
        psi_solps_UL = psiNinterp_RGI(crz_UL)
        psi_solps_LR = psiNinterp_RGI(crz_LR)
        psi_solps_UR = psiNinterp_RGI(crz_UR)
        
        
        psi_solps = np.zeros(rad_range)
        for i in range(rad_range):
            psi_solps[i] = np.mean([psi_solps_LL[i], psi_solps_UL[i], 
                                    psi_solps_LR[i], psi_solps_UR[i]])
        
        psi_solps_2d = np.zeros(rad_range)
        for i in range(rad_range):
            psi_LL_2d = psiNinterp_2d(crLowerLeft[i], czLowerLeft[i])
            psi_UL_2d = psiNinterp_2d(crUpperLeft[i], czUpperLeft[i])
            psi_LR_2d = psiNinterp_2d(crLowerRight[i], czLowerRight[i])
            psi_UR_2d = psiNinterp_2d(crUpperRight[i], czUpperRight[i])
            psi_solps_2d[i] = np.mean([psi_LL_2d, psi_UL_2d, 
                                       psi_LR_2d, psi_UR_2d])
         
        
        psi_solps_RBS = np.zeros(rad_range)
        for i in range(rad_range):
            psi_LL_RBS = psiNinterp_RBS(crLowerLeft[i], czLowerLeft[i])
            psi_UL_RBS = psiNinterp_RBS(crUpperLeft[i], czUpperLeft[i])
            psi_LR_RBS = psiNinterp_RBS(crLowerRight[i], czLowerRight[i])
            psi_UR_RBS = psiNinterp_RBS(crUpperRight[i], czUpperRight[i])
            psi_solps_RBS[i] = np.mean([psi_LL_RBS, psi_UL_RBS, 
                                        psi_LR_RBS, psi_UR_RBS])
        
        "Only work for original case"
        # GF = eq.equilibrium(gfile= self.data['dirdata']['gdir'][0])
        # print(type(GF.psiN._func))
        # psi_solps_GF = np.zeros(rad_range)
        # for i in range(rad_range):
        #     psi_LL_GF = GF.psiN(crLowerLeft[i], czLowerLeft[i])
        #     psi_UL_GF = GF.psiN(crUpperLeft[i], czUpperLeft[i])
        #     psi_LR_GF = GF.psiN(crLowerRight[i], czLowerRight[i])
        #     psi_UR_GF = GF.psiN(crUpperRight[i], czUpperRight[i])
        #     psi_solps_GF[i] = np.mean([psi_LL_GF, psi_UL_GF, 
        #                                 psi_LR_GF, psi_UR_GF])
        
        
        if no_coord_avg_check:
                   
            if self.data['psi']['psival'] == None:
                self.calcpsi()
                crloc = grid_data[:, pol_index]
                czloc = grid_data[:, pol_index]
            else:
                crloc = grid_data[:, pol_index]
                czloc = grid_data[:, pol_index]
            
            psi_solps_cp = np.zeros(rad_range)
            # psi_solps_cp = psiNinterp_RBS(crloc, czloc)
            for i in range(rad_range):
                psi_CP = psiNinterp_RBS(crloc[i], czloc[i])
                psi_solps_cp[i] = psi_CP
        
        
            psival = np.zeros((int(rad_range), 4))
            psival[:, 0] = psi_solps
            psival[:, 1] = psi_solps_2d
            psival[:, 2] = psi_solps_RBS
            psival[:, 3] = psi_solps_cp
        else:
            psival = np.zeros((int(rad_range), 3))
            psival[:, 0] = psi_solps
            psival[:, 1] = psi_solps_2d
            psival[:, 2] = psi_solps_RBS
        
        
        # index_dic = self.calc_sep_index(psi = psi_solps_RBS, rad_range = rad_range)        
        # self.data['DefaultSettings']['SEP'] = index_dic['index'][0]
        
        return psival
        


    def calcpsi_1D(self, pol_loc, no_coord_avg_check):
        
        if self.withshift == False and self.withseries == False:
            
            if self.withshift == False and self.withseries == False:
                
                b2fgeo = self.data['b2fgeo']
                Interp_dic = self.data['gfile']['gcomp']['interp_dic']
                psiN = self.data['psi']['psival']
                
            
            psival = self.calcpsi_1D_method(geo = b2fgeo, psiNinterp_dic = Interp_dic, 
                        pol_loc = pol_loc, no_coord_avg_check = no_coord_avg_check, 
                        grid_data = psiN)
            
            
            self.data['psi']['psi_{}_val'.format(pol_loc)] = psival
        
        elif self.withshift == True and self.withseries == False:
            
            psival_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                
                b2fgeo = self.data['b2fgeo'][aa]
                Interp_dic = self.data['gfile']['gcomp'][aa]['interp_dic'] 
                psiN = self.data['psi']['psival'][aa]
                
                
                psival = self.calcpsi_1D_method(geo = b2fgeo, psiNinterp_dic = Interp_dic, 
                            pol_loc = pol_loc, no_coord_avg_check = no_coord_avg_check, 
                            grid_data = psiN)
                
                psival_dic[aa] = psival
            
            self.data['psi']['psi_{}_val'.format(pol_loc)] = psival_dic
                   
        elif self.withshift == False and self.withseries == True:
            
            
            b2fgeo = self.data['b2fgeo']
            Interp_dic = self.data['gfile']['gcomp']['interp_dic'] 
            psiN = self.data['psi']['psival']
            
            
            psival = self.calcpsi_1D_method(geo = b2fgeo, psiNinterp_dic = Interp_dic, 
                        pol_loc = pol_loc, no_coord_avg_check = no_coord_avg_check, 
                        grid_data = psiN)
                
            self.data['psi']['psi_{}_val'.format(pol_loc)] = psival
        
        elif self.withshift == True and self.withseries == True:
            print('calcpsi_1D function is not there yet!')
        
        else:
            print('calcpsi_1D function has a bug')
            


#---------------------------------------------------------------------------

# load vessel

    def load_vessel_method(self, fdir):
        # try:
        #     WallFile = np.loadtxt('{}/mesh.extra'.format(self.data['dirdata']['tbase']))
        # except:
        #     print('mesh.extra file not found! Using vvfile.ogr instead')
        #     WallFile=None
        
        try:
            VVFILE = np.loadtxt('{}/baserun/vvfile.ogr'.format(fdir))
        except:
            print('load_vessel_method has a bug!')

        return VVFILE
    
    
    def load_vessel(self):
        if self.withshift == False and self.withseries == False:
            filedir = self.data['dirdata']['simutop']
            vessel_file = self.load_vessel_method(fdir = filedir)
            self.data['vessel'] = vessel_file
        
        elif self.withshift == True and self.withseries == False:
            vessel_file_dic = {}
            for aa in self.data['dircomp']['multi_shift']:
                filedir = self.data['dirdata']['simutop'][aa]
                vessel_file = self.load_vessel_method(fdir = filedir)
                vessel_file_dic[aa] = vessel_file
            
            self.data['vessel'] = vessel_file_dic
        
        elif self.withshift == False and self.withseries == True:
            # series_rep = list(self.data['dircomp']['Attempt'].keys())[0]
            filedir = self.data['dirdata']['simutop']
            vessel_file = self.load_vessel_method(fdir = filedir)
            self.data['vessel'] = vessel_file
        
        elif self.withshift == True and self.withseries == True:
            print('load_vessel function is not there yet!')
        
        else:
            print('load_vessel function has a bug')


        


'Way to generate align transport coefficient'
    
"""
m = len(yd)

one_trans = b2tp.InputfileParser(one_list[0], plot= False)
ond = one_trans['1'].T
onki = one_trans['3'].T
onke = one_trans['4'].T
onx= ond[:,0]  #the coordinate here is R-R_sep
fd = ond[:,1]
fki = onki[:,1]
fke = onke[:,1]

d_func = interpolate.interp1d(x, yd, fill_value = 'extrapolate')
ond[:,1] = d_func(onx)
ki_func = interpolate.interp1d(x, yki, fill_value = 'extrapolate')
onki[:,1] = ki_func(onx)
ke_func = interpolate.interp1d(x, yke, fill_value = 'extrapolate')
onke[:,1] = ke_func(onx)


b = b2tp.Generate(cod, CoeffID=1, SpeciesID=2, M=[1])
c = b2tp.WriteInputfile(file='b2.transport.inputfile_align_{}_{}'.format(shift_a, n), points= one_trans ,M_1 = True, M=[1])
"""
            

"""
backup:

try:
    if self.withshift == False and self.withseries == False:
        geo = lcm.read_b2fgmtry(self.data['dirdata']['simutop'] 
                               + '/baserun/b2fgmtry')
    elif self.withshift == True and self.withseries == False:
        geo = lcm.read_b2fgmtry(self.data['dirdata']['simutop'][itername] 
                                   + '/baserun/b2fgmtry')
    elif self.withshift == False and self.withseries == True:
        geo = lcm.read_b2fgmtry(self.data['dirdata']['simutop'] 
                                   + '/baserun/b2fgmtry')
    # print(type(geo))
except:
    print('can not generate geo')
    
    

# g = lcm.loadg(self.data['dirdata']['gdir'][0])
g = lcm.loadg(self.data['dirdata']['gbase'] 
                       + '/MAST__RMP_results/g027205.00275_efitpp')
        

if self.withshift == False and self.withseries == False:
    b2mn = lcm.scrape_b2mn(self.data['dirdata']['simudir']
                          + '/b2mn.dat')
elif self.withshift == True and self.withseries == False:
    b2mn = lcm.scrape_b2mn(self.data['dirdata']['simudir'][itername]
                              + '/b2mn.dat')
elif self.withshift == False and self.withseries == True:
    b2mn = lcm.scrape_b2mn(self.data['dirdata']['simudir'][itername]
                          + '/b2mn.dat')

            elif self.withshift == True and self.withseries == False:
                geo = lcm.read_b2fgmtry(self.data['dirdata']['simutop'][itername] 
                                           + '/baserun/b2fgmtry')
            elif self.withshift == False and self.withseries == True:
                geo = lcm.read_b2fgmtry(self.data['dirdata']['simutop'] 
                                           + '/baserun/b2fgmtry')
            # print(type(geo))
            
            if self.withshift == False and self.withseries == False:
                shift = self.data['dircomp']['shift_value']
                # print(shift)
            elif self.withshift == True and self.withseries == False:
                shift = self.data['dircomp']['shift_dic'][itername]
                # print(shift)
            elif self.withshift == False and self.withseries == True:
                shift = self.data['dircomp']['shift_value']
                # print(shift)
                
            if self.withshift == False and self.withseries == False:
                b2mn = lcm.scrape_b2mn(self.data['dirdata']['simudir']
                                      + '/b2mn.dat')
            elif self.withshift == True and self.withseries == False:
                b2mn = lcm.scrape_b2mn(self.data['dirdata']['simudir'][itername]
                                          + '/b2mn.dat')
            elif self.withshift == False and self.withseries == True:
                b2mn = lcm.scrape_b2mn(self.data['dirdata']['simudir'][itername]
                                      + '/b2mn.dat')



            ---- eliminated geometry input code ----
        
"""


        
            
            
        
    
