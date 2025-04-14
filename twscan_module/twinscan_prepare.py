# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 20:23:48 2025

@author: ychuang
"""




class twscan_assist:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data


    def pair_dic(self, keys, values):
        
        # Use zip() to pair the keys with the values
        zipped_pairs = zip(keys, values)
        
        # Convert the zipped pairs into a dictionary
        result_dic = dict(zipped_pairs)
        
        return result_dic
    
    
    def twinscan_prep(self, ta, keylist_b, scan_style):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == True:
            
            if self.DF.series_flag == 'twin_scan':
                
                color_list = ['red', 'orange', 'green', 'blue', 'purple']
                
                color_dic = self.pair_dic(keys = keylist_b, values = color_list)
                
                # print('check color dic:')
                # print(color_dic)
                
                scan_list = []
                iter_key = []
                
                
                for tb in keylist_b:
                    
                    if scan_style == 'tempscan':
                        
                        it_in = (ta, tb)
                    
                    elif scan_style == 'denscan':
                        
                        it_in = (tb, ta)
                    
                    else:
                        print('twinscan_plot_method, please check the scan_style!')
                    
                    
                    nx = self.data['b2fgeo']['nx']
                    ny = self.data['b2fgeo']['ny']
                    dat_struc = {'nx': nx, 'ny': ny}

                    
                    mid_ne_pro = self.data['midplane_profile'][it_in[0]][it_in[1]]['mid_ne']
                    mid_te_pro = self.data['midplane_profile'][it_in[0]][it_in[1]]['mid_te']
                        
                    
                    if scan_style == 'tempscan':
                        
                        scan_add = '{:.1f} eV'.format(mid_te_pro[0])
                    
                    elif scan_style == 'denscan':
                        
                        scan_add = '{:.2E} '.format(mid_ne_pro[0])
                    
                    else:
                        print('twinscan_plot_method, please check the scan_style!')
                    
                    scan_list.append(scan_add)
                    iter_key.append(it_in)
                
                
                # print('NT scan list: {}'.format(ta))
                # print(scan_list)
                
                
                if scan_style == 'tempscan':
                    
                    midne = self.data['midplane_profile'][ta]['4.115']['mid_ne']

                    scan_title = '{:.2E}'.format(midne[0])
                
                elif scan_style == 'denscan':
                    
                    midte = self.data['midplane_profile']['5.512'][ta]['mid_te']

                    scan_title = '{:.1f}'.format(midte[0])
                
                else:
                    print('twinscan_plot_method, please check the scan_style!')
                
                label_dic = self.pair_dic(keys = keylist_b, values = scan_list)
                

                return iter_key, color_dic, scan_title, label_dic