# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 20:24:37 2025

@author: ychuang
"""

from SOLPS_input.input_setting import Setting_dic
import re



class grab_aptn_method:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    def s_number(self, text):
        # sd = ss.Setting_dic()
        
        
        
        if self.DF.DEV == 'mast':
            
            if self.DF.withshift == False and self.DF.withseries == False:
                name = text.split("/",-1)[-2]
                nu = int(name.split('_')[0])
                
                
            elif self.DF.withshift == False and self.DF.withseries == True:
                
                series_flag = self.DF.series_flag
                
                if series_flag == 'change_den':
                    name = text.split("\\",-1)[-1]
                    nu = re.findall('\d+\.\d+', name)
                    nu.append(name.split('_')[0])
                    # print(nu)
                elif series_flag == 'eireneN':
                    name = text.split("\\",-1)[-1]
                    nu = re.findall('\d+', name)
                    nu.append(name.split('_')[0])
                    # print(nu)
                elif series_flag == 'change_temp':
                    name = text.split("\\",-1)[-1]
                    nu = re.findall('\d+\.\d+', name)
                    nu.append(name.split('_')[0])
                
                elif series_flag == 'two_compare':
                    name = text.split("\\",-1)[-1]
                    nu = re.findall('\d+', name)
                    nu.append(name.split('_')[0])
                    # print(nu)
                
                elif series_flag == 'twin_scan':
                    name = text.split("/",-1)[-1]
                    nu = re.findall('\d+\.\d+', name)
                    nu.append(name.split('_')[0])
                
                else:
                    print('check the series flag')
                
            elif self.DF.withshift == True and self.DF.withseries == False:
                name = text.split("/",-1)[-2]
                nu = int(name.split('_')[0])
                
            elif self.DF.withshift == True and self.DF.withseries == True:
                print('unexpected situation, please check the parameter setting')
            else:
                print('There is a bug in s_number function')
        
        else:
            print('DEV setting is not MAST!')
            
        
        

        return [nu, name]
            



    def mastu_atp_number(self, text, usage):
        # sd = ss.Setting_dic()
        
        if self.DF.DEV == 'mastu':
            
            if self.DF.withshift == False and self.DF.withseries == False:
                
                if usage == 'load_dir':
                    
                    name = text.split("/",-1)[-1]
                    
                    # print(name)
                    
                    nu = int(name.split('_')[0])
                
                elif usage == 'transcoe':
                    
                    name = text.split("/",-1)[-2]
                    
                    # print(name)
                    
                    nu = int(name.split('_')[0])
                
            else:
                print('mastu_atp_number function is not there yet!')
        
        else:
            print('DEV setting is not MAST!')
            
        
        

        return [nu, name]



    def atp_number(self, text):
        
        # sd = ss.Setting_dic()
        
        if self.DF.withshift == False and self.DF.withseries == False:
            name = text.split("/",-1)[-2]
            nu = int(name.split('_')[0])
        elif self.DF.withshift == False and self.DF.withseries == True:
          
            if self.DF.series_flag == 'twin_scan':
                divider = "\\"
                if divider in text:
                    # print('{} is in directory'.format(divider))
                    name = text.split("\\",-1)[-1]
                else:
                    name = text.split("/",-1)[-1]
                    # print(name)
                
                nu_list = re.findall('\d+\.\d+', name)
                nu_tuple = (nu_list[0], nu_list[1])
                nu = [nu_tuple, name.split('_')[0]]
            
            else:
                print('check the series flag')
            
        elif self.DF.withshift == True and self.DF.withseries == False:
            name = text.split("/",-1)[-2]
            nu = int(name.split('_')[0])
        elif self.DF.withshift == True and self.DF.withseries == True:
            print('unexpected situation, please check the parameter setting')
        else:
            print('There is a bug in s_number function')

        return nu



