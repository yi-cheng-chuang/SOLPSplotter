# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 11:54:55 2023

@author: user
"""

class B2plotter:
    def __init__(self, Shot=None, Attempts=None, Parameters=None, 
                 DefaultSettings=None, **kwargs):  
        self.Shot = Shot
        
        "Attempts"
        if isinstance(Attempts, list):
            self.Attempts = Attempts
        else:    
            print('Attempt has to be a list')
            
        
        "Parameters"
        if isinstance(Parameters, list):
            self.Parameters = Parameters
        else:
            print('parameter has to be a list')
            
        if Parameters is None:
            print('There is no parameters input')
        else:
            self.Parameters = Parameters
            
        if isinstance(DefaultSettings, dict):
            self.DefaultSettings = DefaultSettings
        else:
            print('parameter has to be a dictionary')
             
            
        "DefaultSettings"
        if DefaultSettings is None:
            print('There is no input defaultsettings')
        else:
            self.DefaultSettings = DefaultSettings
            
        
                 
    "Add and remove elements from parameter"        
    def add_par(self, new_par):
        if new_par not in self.Parameters:
            self.Parameters.append(new_par)
    
    def remove_par(self, new_par):
        if new_par in self.Parameters:
            self.Parameters.remove(new_par)
    
        
    "Add and remove elements from defaultsettings"
    def add_default(self, new_set):
        keylist = []
        for key, value in self.DefaultSettings.items():
            keylist.append(key)
        if new_set.keys() not in keylist:
            self.DefaultSettings.update(new_set)
    
    def remove_default(self, new_set):
        keylist = []
        for key, value in self.DefaultSettings.items():
            keylist.append(key)
        if new_set.keys() in keylist:
            
        
        if new_par in self.Parameters:
            self.Parameters.remove(new_par)
            
    






