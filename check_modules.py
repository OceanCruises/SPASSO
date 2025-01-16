#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check Python modules for SPASSO run

Created on Thu Jan 16 10:06:06 2025

@author: lrousselet
"""
import importlib

#module name from pyproject.toml file
module_name = ['termcolor','pandas','netCDF4', 'motuclient', 'matplotlib',
               'cmocean','scipy','xarray','basemap','basemap-data-hires',
               'tabulate','simplekml','requests','importlib-metadate',
               'pylatex','copernicusmarine']
list_installed = []
list_notinstalled = []

for mod in module_name:
    spec = importlib.util.find_spec(mod)
    if spec is not None:
        version = importlib.metadata.version(mod)
        list_installed.append(mod+' '+version)
    else:
        list_notinstalled.append(mod)

print('List of installed packages:')
print("\n".join(map(str,list_installed)))

print('\nList of missing packages:')
print("\n".join(map(str,list_notinstalled)))
    
    
    # try:
    #     map(__import__,mod)
    #     print(mod+" is installed.")
    # except ImportError:
    #     print(mod+" is not installed")
