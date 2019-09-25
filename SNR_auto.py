############################
# Date: 24/09/2019
# Title: SHIPS Auto SNR
# Description: Automated script to run SHIPS' SNR map calculations for all stars within directory
# pyKLIP version: 1.1 NOT IMPLEMENTED YET
# Python version: 3 ONLY
############################

import os
import pathlib
rootdir = '/Users/alan/Documents/PhD/Data/SPHERE/IFS'
for subdir, dirs, files in os.walk(rootdir):
    if subdir.startswith("/Users/alan/Documents/PhD/Data/SPHERE/IFS/Q") or subdir.startswith("/Users/alan/Documents/PhD/Data/SPHERE/IFS/HD93403") or subdir.startswith("/Users/alan/Documents/PhD/Data/SPHERE/IFS/P") or subdir.startswith("/Users/alan/Documents/PhD/Data/SPHERE/IFS/D"):
        continue
    #if subdir.startswith("/Users/alan/Documents/PhD/Data/SPHERE/IFS/H") or subdir.startswith("/Users/alan/Documents/PhD/Data/SPHERE/IFS/V"):
    else:
        path = pathlib.PurePath(subdir)
        if path.name != 'IFS':
            #print(path.name)
            os.system('python ships_ifs_run_latest_VIP_SNR_auto.py '+path.name)
            # for file in files:
            #  if file.endswith(".fits"):
            #      path = pathlib.PurePath(subdir)
            #      print(file)
