############################
# Date: 24/09/2019
# Title: SHIPS Auto SNR
# Description: Automated script to run SHIPS' SNR map calculations for all stars within directory
# pyKLIP version: 1.1 NOT IMPLEMENTED YET
# Python version: 3 ONLY
############################

import os
import subprocess
import pathlib
rootdir = '/home/alan/data/Backup_macbook/SPHERE/IRDIS'
for subdir, dirs, files in os.walk(rootdir):
    if subdir.startswith("/home/alan/data/Backup_macbook/SPHERE/IRDIS/H"):
        continue
    #if subdir.startswith("/Users/alan/Documents/PhD/Data/SPHERE/IFS/H") or subdir.startswith("/Users/alan/Documents/PhD/Data/SPHERE/IFS/V"):
    #else:
    elif subdir.startswith("/home/alan/data/Backup_macbook/SPHERE/IRDIS/Q"):
        path = pathlib.PurePath(subdir)
        if path.name != 'IRDIS':
            #print(path.name)
            os.system('python ships_irdis_SNR_auto.py ' +path.name)
            # for file in files:
            #  if file.endswith(".fits"):
            #      #os.system('python ships_ifs_run_latest_VIP_SNR_auto.py ', path.name , file)
            #      subprocess.check_call(["python","ships_ifs_run_latest_VIP_SNR_auto.py",path.name,file])
