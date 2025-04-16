#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 16:14:48 2021

@author: leupinv
"""
import os
import platform
from files_in_out import getListOfFiles,GetFiles
import base_constants as bs

platform.system()

# define starting datafolder 

if platform.system()=='Darwin':
    os.chdir('/Users/leupinv/switchdrive/BBC/WP1/data/EEG/tsk/')

else:
    os.chdir('d:/switchdrive/BBC/WP1/data/EEG/tsk')

target_dir='/Users/leupinv/switchdrive/BBC/WP1/data/EEG/tsk/ana/MNE/source'


# list_dir=os.listdir(target_dir)

# for g in list_dir:
#     g_list=getListOfFiles(g,target_dir)

for g in bs.G_N:
   # g_dir=target_dir+f'{g}/{g}_evoked'
    files=getListOfFiles(target_dir,g_num=None)
    for file in files:
        if g in file:

            filename=file.split('/')[-1]
            
            print(filename)
            dir_name=file[:-(len(filename))]
            filename=filename[4:]
            new_name=f'{g}_n_tsk_{filename}'
            new_dir=dir_name+new_name
            print(new_dir)
            os.rename(file,new_dir)
            
       