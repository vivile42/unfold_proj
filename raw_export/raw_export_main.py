import pandas as pd
import mne
import sys
import platform
import os
import scipy.io
import base.files_in_out as files_in_out
import numpy as np

if platform.system()=='Darwin':
    os.chdir('/Volumes/BBC/BBC/WP1/data/EEG/tsk/')
    #sys.path.append('/Users/leupinv/BBC/WP1/data/Code/python/BBC')
    #if this doesn't work pound line right above this, restart the kernel, rerun the cell.
    #Then uncomment the same line and rerun

else:
    #os.chdir('Z:/BBC/WP1/data/EEG/tsk')
    #sys.path.append('C:/Users/Vivi/switchdrive/BBC/WP1/data/Code/python/BBC')
    os.chdir('Z:/BBC/WP1/data/EEG/tsk')
from base.files_in_out import getListOfFiles,GetFiles
import base.base_constants as b_cs
import raw_export.raw_export_helper as hp


for g_num in b_cs.G_N:
    raw = hp.load_raw_data(g_num)
    onsets, durations, descriptions, event_table = hp.parse_annotations(raw)
    event_table = event_table.sort_values(by="time").reset_index(drop=True)
    new_raw = hp.set_new_annotations(raw, onsets, durations, descriptions)
    hp.save_outputs(new_raw, event_table, g_num, files_in_out)
#%%

#%%
