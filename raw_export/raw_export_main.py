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


for g_num in b_cs.G_N[2]:
    raw = hp.load_raw_data(g_num)
    onsets, durations, descriptions, event_table = hp.parse_annotations(raw)
    event_table = event_table.sort_values(by="time").reset_index(drop=True)
    # Fix HEP cardiac phase to match the corresponding VEP phase (if different)
    corrected_card_phases = event_table["card_phase"].copy()
    # Loop over the event table to find hep-vep pairs
    for i in range(1, len(event_table)):
        if event_table.loc[i, "event"] == "vep" and event_table.loc[i - 1, "event"] == "hep":
            # Check if they are close in time (same trial)
            time_diff = event_table.loc[i, "time"] - event_table.loc[i - 1, "time"]
            print(time_diff)
            if time_diff < 2:  # Adjust threshold as needed
                hep_idx = i - 1
                vep_idx = i
                # Update the HEP card_phase to match the VEP one
                if event_table.loc[hep_idx, "card_phase"] != event_table.loc[vep_idx, "card_phase"]:
                    corrected_card_phases[hep_idx] = event_table.loc[vep_idx, "card_phase"]
    # Apply corrected phases
    event_table["card_phase"] = corrected_card_phases
    new_raw = hp.set_new_annotations(raw, onsets, durations, descriptions)
    hp.save_outputs(new_raw, event_table, g_num, files_in_out)
#%%
