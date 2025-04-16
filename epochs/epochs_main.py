#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 10:39:12 2021

@author: leupinv
"""
import epochs.epochs_helper as hp
import epochs.epochs_constants as cs
import base.files_in_out as files_in_out
import base.base_constants as b_cs



for g_n in b_cs.G_O_ica[1:]:
    for cond in cs.condition[1]:
        files=files_in_out.GetFiles(filepath=cs.datafolder,eeg_format=cs.eeg_format,g_num=g_n)
        files.select_condition(cond)
        #fif_taskfiles=files.condition_files
        files.get_info(end_fix=17)
        epochs=hp.Epoch_HP(files)
        #epochs.get_eog_epochs(thresh=None)
        #epochs.get_exp_epochs()
        #epochs.get_ecg_epochs()
        epochs.run_infoICA()
        #epochs.run_fastICA()
        
        #epochs.save_epochs()

        # try:
        #     epochs.select_ICA_components()
        #
        # except:
        #     if epochs.isica():
        #         pass
        #     else:
        #

        
        #files_in_out.save_report(files,epochs.report)
        #files_in_out.save_report(files,epochs.report,final=True)

        

        
        
# for g_n in b_cs.G_bad_card:
#     for cond in cs.condition[0]:
#         files=files_in_out.GetFiles(filepath=cs.datafolder,eeg_format=cs.eeg_format,g_num=g_n)
#         files.select_condition(cond)
#         #fif_taskfiles=files.condition_files
#         files.get_info(end_fix=17)
#         epochs=hp.Epoch_HP(files)
#         epochs.raw.set_annotations(None)
#         epochs.run_fastICA()

#%%
