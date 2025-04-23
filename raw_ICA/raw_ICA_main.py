#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 10:39:12 2021

@author: leupinv
"""
import raw_ICA.raw_ICA_helper as hp
import raw_ICA.raw_ICA_constants as cs
import base.files_in_out as files_in_out
import base.base_constants as b_cs



for g_n in b_cs.G_N[4:5]:

    for cond in cs.condition[0]:
        files=files_in_out.GetFiles(filepath=cs.datafolder,eeg_format=cs.eeg_format,g_num=g_n)
        files.select_condition(cond)
        #fif_taskfiles=files.condition_files
        files.get_info(end_fix=17)
        epochs=hp.Epoch_HP(files)
        #raw_ICA.get_eog_epochs(thresh=None)
        #raw_ICA.get_exp_epochs()
        #raw_ICA.get_ecg_epochs()
        #raw_ICA.run_infoICA()
        #raw_ICA.run_fastICA()
        
        #raw_ICA.save_epochs()

        epochs.select_ICA_components()
        #
        # except:
        #     if raw_ICA.isica():
        #         pass
        #     else:
        #

        
        #files_in_out.save_report(files,raw_ICA.report)
        #files_in_out.save_report(files,raw_ICA.report,final=True)

        

        
        
# for g_n in b_cs.G_bad_card:
#     for cond in cs.condition[0]:
#         files=files_in_out.GetFiles(filepath=cs.datafolder,eeg_format=cs.eeg_format,g_num=g_n)
#         files.select_condition(cond)
#         #fif_taskfiles=files.condition_files
#         files.get_info(end_fix=17)
#         raw_ICA=hp.Epoch_HP(files)
#         raw_ICA.raw.set_annotations(None)
#         raw_ICA.run_fastICA()

#%%

#%%
