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



for g_n in b_cs.G_N:

    for cond in cs.condition[0]:
        files=files_in_out.GetFiles(filepath=cs.datafolder,eeg_format=cs.eeg_format,g_num=g_n)
        files.select_condition(cond)

        files.get_info(end_fix=17)
        epochs=hp.Epoch_HP(files)
        epochs.select_ICA_components()

