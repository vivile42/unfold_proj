# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 12:18:21 2021

@author: Engi
"""

import epochs.epochs_helper as hp
import epochs.epochs_constants as cs
import base.files_in_out as files_in_out
import base.base_constants as b_cs

g_n= 'g12'
cond='n'
files=files_in_out.GetFiles(filepath=cs.datafolder,eeg_format=cs.eeg_format,g_num=g_n)
files.select_condition(cond)
files.get_info(end_fix=17)
        
epochs=hp.Epoch_HP(files)
#epochs.run_ICA()


# In[test]

epochs.get_eog_epochs(ch_name='C30',thresh=70e-6)

epochs.raw.plot(events=epochs.eog_epochs.events,n_channels=64)

print(epochs.eog_epochs.events)


# In[save]
epochs.get_exp_epochs()
epochs.get_ecg_epochs()
epochs.save_epochs()