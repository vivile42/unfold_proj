import mne
import numpy as np
import pandas as pd

def load_raw_data(g_num, datafolder='raw'):
    filepath = f"{datafolder}/{g_num}/{g_num}_n_tsk_ICA_rec-raw.fif"
    raw = mne.io.read_raw_fif(filepath, preload=True)
    return raw

def parse_annotations(raw):
    sfreq = raw.info["sfreq"]
    annotations = raw.annotations
    sample_indices = np.round(annotations.onset * sfreq).astype(int)  # Unused but might be useful later

    event_data = {
        "time": [],
        "event": [],
        "awareness": [],
        "duration": [],
        "card_phase": [],
        "rsp_phase": [],
    }
    filtered_onsets, filtered_durations, filtered_descriptions = [], [], []

    for i, desc in enumerate(annotations.description):
        tags = desc.split("/")
        if "bad" in tags[0].lower():
            continue
        if len(tags)<=1:
            continue

        onset = annotations.onset[i]
        duration = annotations.duration[i]

        if "normal" in tags[1].lower() and "correct" in tags[2].lower():
            if "hep" in tags[0].lower() and tags[4] in ["R", "R2"]:
                awareness = "aware" if tags[3] == "aware" else "unaware"
                card_phase = "dia" if tags[-3] == "R2" else "sys"
                rsp_phase = "inh" if tags[-2] == "inh" else "exh"
                condition_name = "hep"
            elif "vep" in tags[0].lower() and tags[3].lower() in ["aware", "unaware"]:
                awareness = "unaware" if "unaware" in tags[3].lower() else "aware"
                card_phase = "dia" if tags[-2] == "dia" else "sys"
                rsp_phase = "inh" if tags[-1] == "inh" else "exh"
                condition_name = "vep"
            else:
                continue

            filtered_onsets.append(onset)
            filtered_durations.append(duration)
            filtered_descriptions.append(condition_name)

            event_data["time"].append(onset)
            event_data["event"].append(condition_name)
            event_data["awareness"].append(awareness)
            event_data["duration"].append(duration)
            event_data["card_phase"].append(card_phase)
            event_data["rsp_phase"].append(rsp_phase)
    event_data=pd.DataFrame(event_data)
    # Fix HEP cardiac phase to match the corresponding VEP phase (if different)
    corrected_card_phases = event_data["card_phase"].copy()
    # Loop over the event table to find hep-vep pairs
    for i in range(1, len(event_data)):
        if event_data.loc[i, "event"] == "vep" and event_data.loc[i - 1, "event"] == "hep":
            # Check if they are close in time (same trial)
            time_diff = event_data.loc[i, "time"] - event_data.loc[i - 1, "time"]
            if time_diff < 2:  # Adjust threshold as needed
                hep_idx = i - 1
                vep_idx = i
                # Update the HEP card_phase to match the VEP one
                if event_data.loc[hep_idx, "card_phase"] != event_data.loc[vep_idx, "card_phase"]:
                    corrected_card_phases[hep_idx] = event_data.loc[vep_idx, "card_phase"]
    # Apply corrected phases
    event_data["card_phase"] = corrected_card_phases

    return filtered_onsets, filtered_durations, filtered_descriptions, event_data

def set_new_annotations(raw, onsets, durations, descriptions):
    new_annotations = mne.Annotations(onsets, durations, descriptions)
    new_raw = raw.copy()
    new_raw.set_annotations(new_annotations)
    new_raw.pick_types(eeg=True)
    return new_raw

def save_outputs(new_raw, event_table, g_num, files_in_out):
    eeg_format = 'fif'
    eeg_exp = 'tsk'
    datafolder = 'raw'

    files = files_in_out.GetFiles(filepath=datafolder, eeg_format=eeg_format, g_num=g_num)
    files.select_condition('n')
    files.get_info(end_fix=29)
    files.current_filename = files.current_filename[:9]

    output_file_eeg = files.out_filename(type_sig='deconv', file_end='eeg_data_sys.set')
    output_file_table = files.out_filename(type_sig='deconv', file_end='event_table_sys.csv')

    mne.export.export_raw(output_file_eeg, raw=new_raw, overwrite=True)
    event_table.to_csv(output_file_table, index=False)

