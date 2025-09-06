Code accompanying Leupin and Britz, Interoceptive ability determines whether pre-stimulus heartbeat-evoked activity can predict awareness at the visual threshold.
We here provide the code for preprocessing, deconvolution analysis, erp analysis and statistics. EEG data can be made available upon request.

# Code Organization

Each folder is generally organized with a **main**, **helper**, and **constants** script.  
The base folder contains helper functions to filter through the data directories.

- **Main scripts** → contain the code that must be run.  
- **Helper scripts** → contain helper functions and classes used in the main code.  
- **Constants files** → contain constants called in the scripts.  

This project uses:
- **Python** (`.py`) for preprocessing and data organization  
- **MATLAB** (`.m`) for deconvolution analysis with the *Unfold* toolbox  
- **R** (`.Rmd`) for statistical modeling  

---

## Preprocessing

### Step 1 – Markers
`markers/markers_main.py`  
*Analyzes cardiac and respiratory signals and generates markers to classify each stimulus according to the behavioral response and the cardiac/respiratory phase.*

---

### Step 2 – Deconvolution preprocessing
*Separate preprocessing in Python to apply ICA to the raw file and export it to MATLAB. Deconvolution analyses are then performed in MATLAB using the **Unfold** toolbox.*

- `raw_ICA/raw_ICA_main.py` → Compute ICA and apply to raw file before export  
- `raw_export/raw_export_main.py` → Export raw FIF file to MATLAB `.set`  

---

### Step 3 – Evoked preprocessing
*Separate preprocessing to compute the HEP.*

- `epochs/epochs.py` → Epoching  
- `evoked/autoreject_main.py` → Automatic epoch cleaning using **Autoreject**  
- `evoked/evoked_MNE_main.py` → Compute evoked waveforms  

---

## Deconvolution Analysis
*Performed in MATLAB with the **Unfold** toolbox.*

- `unfold_scripts/unfold_base.m` → Deconvolution for simple R-peak marker (baseline)  
- `unfold_scripts/unfold_ana_aware.m` → Deconvolution for aware vs. unaware  
- `unfold_scripts/unfold_cardfit.m` → Deconvolution for systole vs. diastole  

---

## Statistics

- `stats/behavioral/Behav_stats.ipynb` → Creates behavioral dataframe  
- `stats/behavioral/behavioral_GLM_unfold.Rmd` → GLM modeling of behavioral data  
- `stats/HCT analysis/perceivers_stats.ipynb` → Heartbeat counting task descriptives, fig 2b, fig5
- `stats/HCT analysis/ROI_unfold_hbc.Rmd` → Bayes Factor (BF) computation for the mass univariate comparison in high vs low perceivers
- `stats/HCT analysis/ROIs.ipynb` -> compute ROI df
- `stats/figure_3.ipynb` -> necessary computation to reproduce figure 3
- `stats/figure4.ipynb` -> necessary computation to reproduce figure 4
- `stats/stats_helper.py` -> helping functions used throughout statistics jupyter notebooks
  

