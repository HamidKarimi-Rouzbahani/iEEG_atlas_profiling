# Normative HCTSA & Feature Extraction Pipeline for Surgical Outcome Prediction in Epilepsy

This repository contains the end-to-end data processing, feature extraction, interactive pruning, normative atlas comparisons, and leave-one-subject-out (LOSO) machine learning classification pipeline for evaluating electrophysiological predictors of surgical treatment outcome in patients with focal epilepsy.

---

## 📌 Repository Overview

The codebase processes multi-centre intracranial electroencephalography (iEEG/sEEG/ECoG) data through a standardized 5-step workflow:

1. **Preprocessing & Bipolar GM Localization** (`C0`, `C1`): Resampling, anti-alias filtering, power-line interference reduction, MNI gray matter coordinate mapping, distance-constrained bipolar derivation ($\le 10\text{ mm}$), 1-minute epoching, and clinical metadata merging.
2. **HCTSA Data Formatting & Interactive Pruning** (`C2`, `C3`): Formatting time series for `hctsa` analysis and interactive visual artifact screening.
3. **High-Dimensional Feature Computation** (`C4`, `G1`): Computing $7,770$ non-linear time-series operations using `hctsa` alongside automated HFO-spike detection.
4. **Metadata Alignment** (`C5`, `C6`): Synchronizing spatial coordinates, channel status (SOZ, resected, spared), and clinical attributes post-pruning.
5. **LOSO Classification & Outcome Modeling** (`G8`, `G5`): Leave-One-Subject-Out Cross-Validation using Random Forest ensembles across Clinical, PowCoh, and HCTSA feature spaces.

---

## ⚙️ Environment & Dependencies

### MATLAB Requirements
* **MATLAB** R2021a or newer
* **Parallel Computing Toolbox** (for `parpool` acceleration)

### Required Toolboxes & External Libraries
* **FieldTrip Toolbox** (e.g., `fieldtrip-20250523` or compatible) for signal processing and EDF I/O
* **`hctsa` Toolbox**: Highly Comparative Time-Series Analysis framework (`hctsa-main`)
* **MNI to Atlas Converter** (`mni2atlas`) and **NIfTI Toolbox** (`NIfTI_20140122`) for anatomical localization
* **Spike-Ripple Detector Toolbox** (`Spike-Ripple-Detector-Method-master`)

---

## 📂 File Structure & Execution Pipeline

.
├── C0_PreProcess_MultiCentre_dataset_gray_BiP_more_epochs.m
├── C1_adding_clinical_discrpt_BiPo_more_epochs.m
├── C2_Loading_concatnting_data_epochs_ds4100_for_HCTSA.m
├── C3_Interactive_pruning_data.m
├── C4_Extracting_hcta_features.m
├── C5_Loading_concatnting_meta_epchs_ds4100_for_after_HCTSA.m
├── C6_Fixing_concatenated_meta_data_upon_pruning.m
├── G1_Extracting_HFO_spiking_Features_fixed_per_part.m
├── G5_plotting_feature_importances_tmp.m
├── G8_Predicting_All_Vars_from_clin_PowCoh_HCTSA_ups_SZNLOPO_cmplt.m
└── README.md


---

## 🚀 Execution Workflow

### Step 1: Preprocessing & Bipolar Channel Derivation
* **Script:** `C0_PreProcess_MultiCentre_dataset_gray_BiP_more_epochs.m`
* **Input:** Raw `.edf` files, BIDS `_channels.tsv`, `_electrodes.tsv`.
* **Functionality:**
  * Resamples signals to $200\text{ Hz}$ with an anti-aliasing low-pass filter at $80\text{ Hz}$.
  * Applies DFT notch filters at $60\text{ Hz}$, $120\text{ Hz}$, and $180\text{ Hz}$ and an FIR high-pass filter at $0.5\text{ Hz}$.
  * Maps MNI coordinates to the Harvard-Oxford atlas to exclude white matter contacts.
  * Forms bipolar channels between adjacent gray matter contacts within $\le 10\text{ mm}$ Euclidean distance.
  * Segments continuous data into $60$-second epochs.
* **Output:** `HUP<ID>_SEEG_InterIctal_Processed_Trial<XX>_Epoch<YY>_Bipolar_GM_only.mat`.

### Step 2: Clinical Metadata Merging
* **Script:** `C1_adding_clinical_discrpt_BiPo_more_epochs.m`
* **Input:** BIDS `participants.tsv` and epoched `.mat` files from Step 1.
* **Functionality:** Appends patient demographic data (age, sex, handedness, age of onset, disease duration, surgical target, lesion status, Engel outcome) to each electrode channel structure.

### Step 3: HCTSA Data Structuring
* **Script:** `C2_Loading_concatnting_data_epochs_ds4100_for_HCTSA.m`
* **Functionality:** Aggregates all epoched signals across patients, trials, and channels into a unified `timeSeriesData` cell array formatted for `hctsa`.
* **Output:** `HCTSA_AllPreprocessedSignals_epochs.mat`.

### Step 4: Interactive Artifact Pruning (Optional / QA)
* **Script:** `C3_Interactive_pruning_data.m`
* **Functionality:** GUI-based batch visualization allowing manual inspection and removal of noisy/corrupted channels.
* **Output:** `HCTSA_AllPreprocessedSignals_epochs_pruned.mat`.

### Step 5: Feature Extraction
* **Script A (`hctsa`):** `C4_Extracting_hcta_features.m`
  * Initializes and computes $7,770$ non-linear features (`TS_Init`, `TS_Compute`).
  * **Output:** `HCTSA.mat`.
* **Script B (HFO-Spikes):** `G1_Extracting_HFO_spiking_Features_fixed_per_part.m`
  * Runs automated spike-ripple detection across channels.
  * **Output:** `DS4100_HFO_Spikes_Features.mat`.

### Step 6: Metadata Synchronization
* **Scripts:** `C5_Loading_concatnting_meta_epchs_ds4100_for_after_HCTSA.m` & `C6_Fixing_concatenated_meta_data_upon_pruning.m`
* **Functionality:** Concatenates channel metadata across all subjects and updates channel indices to match pruned time series.
* **Output:** `ConcatenatedChannelMetadata_epochs_pruned.mat`.

### Step 7: LOSO Machine Learning & Outcome Prediction
* **Script:** `G8_Predicting_All_Vars_from_clin_PowCoh_HCTSA_ups_SZNLOPO_cmplt.m`
* **Methodology:**
  * Computes regional abnormality profiles ($\text{SOZ} - \text{Non-SOZ}$) per patient.
  * Reduces $7,770$ HCTSA features to 10 principal components via PCA per fold.
  * Evaluates 7 model architectures (`Clin_Only`, `HCTSA_Only`, `PowCoh_Only`, `HCTSA_PowCoh`, `Clin_HCTSA`, `Clin_PowCoh`, `All_Features`) using Random Forest classifiers (`TreeBagger`, 150 trees) under Leave-One-Subject-Out (LOSO) cross-validation.
  * Projects PCA out-of-bag feature importances back to raw feature space.
* **Output:** Classification metrics, ROC curves, and MAT files in `Classification_Results_LOSO_cmplt/`.

### Step 8: Visualization & Feature Importance Analysis
* **Script:** `G5_plotting_feature_importances_tmp.m`
* **Functionality:** Generates heatmap matrices and marginal traces for top-ranked features across within-subject and across-subject models.

---

## 📊 Summary of Model Predictor Domains

| Model Name | Included Feature Sets |
| :--- | :--- |
| **`Clin_Only`** | Age, Age of Onset, Duration, Sex, Handedness, Implant, Target, Lesion Status |
| **`HCTSA_Only`** | Top 10 PCs derived from 7,770 time-series operations |
| **`PowCoh_Only`** | Absolute spectral power & coherence across 5 frequency bands ($\delta, \theta, \alpha, \beta, \gamma$) |
| **`Clin_HCTSA`** | Clinical predictors + HCTSA 10 PCs |
| **`All_Features`** | Integrated Clinical + PowCoh + HCTSA feature space |

---

## ✉️ Contact & Citation
For questions regarding dataset access, execution bugs, or manuscript replication details, please open an issue in this repository.
