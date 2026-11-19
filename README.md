# Normative HCTSA & Feature Extraction Pipeline for Surgical Outcome Prediction in Epilepsy

This repository contains the end-to-end data processing, feature extraction, interactive pruning, normative atlas comparisons, glass-brain and exemplar signal visualization, and leave-one-subject-out (LOSO) machine learning classification pipeline for evaluating electrophysiological predictors of surgical treatment outcome in patients with focal epilepsy.

## 📌 Repository Overview

The codebase processes multi-centre intracranial electroencephalography (iEEG/sEEG/ECoG) data through a standardized 10-step workflow:

1. Preprocessing & Bipolar GM Localization (C0, C1): Resampling, anti-alias filtering, power-line interference reduction, MNI gray matter coordinate mapping, distance-constrained bipolar derivation (≤10 mm), 1-minute epoching, and clinical metadata merging.
2. HCTSA Data Formatting & Interactive Pruning (C2, C3): Formatting time series for hctsa analysis and interactive visual artifact screening.
3. High-Dimensional Feature Computation (C4, G1_*): Computing 7,770 non-linear time-series operations using hctsa, 10-band spectral power/coherence (PowCoh), and automated HFO-spike detection.
4. Metadata Alignment (C5, C6): Synchronizing spatial coordinates, channel status (SOZ, resected, spared), and clinical attributes post-pruning.
5. Anatomical Mapping & Normative Atlas Comparison (G2_*, G3_*): Mapping contact coordinates to brain atlas regions and z-scoring patient feature distributions against healthy normative baselines (Frauscher atlas).
6. Cohort Feature Abnormality & Distribution Plotting (G4_*): Evaluating contact-level feature abnormality differentials (Resected SOZ vs. Spared/Baseline) across HCTSA, PowCoh, and HFO-spike domains.
7. Publication Glass Brain Visualization (brains_figure.py, brains_figure_all_subjects_abnormalities_HCTA_top_families.py): Generating horizontal multi-view glass-brain anatomical projections and multi-feature abnormality grids using standard brain templates.
8. Individual & Cohort Abnormality Topology with Bayes Factors (brains_figure_abnormality_indiv_subjects.py): Rendering combined PowCoh/HCTSA glass brains, resected-vs.-spared contact abnormality scatter plots with one-sample Savage-Dickey Bayes factors, and cohort-level summaries across surgical outcome groups.
9. Clinical Covariate Analysis (G6_Clinical_vars_cors.m): Computing statistical correlation matrices across demographic, diagnostic, and surgical variables.
10. LOSO Classification, Feature Importance & Exemplar Trace Visualization (G8_*, G5, Simulated_signals_plotting.py): Leave-One-Subject-Out Cross-Validation using Random Forest ensembles across Clinical, PowCoh, and HCTSA feature spaces, feature importance heatmaps, and synthetic exemplar time-series illustrations for top features.

## ⚙️ Environment & Dependencies

### Environments
- MATLAB R2021a or newer (Parallel Computing Toolbox recommended for parpool acceleration)
- Python 3.8+ (for glass-brain and exemplar signal plotting scripts)

### Required Toolboxes & External Libraries
- FieldTrip Toolbox (e.g., fieldtrip-20250523) for signal processing and EDF I/O
- hctsa Toolbox: Highly Comparative Time-Series Analysis framework (hctsa-main)
- MNI to Atlas Converter (mni2atlas) and NIfTI Toolbox (NIfTI_20140122) for anatomical localization
- Spike-Ripple Detector Toolbox (Spike-Ripple-Detector-Method-master)
- Python Libraries: h5py, nilearn, matplotlib, numpy, scipy

## 📂 File Structure & Pipeline Scripts

.
├── C0_PreProcess_MultiCentre_dataset_gray_BiP_more_epochs.m
├── C1_adding_clinical_discrpt_BiPo_more_epochs.m
├── C2_Loading_concatnting_data_epochs_ds4100_for_HCTSA.m
├── C3_Interactive_pruning_data.m
├── C4_Extracting_hcta_features.m
├── C5_Loading_concatnting_meta_epchs_ds4100_for_after_HCTSA.m
├── C6_Fixing_concatenated_meta_data_upon_pruning.m
├── G1_Extracting_HFO_spiking_Features_fixed_per_part.m
├── G1_Extracting_power_and_coherence_Features_fixed_per_part.m
├── G2_Mapping_contacts_to_atlas_regions_ZscoringDisbald.m
├── G2_Mapping_contacts_to_atlas_regions_ZscoringDisbald_PowCoh.m
├── G3_Comparing_the_DS4100_to_Frauscher_Nrmlzd_Zscore_per.m
├── G3_Comparing_the_DS4100_to_Frauscher_Nrmlzd_Zscore_per_PowCoh.m
├── G4_Classif_AcrsSubj_par_feat_imprtc_HFOSpk_resct.m
├── G4_Classif_AcrsSubj_par_feat_imprtc_HFOSpk_plot.m
├── G4_Classif_AcrsSubj_par_feat_imprtc_both_features.m
├── G4_Classif_AcrsSubj_par_feat_imprtc_both_features_plot.m
├── G4_Classif_AcrsSubj_par_feat_imprtc_both_features_cors_plt.m
├── G4_Classif_AcrsSubj_par_feat_imprtc_and_dist_plt.m
├── G4_Classif_AcrsSubj_par_feat_imprtc_HFOSpk_plot.m
├── G5_plotting_feature_importances_tmp.m
├── G6_Clinical_vars_cors.m
├── G8_Pred_All_Vars_from_clin_PowCoh_HCTSA_upsm_SZNon_feature_dist.m
├── G8_Predicting_All_Vars_from_clin_PowCoh_HCTSA_ups_SZNLOPO_cmplt.m
├── G8_Predicting_All_Vars_from_clin_PowCoh_HCTSA_SZNLOPO_cmplt_plt.m
├── Simulated_signals_plotting.py
├── brains_figure.py
├── brains_figure_abnormality_indiv_subjects.py
├── brains_figure_all_subjects_abnormalities_HCTA_top_families.py
└── README.md

## 🚀 Execution Workflow

### Step 1: Preprocessing & Bipolar Channel Derivation
- Script: C0_PreProcess_MultiCentre_dataset_gray_BiP_more_epochs.m
- Input: Raw .edf files, BIDS _channels.tsv, _electrodes.tsv.
- Functionality:
  - Resamples signals to 200 Hz with an anti-aliasing low-pass filter at 80 Hz.
  - Applies DFT notch filters (60, 120, 180 Hz) and an FIR high-pass filter at 0.5 Hz.
  - Maps MNI coordinates to the Harvard-Oxford atlas to exclude white matter contacts.
  - Forms bipolar channels between adjacent gray matter contacts within ≤10 mm Euclidean distance.
  - Segments continuous data into 60-second epochs.
- Output: HUP<ID>_SEEG_InterIctal_Processed_Trial<XX>_Epoch<YY>_Bipolar_GM_only.mat.

### Step 2: Clinical Metadata Merging
- Script: C1_adding_clinical_discrpt_BiPo_more_epochs.m
- Input: BIDS participants.tsv and epoched .mat files from Step 1.
- Functionality: Appends patient demographic data (age, sex, handedness, age of onset, disease duration, surgical target, lesion status, Engel outcome) to each channel structure.

### Step 3: HCTSA Data Structuring & Interactive Pruning
- Scripts: C2_Loading_concatnting_data_epochs_ds4100_for_HCTSA.m & C3_Interactive_pruning_data.m
- Functionality: Aggregates epoched signals across patients into a unified timeSeriesData cell array formatted for hctsa, followed by interactive GUI-based artifact screening.
- Output: HCTSA_AllPreprocessedSignals_epochs_pruned.mat.

### Step 4: Multi-Domain Feature Extraction
- Script 4A (hctsa): C4_Extracting_hcta_features.m
  - Computes 7,770 non-linear time-series operations using TS_Init and TS_Compute.
  - Output: HCTSA.mat.
- Script 4B (HFO-Spikes): G1_Extracting_HFO_spiking_Features_fixed_per_part.m
  - Runs automated spike-ripple detection (spike_ripple_detector) to extract interictal discharge rates and spike timestamps per contact.
  - Output: DS4100_HFO_Spikes_Features.mat.
- Script 4C (Power & Coherence): G1_Extracting_power_and_coherence_Features_fixed_per_part.m
  - Computes absolute spectral power and pairwise mean coherence across 5 standard bands (δ, θ, α, β, γ).
  - Output: Data_AllContacts_BothDatasets_epochs_PowCoh_per_part.mat.

### Step 5: Metadata Synchronization & Normative Atlas Mapping
- Metadata Alignment: C5_Loading_concatnting_meta_epchs_ds4100_for_after_HCTSA.m & C6_Fixing_concatenated_meta_data_upon_pruning.m
- Atlas Mapping & Normative Z-Scoring:
  - HCTSA Domain: G2_Mapping_contacts_to_atlas_regions_ZscoringDisbald.m & G3_Comparing_the_DS4100_to_Frauscher_Nrmlzd_Zscore_per.m
  - PowCoh Domain: G2_Mapping_contacts_to_atlas_regions_ZscoringDisbald_PowCoh.m & G3_Comparing_the_DS4100_to_Frauscher_Nrmlzd_Zscore_per_PowCoh.m
- Functionality: Maps bipolar contact locations to standardized brain parcellations and computes per-contact feature abnormality z-scores relative to healthy normative control baselines (Frauscher atlas).

### Step 6: Contact-Level Abnormality Evaluation & Visualizations (G4 Series)
- Core Computational Scripts:
  - G4_Classif_AcrsSubj_par_feat_imprtc_both_features.m: Evaluates multivariate HCTSA + PowCoh feature abnormalities between resected SOZ contacts (resect,soz) and spared contacts (soz, resect, n/a).
  - G4_Classif_AcrsSubj_par_feat_imprtc_HFOSpk_resct.m: Evaluates contact-level abnormality profiles specifically for HFO-spike rates.
- Plotting & Analysis Variants:
  - G4_Classif_AcrsSubj_par_feat_imprtc_both_features_plot.m: Generates individual and cohort-level feature abnormality comparison plots.
  - G4_Classif_AcrsSubj_par_feat_imprtc_both_features_cors_plt.m: Plots correlation structure between HCTSA and PowCoh abnormality metrics across contacts.
  - G4_Classif_AcrsSubj_par_feat_imprtc_and_dist_plt.m: Plots contact abnormality distributions as a function of Euclidean distance from the resection zone boundary.
  - G4_Classif_AcrsSubj_par_feat_imprtc_HFOSpk_plot.m: Plots HFO-spike rate abnormality distributions across outcome groups.

### Step 7: Anatomical & Top-Feature Glass Brain Visualization
- Scripts: brains_figure.py & brains_figure_all_subjects_abnormalities_HCTA_top_families.py
- Functionality: Python scripts using nilearn and matplotlib to render glass-brain visualizations (RdYlGn_r colormap) in x, y, and z projections. brains_figure_all_subjects_abnormalities_HCTA_top_families.py normalizes patient-level normative z-scores with 5th-95th percentile scaling and generates separate and side-by-side multi-row grids for top PowCoh and HCTSA features across the DS4100 cohort, exporting PNG, PDF, and SVG files.

### Step 8: Individual & Cohort Abnormality Topology with Bayes Factors
- Script: brains_figure_abnormality_indiv_subjects.py
- Functionality: Loads HCTSA (DS4100 vs. Frauscher) and PowCoh HDF5 datasets alongside LOSO classification results to identify top features, compute cohort-level Mann-Whitney U and Spearman summary statistics for normalized Δ|Z| (resected vs. spared), and render per-subject combined PowCoh/HCTSA glass brains, resected-vs.-spared contact abnormality scatter plots annotated with one-sample Savage-Dickey Bayes factors (BF10), and pooled Success vs. Failure cohort distribution figures (saved as PNG, PDF, and SVG).

### Step 9: Clinical Covariates & Feature Importance Analysis
- Clinical Correlations: G6_Clinical_vars_cors.m
  - Computes correlation matrices and independence tests across clinical predictors (age, age of onset, disease duration, lesion status, surgical target, therapy).
- Feature Importance Heatmaps: G5_plotting_feature_importances_tmp.m
  - Visualizes top-ranked Random Forest feature importances across subjects as normalized 2D heatmaps with marginal average traces.

### Step 10: LOSO Machine Learning, Outcome Prediction & Exemplar Signal Traces
- Core Classification Script: G8_Predicting_All_Vars_from_clin_PowCoh_HCTSA_ups_SZNLOPO_cmplt.m
  - Trains Random Forest models (TreeBagger, 150 trees) under Leave-One-Subject-Out (LOSO) cross-validation with training fold upsampling.
  - Evaluates 7 model architectures (Clin_Only, HCTSA_Only, PowCoh_Only, HCTSA_PowCoh, Clin_HCTSA, Clin_PowCoh, All_Features) across outcome targets (outcome, lesion_status, target).
- Variant & Plotting Scripts:
  - G8_Pred_All_Vars_from_clin_PowCoh_HCTSA_upsm_SZNon_feature_dist.m: Similar LOSO pipeline variant evaluating feature distributions between SOZ and Non-SOZ contacts.
  - G8_Predicting_All_Vars_from_clin_PowCoh_HCTSA_SZNLOPO_cmplt_plt.m: Generates summary ROC curves, accuracy bar charts, and variance plots from completed LOSO MAT results files.
- Exemplar Time-Series Illustrations (Simulated_signals_plotting.py): Synthesizes normal (green) versus abnormal (red) exemplar EEG/iEEG traces and proxy feature values for two 5-row HCTSA feature sets and the top 10 PowCoh features, exporting figures as PDF, PNG, and SVG.

## 📊 Summary of Model Predictor Domains

| Model Name | Included Feature Sets |
| :--- | :--- |
| Clin_Only | Age, Age of Onset, Duration, Sex, Handedness, Implant, Target, Lesion Status |
| HCTSA_Only | Top 10 PCs derived from 7,770 time-series operations |
| PowCoh_Only | Absolute spectral power & coherence across 5 frequency bands (δ, θ, α, β, γ) |
| Clin_HCTSA | Clinical predictors + HCTSA 10 PCs |
| Clin_PowCoh | Clinical predictors + PowCoh features |
| HCTSA_PowCoh | Integrated HCTSA 10 PCs + PowCoh features |
| All_Features | Integrated Clinical + PowCoh + HCTSA feature space |
```[cite: 4, 5, 6]
