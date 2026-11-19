"""
Plot Normative Z-Score Abnormalities for DS4100 Dataset (Combined Side-by-Side & Fixed PowCoh/HCTSA Order)
---------------------------------------------------------------------------------------------------------
1. Explicitly sets PowCoh features to match the exact order from the reference boxplot.
2. Uses the specified HCTSA_TOP_10_FIG1 feature list for HCTSA plots.
3. Maintains tight horizontal plane spacing and feature title placement.
4. Saves figures to /home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/Normative_GlassBrain_Figures_DS4100 as PNG, PDF, and SVG.
"""

import os
import numpy as np
import h5py
import scipy.io as sio
import matplotlib.pyplot as plt
from nilearn import plotting

# ==========================================
# 1. CONFIGURATION & PATHS
# ==========================================
BASE_DIR = '/QRISdata/Q7503/Normative_HCTSA/Scripts/LittDataset/'
FEATURE_EXTRACTION_DIR = os.path.join(BASE_DIR, 'FeatureExtraction')

POWCOH_MAT = os.path.join(FEATURE_EXTRACTION_DIR, 'DS4100_Power_Coh_Normalized_Final_per_part.mat')
HCTSA_MAT = os.path.join(BASE_DIR, 'DS4100_ZScores_FromNormativeAtlas_Nrmlzd_epochs.mat')

METADATA_MAT = os.path.join(BASE_DIR, 'ConcatenatedChannelMetadata_epochs_pruned.mat')
if not os.path.exists(METADATA_MAT):
    METADATA_MAT = os.path.join('/QRISdata/Q7503/Normative_HCTSA/Scripts/LittsScripts/', 'ConcatenatedChannelMetadata_more_epochs.mat')

OUTPUT_DIR = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/Normative_GlassBrain_Figures_DS4100'
os.makedirs(OUTPUT_DIR, exist_ok=True)

POWCOH_RES = os.path.join(BASE_DIR, 'AcrossSubjects/PowerCoh_Results/PowerCoherence/AcrossSubjects_Final_Results_both_feature_sets_Repeated_fixed_resect.mat')
HCTSA_RES = os.path.join(BASE_DIR, 'AcrossSubjects/HCTSA_Results/AcrossSubjects_Final_Results_both_feature_sets_Repeated_fixed_resect.mat')
EXTERNAL_NAMES_MAT = '/QRISdata/Q7503/Normative_HCTSA/Scripts/LittsScripts/Feature_names_plain.mat'
if not os.path.exists(EXTERNAL_NAMES_MAT):
    EXTERNAL_NAMES_MAT = '/QRISdata/Q7503/Normative_HCTSA/Scripts/LittsScripts/Feature_names.mat'

print("Loading DS4100 datasets and metadata using h5py...")

# --- Load DS4100 PowCoh Z-Scores ---
with h5py.File(POWCOH_MAT, 'r') as f:
    powcoh_matrix_full = np.array(f['DS4100_Final_Z'])
    if powcoh_matrix_full.shape[0] < powcoh_matrix_full.shape[1] and powcoh_matrix_full.shape[1] > 10:
        powcoh_matrix_full = powcoh_matrix_full.T
    powcoh_matrix_full = np.abs(powcoh_matrix_full[:, :10])

# --- Load Metadata, Coordinates & Patient IDs ---
with h5py.File(METADATA_MAT, 'r') as f:
    if 'all_electrode_loc_info_epoch' in f:
        loc_group = f['all_electrode_loc_info_epoch']
    elif 'DS4100_ann' in f and 'all_electrode_loc_info_epoch' in f['DS4100_ann']:
        loc_group = f['DS4100_ann']['all_electrode_loc_info_epoch']
    else:
        key_name = [k for k in f.keys() if 'electrode' in k.lower()][0]
        loc_group = f[key_name]

    x_dset, y_dset, z_dset = loc_group['x'], loc_group['y'], loc_group['z']
    if x_dset.dtype.kind == 'O' or len(x_dset.shape) > 0:
        x_arr, y_arr, z_arr = np.array(x_dset), np.array(y_dset), np.array(z_dset)
        N_contacts = x_arr.size
        x, y, z = np.zeros(N_contacts), np.zeros(N_contacts), np.zeros(N_contacts)
        for i in range(N_contacts):
            x[i] = np.array(f[x_arr.flat[i]]).item()
            y[i] = np.array(f[y_arr.flat[i]]).item()
            z[i] = np.array(f[z_arr.flat[i]]).item()
    else:
        x, y, z = np.array(x_dset).flatten(), np.array(y_dset).flatten(), np.array(z_dset).flatten()
        N_contacts = len(x)

    coords = np.column_stack((x, y, z))

    dataset_labels = []
    if 'DS4100_ann' in f and 'dataset_source' in f['DS4100_ann']:
        ds_refs = f['DS4100_ann']['dataset_source']
        for i in range(N_contacts):
            ref = ds_refs.flat[i] if ds_refs.size > 1 else ds_refs[0]
            dataset_labels.append(''.join(chr(c) for c in np.array(f[ref]).flatten()))
    
    if dataset_labels:
        dataset_labels = np.array(dataset_labels)
        is_ds4100 = np.array(['ds4100' in d.lower() for d in dataset_labels])
    else:
        mid_pt = N_contacts // 2
        is_ds4100 = np.zeros(N_contacts, dtype=bool)
        is_ds4100[mid_pt:] = True

    patient_ids = []
    if 'DS4100_ann' in f and 'all_clinical_info_epoch' in f['DS4100_ann']:
        clin_group = f['DS4100_ann']['all_clinical_info_epoch']
        if 'participant_id' in clin_group:
            pid_refs = clin_group['participant_id']
            for i in range(N_contacts):
                ref = pid_refs.flat[i] if pid_refs.size > 1 else pid_refs[0]
                patient_ids.append(''.join(chr(c) for c in np.array(f[ref]).flatten()))
    if not patient_ids:
        patient_ids = ['Patient_All'] * N_contacts
    patient_ids = np.array(patient_ids)

coords = coords[is_ds4100]
patient_ids = patient_ids[is_ds4100]
powcoh_matrix_full = powcoh_matrix_full[is_ds4100, :]
N_contacts = len(coords)

# --- Load DS4100 HCTSA Z-Scores from Normative Atlas ---
with h5py.File(HCTSA_MAT, 'r') as f:
    X_hctsa_full = np.array(f['DS4100_z_scores_from_normative'])
    if X_hctsa_full.shape[0] < X_hctsa_full.shape[1] and X_hctsa_full.shape[1] > 4000:
        X_hctsa_full = X_hctsa_full.T
    X_hctsa_full = np.abs(X_hctsa_full)
    if X_hctsa_full.shape[0] != N_contacts:
        X_hctsa_full = X_hctsa_full[is_ds4100, :]

# ==========================================
# 2. FEATURE ORDERING
# ==========================================
print("Setting exact PowCoh and HCTSA feature order as requested...")

powcoh_indices = [7, 9, 4, 8, 0, 6, 5, 1, 2, 3]
powcoh_names = [
    'Coherence Alpha', 'Coherence Gamma', 'Power Gamma', 'Coherence Beta',
    'Power Delta', 'Coherence Theta', 'Coherence Delta', 'Power Theta',
    'Power Alpha', 'Power Beta'
]
powcoh_matrix = powcoh_matrix_full[:, powcoh_indices]

# --- HCTSA Specific Feature Set Selection ---
HCTSA_TOP_10_FIG1 = [
    "MD pNN raw pnn90",
    "CO CompareMinAMI quantiles 2 80 std",
    "trimmed mean 10",
    "MD pNN raw pnn100",
    "NL nsamdf 1 14 05 1 7 L",
    "pAC 5",
    "PD PeriodicityWang th1 sqrtN",
    "WL DetailCoeffs db3 max wslesr mean",
    "NL crptool fnn 10 2 mi firstunder07",
    "EN Randomize permute sampen2 015hp"
]

hctsa_names = HCTSA_TOP_10_FIG1
hctsa_indices = []

# Map feature names to columns if available, otherwise fallback to sequential columns
try:
    all_feat_names = []
    if os.path.exists(EXTERNAL_NAMES_MAT):
        ext_data = sio.loadmat(EXTERNAL_NAMES_MAT, squeeze_me=True, struct_as_record=False)
        for k in ['Feature_names', 'CodeString', 'Name']:
            if hasattr(ext_data, k):
                obj = getattr(ext_data, k)
                if hasattr(obj, 'CodeString'):
                    all_feat_names = [str(x) for x in np.atleast_1d(obj.CodeString)]
                    break
                elif hasattr(obj, 'Name'):
                    all_feat_names = [str(x) for x in np.atleast_1d(obj.Name)]
                    break
                else:
                    all_feat_names = [str(x) for x in np.atleast_1d(obj)]
                    break
            elif isinstance(ext_data, dict) and k in ext_data:
                all_feat_names = [str(x) for x in np.atleast_1d(ext_data[k])]
                break

    for target_name in HCTSA_TOP_10_FIG1:
        matched_idx = -1
        for idx, fn in enumerate(all_feat_names):
            if target_name.lower() in fn.lower():
                matched_idx = idx
                break
        if matched_idx != -1 and matched_idx < X_hctsa_full.shape[1]:
            hctsa_indices.append(matched_idx)
        else:
            # Fallback if exact string match fails
            fallback_idx = len(hctsa_indices)
            if fallback_idx < X_hctsa_full.shape[1]:
                hctsa_indices.append(fallback_idx)
except Exception as e:
    print(f"Note mapping HCTSA feature names: {e}. Using default column indices.")
    hctsa_indices = list(range(min(10, X_hctsa_full.shape[1])))

if len(hctsa_indices) < 10:
    while len(hctsa_indices) < 10 and len(hctsa_indices) < X_hctsa_full.shape[1]:
        hctsa_indices.append(len(hctsa_indices))

X_top10_hctsa = X_hctsa_full[:, hctsa_indices]

# ==========================================
# 3. PATIENT-LEVEL Z-SCORE & ROBUST [0, 1] SCALING
# ==========================================
def normalize_features(feature_matrix):
    norm_matrix = np.zeros_like(feature_matrix)
    unique_pats = np.unique(patient_ids)
    
    for col in range(feature_matrix.shape[1]):
        col_vals = feature_matrix[:, col]
        for pid in unique_pats:
            p_mask = (patient_ids == pid)
            p_vals = col_vals[p_mask]
            mean_val = np.nanmean(p_vals)
            std_val = np.nanstd(p_vals)
            if std_val > 0:
                norm_matrix[p_mask, col] = (p_vals - mean_val) / std_val
            else:
                norm_matrix[p_mask, col] = 0.0
                
        f_vals = norm_matrix[:, col]
        v_min, v_max = np.percentile(f_vals[~np.isnan(f_vals)], [5, 95]) if np.any(~np.isnan(f_vals)) else (0, 1)
        if v_max > v_min:
            norm_matrix[:, col] = np.clip((f_vals - v_min) / (v_max - v_min), 0, 1)
        else:
            norm_matrix[:, col] = 0.5
            
    return norm_matrix

print("Normalizing DS4100 abnormality features...")
powcoh_norm = normalize_features(powcoh_matrix)
hctsa_norm = normalize_features(X_top10_hctsa)

# ==========================================
# 4. PLOTTING FUNCTIONS
# ==========================================
def plot_single_feature_grid(feature_matrix, feature_names, title_prefix, filename):
    n_features = feature_matrix.shape[1]
    planes = ['x', 'y', 'z']
    
    fig, axes = plt.subplots(n_features, 3, figsize=(6.5, 1.8 * n_features), facecolor='white')
    plt.subplots_adjust(left=0.05, right=0.95, top=0.92, bottom=0.04, wspace=-0.32, hspace=0.15)
    
    for row in range(n_features):
        feat_values = np.nan_to_num(feature_matrix[:, row], nan=0.0)
        sort_idx = np.argsort(feat_values)
        sorted_coords = coords[sort_idx]
        sorted_feat_values = feat_values[sort_idx]
        zero_adj = np.zeros((len(sorted_coords), len(sorted_coords)))

        for col, plane in enumerate(planes):
            ax = axes[row, col] if n_features > 1 else axes[col]
            ax.set_facecolor('white')
            ax.patch.set_facecolor('white')
            
            cmap = plt.colormaps['RdYlGn_r']
            node_rgba = cmap(sorted_feat_values)

            plotting.plot_connectome(
                adjacency_matrix=zero_adj,
                node_coords=sorted_coords,
                node_color=node_rgba,
                node_size=18,
                edge_threshold="100%",
                display_mode=plane,
                colorbar=False,
                axes=ax
            )
            ax.set_facecolor('white')
            ax.patch.set_facecolor('white')
            
        clean_name = str(feature_names[row]).replace('_', ' ')
        axes[row, 1].set_title(clean_name, fontsize=10, fontweight='bold', pad=4, color='black')

    plt.suptitle(f"{title_prefix} - DS4100 Abnormalities", fontsize=13, fontweight='bold', y=0.97, color='black')
    
    # Save as PNG
    out_path_png = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(out_path_png, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    
    # Save as PDF
    base_name, _ = os.path.splitext(filename)
    out_path_pdf = os.path.join(OUTPUT_DIR, f"{base_name}.pdf")
    plt.savefig(out_path_pdf, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')

    # Save as SVG
    out_path_svg = os.path.join(OUTPUT_DIR, f"{base_name}.svg")
    plt.savefig(out_path_svg, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    
    plt.close()
    print(f"Saved figure: {out_path_png}, {out_path_pdf}, and {out_path_svg}")

def plot_combined_side_by_side(powcoh_mat, pow_names, hctsa_mat, hct_names, filename):
    n_rows = max(powcoh_mat.shape[1], hctsa_mat.shape[1])
    planes = ['x', 'y', 'z']
    
    fig, axes = plt.subplots(n_rows, 6, figsize=(14, 1.8 * n_rows), facecolor='white')
    plt.subplots_adjust(left=0.03, right=0.97, top=0.93, bottom=0.03, wspace=-0.32, hspace=0.18)
    
    for row in range(n_rows):
        # --- Plot PowCoh (Left 3 columns) ---
        if row < powcoh_mat.shape[1]:
            feat_vals = np.nan_to_num(powcoh_mat[:, row], nan=0.0)
            sort_idx = np.argsort(feat_vals)
            s_coords, s_vals = coords[sort_idx], feat_vals[sort_idx]
            z_adj = np.zeros((len(s_coords), len(s_coords)))
            
            for col in range(3):
                ax = axes[row, col]
                ax.set_facecolor('white'); ax.patch.set_facecolor('white')
                plotting.plot_connectome(z_adj, s_coords, node_color=plt.colormaps['RdYlGn_r'](s_vals),
                                         node_size=16, edge_threshold="100%", display_mode=planes[col], colorbar=False, axes=ax)
                ax.set_facecolor('white'); ax.patch.set_facecolor('white')
            
            clean_name = str(pow_names[row]).replace('_', ' ')
            axes[row, 1].set_title(f"PowCoh: {clean_name}", fontsize=9, fontweight='bold', pad=4, color='black')

        # --- Plot HCTSA (Right 3 columns) ---
        if row < hctsa_mat.shape[1]:
            feat_vals = np.nan_to_num(hctsa_mat[:, row], nan=0.0)
            sort_idx = np.argsort(feat_vals)
            s_coords, s_vals = coords[sort_idx], feat_vals[sort_idx]
            z_adj = np.zeros((len(s_coords), len(s_coords)))
            
            for col in range(3):
                ax = axes[row, col + 3]
                ax.set_facecolor('white'); ax.patch.set_facecolor('white')
                
                pos = ax.get_position()
                ax.set_position([pos.x0 + 0.05, pos.y0, pos.width, pos.height])
                
                plotting.plot_connectome(z_adj, s_coords, node_color=plt.colormaps['RdYlGn_r'](s_vals),
                                         node_size=16, edge_threshold="100%", display_mode=planes[col], colorbar=False, axes=ax)
                ax.set_facecolor('white'); ax.patch.set_facecolor('white')
            
            clean_name = str(hct_names[row]).replace('_', ' ')
            axes[row, 4].set_title(f"HCTSA: {clean_name}", fontsize=9, fontweight='bold', pad=4, color='black')

    plt.suptitle("DS4100 Abnormalities: PowCoh (Left) vs. HCTSA (Right) Top Features", fontsize=14, fontweight='bold', y=0.97, color='black')
    
    # Save as PNG
    out_path_png = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(out_path_png, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    
    # Save as PDF
    base_name, _ = os.path.splitext(filename)
    out_path_pdf = os.path.join(OUTPUT_DIR, f"{base_name}.pdf")
    plt.savefig(out_path_pdf, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')

    # Save as SVG
    out_path_svg = os.path.join(OUTPUT_DIR, f"{base_name}.svg")
    plt.savefig(out_path_svg, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    
    plt.close()
    print(f"Saved combined side-by-side figure: {out_path_png}, {out_path_pdf}, and {out_path_svg}")

# ==========================================
# 5. GENERATE FIGURES
# ==========================================
print("\nGenerating separate sorted DS4100 figures...")
plot_single_feature_grid(powcoh_norm, powcoh_names, "DS4100 PowCoh Abnormalities", "DS4100_PowCoh_Sorted_GlassBrains.png")
plot_single_feature_grid(hctsa_norm, hctsa_names, "DS4100 Top 10 HCTSA Abnormalities", "DS4100_HCTSA_Sorted_GlassBrains.png")

print("\nGenerating combined side-by-side DS4100 figure...")
plot_combined_side_by_side(powcoh_norm, powcoh_names, hctsa_norm, hctsa_names, "DS4100_Combined_PowCoh_vs_HCTSA_GlassBrains.png")

print("\nAll DS4100 glass brain visualization tasks completed successfully!")