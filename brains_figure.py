import os
import numpy as np
import h5py
import scipy.io as sio
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from nilearn import plotting

# ==============================================================================
# Helper Function: Compute One-Sample Bayes Factor (Savage-Dickey Ratio)
# ==============================================================================
def compute_bayes_factor_one_sample(diff_array):
    diff_array = diff_array[~np.isnan(diff_array)]
    N = len(diff_array)
    if N < 2 or np.std(diff_array) == 0:
        return 1.0
    t_stat = np.mean(diff_array) / (np.std(diff_array, ddof=1) / np.sqrt(N))
    r = np.sqrt(2) / 2.0  # Default JZS Cauchy prior scale factor
    
    def integrand(g):
        return (1.0 + N * g * (r**2))**(-0.5) * \
               (1.0 + (t_stat**2) / ((1.0 + N * g * (r**2)) * (N - 1)))**(-N / 2.0) * \
               (g**(-1.5)) * np.exp(-1.0 / (2.0 * g))
               
    try:
        from scipy.integrate import quad
        marg_lik_h1, _ = quad(integrand, 0, np.inf)
    except Exception:
        marg_lik_h1 = 1e-300
        
    marg_lik_h0 = (1.0 + (t_stat**2) / (N - 1))**(-N / 2.0)
    bf10 = marg_lik_h1 / marg_lik_h0
    return max(bf10, 1e-300)

# ==============================================================================
# 1. Setup Paths & Load Datasets
# ==============================================================================
base_folder         = '/scratch/user/uqhkarim/Normative_HCTSA/Scripts/LittDataset/'
data_path           = os.path.join(base_folder, 'Data_AllContacts_BothDatasets_epochs.mat')
powcoh_path         = os.path.join(base_folder, 'FeatureExtraction/Data_AllContacts_BothDatasets_epochs_PowCoh_per_part.mat')
external_names_file = '/home/uqhkarim/Normative_HCTSA/Scripts/LittsScripts/Feature_names.mat'
results_mat_file    = os.path.join(base_folder, 'Classification_Results_LOSO_cmplt/outcome/Classification_outcome_LOSO_Results_complete.mat')
output_fig_dir      = os.path.join(base_folder, 'SinglePanel_GlassBrains_Figures/')
os.makedirs(output_fig_dir, exist_ok=True)

print("Loading Master HDF5 dataset...")
with h5py.File(data_path, 'r') as f:
    ann = f['DS4100_annotations']
    
    # Load HCTSA Raw Matrices
    DS4100_raw    = np.array(f['DS4100'])        # Shape: (7770, 11380) or (11380, 7770)
    Frauscher_raw = np.array(f['Frauscher'])    # Shape: (7770, 1772) or (1772, 7770)
    
    # Standardize HCTSA orientation to: (11380_epochs, 7770_features)
    DS4100_mat_hctsa = DS4100_raw.T if DS4100_raw.shape[0] == 7770 else DS4100_raw
    Frauscher_mat_hctsa = Frauscher_raw.T if Frauscher_raw.shape[0] == 7770 else Frauscher_raw

    # Dereference Electrode Coordinates (x, y, z)
    loc_group = ann['all_electrode_loc_info_epoch']
    def extract_coord(key):
        refs = np.array(loc_group[key]).squeeze()
        return np.array([float(np.array(f[ref][()] if isinstance(ref, h5py.h5r.Reference) else ref).squeeze()) for ref in refs])
        
    coords_epoch = np.column_stack([extract_coord('x'), extract_coord('y'), extract_coord('z')])
    
    # Dereference Metadata Strings
    def extract_str(ref_array):
        refs = np.array(ref_array).squeeze()
        out = []
        for ref in refs:
            obj = f[ref] if isinstance(ref, h5py.h5r.Reference) else ref
            out.append("".join([chr(c[0]) for c in np.array(obj)]).strip())
        return np.array(out)

    status_desc    = np.array([s.lower() for s in extract_str(ann['all_channel_status_description'])])
    pat_ids_epoch  = np.array([pid.replace('sub-', '') for pid in extract_str(ann['all_clinical_info_epoch']['participant_id'])])
    outcomes_epoch = np.array([out.upper() for out in extract_str(ann['all_clinical_info_epoch']['outcome'])])

# Align Channel Epoch Rows
num_epochs = DS4100_mat_hctsa.shape[0]
if len(pat_ids_epoch) != num_epochs:
    step = len(pat_ids_epoch) // num_epochs
    pat_ids  = pat_ids_epoch[::step][:num_epochs]
    coords   = coords_epoch[::step][:num_epochs]
    status   = status_desc[::step][:num_epochs]
    outcomes = outcomes_epoch[::step][:num_epochs]
else:
    pat_ids, coords, status, outcomes = pat_ids_epoch, coords_epoch, status_desc, outcomes_epoch

is_resected = np.isin(status, ['resected', 'resect,soz', 'soz'])

# ==============================================================================
# 2. Load & Prepare PowCoh Matrices
# ==============================================================================
print("Loading PowCoh dataset...")
with h5py.File(powcoh_path, 'r') as f_pc:
    ds_p = np.array(f_pc['DS4100_Power'])
    ds_c = np.array(f_pc['DS4100_Coherence_Mean'])
    fr_p = np.array(f_pc['Frauscher_Power'])
    fr_c = np.array(f_pc['Frauscher_Coherence_Mean'])
    
    # Orient to (11380 epochs, 5 bands/pairs)
    ds_p = ds_p if ds_p.shape[0] == num_epochs else ds_p.T
    ds_c = ds_c if ds_c.shape[0] == num_epochs else ds_c.T
    fr_p = fr_p if fr_p.shape[0] == Frauscher_mat_hctsa.shape[0] else fr_p.T
    fr_c = fr_c if fr_c.shape[0] == Frauscher_mat_hctsa.shape[0] else fr_c.T
    
    DS4100_mat_powcoh    = np.hstack([ds_p, ds_c])
    Frauscher_mat_powcoh = np.hstack([fr_p, fr_c])

powcoh_feat_names = np.array([
    'Power Delta', 'Power Theta', 'Power Alpha', 'Power Beta', 'Power Gamma',
    'Coh Delta', 'Coh Theta', 'Coh Alpha', 'Coh Beta', 'Coh Gamma'
])

# ==============================================================================
# 3. Determine Top Feature Indices from MATLAB Classification Results
# ==============================================================================
hctsa_feat_names = []
if os.path.exists(external_names_file):
    try:
        ext_data = sio.loadmat(external_names_file, struct_as_record=False, squeeze_me=True)
        ft_obj = ext_data['Feature_names']
        for col_field in ['Name', 'CodeString', 'name', 'codestring']:
            if hasattr(ft_obj, col_field):
                raw_col = getattr(ft_obj, col_field)
                hctsa_feat_names = [str(item.squeeze()) if isinstance(item, np.ndarray) else str(item) for item in raw_col]
                break
    except Exception:
        pass

hctsa_feat_names = np.array(hctsa_feat_names, dtype=str) if len(hctsa_feat_names) == DS4100_mat_hctsa.shape[1] else None

top_features = {'HCTSA': {}, 'PowCoh': {}}

if os.path.exists(results_mat_file):
    print(f"Loading classification results: {results_mat_file}")
    res_data = sio.loadmat(results_mat_file, struct_as_record=False, squeeze_me=True)
    
    # Extract HCTSA Top Feature
    if hasattr(res_data['cv_results'], 'HCTSA_Only'):
        h_struct = getattr(res_data['cv_results'], 'HCTSA_Only')
        if hasattr(h_struct, 'feature_importance'):
            imps_h = np.array(h_struct.feature_importance).squeeze()
            idx_h = int(np.argmax(imps_h))
            name_h = hctsa_feat_names[idx_h] if (hctsa_feat_names is not None and idx_h < len(hctsa_feat_names)) else "CO_fzcglscf_2_5"
            top_features['HCTSA'] = {'idx': idx_h, 'name': name_h, 'mat_ds': DS4100_mat_hctsa, 'mat_fr': Frauscher_mat_hctsa}

    # Extract PowCoh Top Feature
    if hasattr(res_data['cv_results'], 'PowCoh_Only'):
        p_struct = getattr(res_data['cv_results'], 'PowCoh_Only')
        if hasattr(p_struct, 'feature_importance'):
            imps_p = np.array(p_struct.feature_importance).squeeze()
            idx_p = int(np.argmax(imps_p))
            name_p = powcoh_feat_names[idx_p] if idx_p < len(powcoh_feat_names) else f"PowCoh_Feat_{idx_p}"
            top_features['PowCoh'] = {'idx': idx_p, 'name': name_p, 'mat_ds': DS4100_mat_powcoh, 'mat_fr': Frauscher_mat_powcoh}

# Fallbacks if MAT result file is missing
if 'idx' not in top_features['HCTSA']:
    top_features['HCTSA'] = {'idx': 590, 'name': "CO_fzcglscf_2_5", 'mat_ds': DS4100_mat_hctsa, 'mat_fr': Frauscher_mat_hctsa}

if 'idx' not in top_features['PowCoh']:
    top_features['PowCoh'] = {'idx': 9, 'name': "Coh Gamma", 'mat_ds': DS4100_mat_powcoh, 'mat_fr': Frauscher_mat_powcoh}

print(f"Top HCTSA Feature Selected : '{top_features['HCTSA']['name']}' (Index #{top_features['HCTSA']['idx']})")
print(f"Top PowCoh Feature Selected: '{top_features['PowCoh']['name']}' (Index #{top_features['PowCoh']['idx']})")

# Precompute Normative Z-Scores across cohort for both feature domains
norm_z_dict = {}
unique_pats = np.unique(pat_ids)

for f_type in ['PowCoh', 'HCTSA']:
    f_info = top_features[f_type]
    feat_idx, mat_ds, mat_fr = f_info['idx'], f_info['mat_ds'], f_info['mat_fr']
    
    mu_fr0 = np.nanmean(mat_fr[:, feat_idx])
    sd_fr0 = np.nanstd(mat_fr[:, feat_idx])
    if sd_fr0 == 0 or np.isnan(sd_fr0): sd_fr0 = 1.0
    
    z_scores = np.abs((mat_ds[:, feat_idx] - mu_fr0) / sd_fr0)
    norm_z_dict[f_type] = np.nan_to_num(z_scores, nan=0.0, posinf=10.0, neginf=0.0)

    # Print Group Summary Report to Console
    delta_z_norm_list, outcome_binary_list = [], []
    for pid in unique_pats:
        p_mask = (pat_ids == pid)
        p_res  = is_resected[p_mask]
        p_z    = norm_z_dict[f_type][p_mask]
        
        if np.sum(p_res) > 0 and np.sum(~p_res) > 0:
            z_min, z_max = np.min(p_z), np.max(p_z)
            p_z_norm = (p_z - z_min) / (z_max - z_min) if z_max > z_min else np.zeros_like(p_z)
            dz_norm = np.mean(p_z_norm[p_res]) - np.mean(p_z_norm[~p_res])
            pat_out = outcomes[p_mask][0]
            is_succ = 1 if pat_out in ['S', 'SUCCESS', '1A', '1B', 'ENGEL 1', 'ENGEL I'] else 0
            delta_z_norm_list.append(dz_norm)
            outcome_binary_list.append(is_succ)

    dz_arr, out_arr = np.array(delta_z_norm_list), np.array(outcome_binary_list)
    succ_dz, fail_dz = dz_arr[out_arr == 1], dz_arr[out_arr == 0]
    u_stat, p_val = stats.mannwhitneyu(succ_dz, fail_dz, alternative='greater')
    r_s, c_pval   = stats.spearmanr(dz_arr, out_arr)

    print("\n" + "="*75)
    print(f"  GROUP STATISTICAL SUMMARY [{f_type}: {f_info['name']}]")
    print("="*75)
    print(f"  Success Group (N={len(succ_dz)}) Mean Normalized Delta |Z| : {np.mean(succ_dz):.4f} +/- {np.std(succ_dz):.4f}")
    print(f"  Failure Group (N={len(fail_dz)}) Mean Normalized Delta |Z| : {np.mean(fail_dz):.4f} +/- {np.std(fail_dz):.4f}")
    print(f"  Mann-Whitney U Test (Success > Failure)               : U = {u_stat:.1f}, p = {p_val:.4e}")
    print(f"  Spearman Correlation (Normalized Delta |Z| vs Outcome) : r_s = {r_s:.4f}, p = {c_pval:.4e}")
    print("="*75 + "\n")

# ==============================================================================
# 4. Render Merged Horizontal Glass Brain Figures & Companion Abnormality Plots
# ==============================================================================
pats_to_process = unique_pats  # Process test patients
norm_01 = mcolors.Normalize(vmin=0.0, vmax=1.0)
cmap_rdylgn = plt.get_cmap('RdYlGn_r')

for target_pat in pats_to_process:
    p_mask = (pat_ids == target_pat)
    
    # Deduplicate unique physical electrode coordinates
    _, u_idx = np.unique(coords[p_mask], axis=0, return_index=True)
    p_coords = coords[p_mask][u_idx]
    p_resect = is_resected[p_mask][u_idx]
    
    pat_out   = outcomes[p_mask][0]
    out_label = "SUCCESS (Engel I)" if pat_out in ['S', 'SUCCESS', '1A', '1B'] else f"FAILURE ({pat_out})"
    
    # --------------------------------------------------------------------------
    # FIGURE 1: Merged Glass Brain Plot (PowCoh Left | HCTSA Right)
    # --------------------------------------------------------------------------
    fig = plt.figure(figsize=(20, 5.5), dpi=300)
    outer_gs = gridspec.GridSpec(1, 3, width_ratios=[1.0, 1.0, 0.04], wspace=0.05, top=0.88)
    
    gs_left  = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer_gs[0], width_ratios=[1.0, 0.75], height_ratios=[1, 1], wspace=0.0001, hspace=0.02)
    gs_right = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer_gs[1], width_ratios=[1.0, 0.75], height_ratios=[1, 1], wspace=0.0001, hspace=0.02)
    cax      = fig.add_subplot(outer_gs[2])

    domains = [('PowCoh', gs_left), ('HCTSA', gs_right)]
    patient_z_normalized = {}
    
    for f_type, gs_sub in domains:
        p_z_raw = norm_z_dict[f_type][p_mask][u_idx]
        
        # Intra-patient min-max normalization [0, 1]
        z_min, z_max = np.min(p_z_raw), np.max(p_z_raw)
        p_z = (p_z_raw - z_min) / (z_max - z_min) if z_max > z_min else np.zeros_like(p_z_raw)
        patient_z_normalized[f_type] = p_z
        
        res_mask, spa_mask = p_resect, ~p_resect
        res_coords, res_z = p_coords[res_mask], p_z[res_mask]
        spa_coords, spa_z = p_coords[spa_mask], p_z[spa_mask]
        
        ax_axial = fig.add_subplot(gs_sub[:, 0])
        ax_sag   = fig.add_subplot(gs_sub[0, 1])
        ax_cor   = fig.add_subplot(gs_sub[1, 1])
        
        disp_ax  = plotting.plot_glass_brain(None, display_mode='z', axes=ax_axial, annotate=False, black_bg=False)
        disp_sag = plotting.plot_glass_brain(None, display_mode='x', axes=ax_sag, annotate=False, black_bg=False)
        disp_cor = plotting.plot_glass_brain(None, display_mode='y', axes=ax_cor, annotate=False, black_bg=False)
        
        # Overlay Contact Markers using 'RdYlGn_r'
        for disp in [disp_ax, disp_sag, disp_cor]:
            if len(spa_coords) > 0:
                disp.add_markers(
                    marker_coords=spa_coords,
                    marker_color=spa_z,
                    cmap='RdYlGn_r',
                    marker_size=120,
                    norm=norm_01,
                    edgecolor='white',
                    linewidth=1.2
                )
            if len(res_coords) > 0:
                disp.add_markers(
                    marker_coords=res_coords,
                    marker_color=res_z,
                    cmap='RdYlGn_r',
                    marker_size=120,
                    norm=norm_01,
                    edgecolor='black',
                    linewidth=3.2
                )

        # Removed redundant mean/delta values from the sub-panel titles
        ax_axial.set_title(f"[{f_type}] {top_features[f_type]['name']}", fontsize=11, fontweight='bold', pad=10)

    # Common Colorbar
    sm = plt.cm.ScalarMappable(cmap='RdYlGn_r', norm=norm_01)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label('Normalized Feature Abnormality ($[0, 1]$)', fontsize=11, fontweight='bold')
    cbar.set_ticks([0.0, 0.25, 0.50, 0.75, 1.0])

    fig.suptitle(f"Patient {target_pat} | Outcome: {out_label}", fontsize=14, fontweight='bold', y=1.02)

    save_path = os.path.join(output_fig_dir, f"Patient_{target_pat}_GlassBrain_Combined.png")
    plt.savefig(save_path, bbox_inches='tight', dpi=300)
    plt.show()
    
    print(f"Rendered Merged Glass Brain for Patient {target_pat}: Saved to: {save_path}")

    # --------------------------------------------------------------------------
    # FIGURE 2: Companion Contact-Level Abnormality Scatter Plot (PowCoh & HCTSA)
    # --------------------------------------------------------------------------
    fig_dot = plt.figure(figsize=(10, 5), dpi=300)
    
    # Removed colorbar axis, keeping 2 subplots for PowCoh and HCTSA
    gs_dot  = gridspec.GridSpec(1, 2, wspace=0.35)
    
    ax_pow = fig_dot.add_subplot(gs_dot[0])
    ax_hct = fig_dot.add_subplot(gs_dot[1])

    for ax, f_type in zip([ax_pow, ax_hct], ['PowCoh', 'HCTSA']):
        p_z = patient_z_normalized[f_type]
        res_z = p_z[p_resect]
        spa_z = p_z[~p_resect]
        
        m_res = np.mean(res_z) if len(res_z) > 0 else 0.0
        m_spa = np.mean(spa_z) if len(spa_z) > 0 else 0.0
        
        # Calculate One-Sample Bayes Factor (Resected vs. Spared difference)
        diff_vals = []
        if len(res_z) > 0 and len(spa_z) > 0:
            diff_vals = np.array([r - s for r in res_z for s in spa_z])
            bf10 = compute_bayes_factor_one_sample(diff_vals)
        else:
            bf10 = 1.0
            
        bf_str = f"BF = {bf10:.2f}" if bf10 < 1000 else f"BF = {bf10:.2e}"
        
        # Background Mean +/- Std Shading
        if len(res_z) > 0:
            std_res = np.std(res_z)
            ax.add_patch(plt.Rectangle((0.75, max(0.0, m_res - std_res)), 0.5, min(1.0, 2 * std_res), color='#e0e0e0', zorder=1))
            ax.plot([0.75, 1.25], [m_res, m_res], color='black', linewidth=2.0, zorder=2)
            
        if len(spa_z) > 0:
            std_spa = np.std(spa_z)
            ax.add_patch(plt.Rectangle((1.75, max(0.0, m_spa - std_spa)), 0.5, min(1.0, 2 * std_spa), color='#e0e0e0', zorder=1))
            ax.plot([1.75, 2.25], [m_spa, m_spa], color='black', linewidth=2.0, zorder=2)

        # Scatter Contact Dots with Matching Brain Colors & Jitter
        np.random.seed(42)
        if len(res_z) > 0:
            x_jitt_res = 1.0 + (np.random.rand(len(res_z)) - 0.5) * 0.25
            colors_res = cmap_rdylgn(norm_01(res_z))
            ax.scatter(x_jitt_res, res_z, c=colors_res, s=80, edgecolors='black', linewidths=1.2, zorder=3)
            
        if len(spa_z) > 0:
            x_jitt_spa = 2.0 + (np.random.rand(len(spa_z)) - 0.5) * 0.25
            colors_spa = cmap_rdylgn(norm_01(spa_z))
            ax.scatter(x_jitt_spa, spa_z, c=colors_spa, s=80, edgecolors='white', linewidths=0.8, zorder=3)

        # Bracket Annotation for Bayes Factor (BF=...)
        y_max = max(np.max(p_z) if len(p_z) > 0 else 1.0, 0.9)
        bracket_y = min(1.05, y_max + 0.08)
        ax.plot([1.0, 1.0, 2.0, 2.0], [bracket_y - 0.02, bracket_y, bracket_y, bracket_y - 0.02], color='#888888', linewidth=1.5, zorder=4)
        ax.text(1.5, bracket_y + 0.03, bf_str, ha='center', va='bottom', fontsize=11, fontweight='bold', color='#444444')

        ax.set_xticks([1.0, 2.0])
        ax.set_xticklabels(['resected', 'spared'], fontsize=11, fontweight='bold')
        
        # Updated Y-Axis Label to match first figure
        ax.set_ylabel('normalized feature abnormality ([0, 1])', fontsize=10, fontweight='bold')
        ax.set_title(f"[{f_type}] {top_features[f_type]['name']}", fontsize=11, fontweight='bold')
        ax.set_ylim([-0.05, 1.2])
        ax.set_xlim([0.5, 2.5])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig_dot.suptitle(f"Patient {target_pat} | Abnormality Distribution ({out_label})", fontsize=13, fontweight='bold', y=0.98)
    
    save_path_dot = os.path.join(output_fig_dir, f"Patient_{target_pat}_Abnormality_Distribution.png")
    plt.savefig(save_path_dot, bbox_inches='tight', dpi=300)
    plt.show()
    
    print(f"Rendered Abnormality Distribution Plot for Patient {target_pat}: Saved to: {save_path_dot}")

# ==============================================================================
# 5. Render Cohort-Level Abnormality Distribution Plots (Success & Failure)
# ==============================================================================
print("\nGenerating Cohort-Level Summary Abnormality Scatter Plots...")

# Pool intra-patient normalized Z-scores across all subjects by outcome group
cohort_data = {
    'SUCCESS': {'PowCoh': {'res': [], 'spa': []}, 'HCTSA': {'res': [], 'spa': []}},
    'FAILURE': {'PowCoh': {'res': [], 'spa': []}, 'HCTSA': {'res': [], 'spa': []}}
}

for pid in unique_pats:
    p_mask = (pat_ids == pid)
    _, u_idx = np.unique(coords[p_mask], axis=0, return_index=True)
    p_resect = is_resected[p_mask][u_idx]
    
    pat_out = outcomes[p_mask][0]
    out_grp = 'SUCCESS' if pat_out in ['S', 'SUCCESS', '1A', '1B', 'ENGEL 1', 'ENGEL I'] else 'FAILURE'
    
    for f_type in ['PowCoh', 'HCTSA']:
        p_z_raw = norm_z_dict[f_type][p_mask][u_idx]
        z_min, z_max = np.min(p_z_raw), np.max(p_z_raw)
        p_z = (p_z_raw - z_min) / (z_max - z_min) if z_max > z_min else np.zeros_like(p_z_raw)
        
        cohort_data[out_grp][f_type]['res'].extend(p_z[p_resect])
        cohort_data[out_grp][f_type]['spa'].extend(p_z[~p_resect])

# Render cohort-wide figures for Success and Failure
for out_grp in ['SUCCESS', 'FAILURE']:
    fig_cohort = plt.figure(figsize=(10, 5), dpi=300)
    gs_cohort  = gridspec.GridSpec(1, 2, wspace=0.35)
    
    ax_pow = fig_cohort.add_subplot(gs_cohort[0])
    ax_hct = fig_cohort.add_subplot(gs_cohort[1])

    for ax, f_type in zip([ax_pow, ax_hct], ['PowCoh', 'HCTSA']):
        res_z = np.array(cohort_data[out_grp][f_type]['res'])
        spa_z = np.array(cohort_data[out_grp][f_type]['spa'])
        
        m_res = np.mean(res_z) if len(res_z) > 0 else 0.0
        m_spa = np.mean(spa_z) if len(spa_z) > 0 else 0.0
        
        # Calculate Bayes Factor testing Resected vs. Spared across pooled cohort contacts
        diff_vals = []
        if len(res_z) > 0 and len(spa_z) > 0:
            diff_vals = np.array([r - s for r in res_z for s in spa_z])
            bf10 = compute_bayes_factor_one_sample(diff_vals)
        else:
            bf10 = 1.0
            
        bf_str = f"BF = {bf10:.2f}" if bf10 < 1000 else f"BF = {bf10:.2e}"
        
        # Background Mean +/- Std Shading
        if len(res_z) > 0:
            std_res = np.std(res_z)
            ax.add_patch(plt.Rectangle((0.75, max(0.0, m_res - std_res)), 0.5, min(1.0, 2 * std_res), color='#e0e0e0', zorder=1))
            ax.plot([0.75, 1.25], [m_res, m_res], color='black', linewidth=2.0, zorder=2)
            
        if len(spa_z) > 0:
            std_spa = np.std(spa_z)
            ax.add_patch(plt.Rectangle((1.75, max(0.0, m_spa - std_spa)), 0.5, min(1.0, 2 * std_spa), color='#e0e0e0', zorder=1))
            ax.plot([1.75, 2.25], [m_spa, m_spa], color='black', linewidth=2.0, zorder=2)

        # Scatter Contact Dots with Matching Brain Colors & Jitter
        np.random.seed(42)
        if len(res_z) > 0:
            x_jitt_res = 1.0 + (np.random.rand(len(res_z)) - 0.5) * 0.25
            colors_res = cmap_rdylgn(norm_01(res_z))
            ax.scatter(x_jitt_res, res_z, c=colors_res, s=80, edgecolors='black', linewidths=1.2, zorder=3)
            
        if len(spa_z) > 0:
            x_jitt_spa = 2.0 + (np.random.rand(len(spa_z)) - 0.5) * 0.25
            colors_spa = cmap_rdylgn(norm_01(spa_z))
            ax.scatter(x_jitt_spa, spa_z, c=colors_spa, s=80, edgecolors='white', linewidths=0.8, zorder=3)

        # Bracket Annotation for Bayes Factor (BF=...)
        all_pooled = np.concatenate([res_z, spa_z]) if (len(res_z) > 0 and len(spa_z) > 0) else np.array([1.0])
        y_max = max(np.max(all_pooled) if len(all_pooled) > 0 else 1.0, 0.9)
        bracket_y = min(1.05, y_max + 0.08)
        ax.plot([1.0, 1.0, 2.0, 2.0], [bracket_y - 0.02, bracket_y, bracket_y, bracket_y - 0.02], color='#888888', linewidth=1.5, zorder=4)
        ax.text(1.5, bracket_y + 0.03, bf_str, ha='center', va='bottom', fontsize=11, fontweight='bold', color='#444444')

        ax.set_xticks([1.0, 2.0])
        ax.set_xticklabels(['resected', 'spared'], fontsize=11, fontweight='bold')
        ax.set_ylabel('normalized feature abnormality ([0, 1])', fontsize=10, fontweight='bold')
        ax.set_title(f"[{f_type}] {top_features[f_type]['name']}", fontsize=11, fontweight='bold')
        ax.set_ylim([-0.05, 1.2])
        ax.set_xlim([0.5, 2.5])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    grp_title = "SUCCESS (Engel I)" if out_grp == 'SUCCESS' else "FAILURE (Engel II-IV)"
    fig_cohort.suptitle(f"Cohort Group Level | Abnormality Distribution ({grp_title})", fontsize=13, fontweight='bold', y=0.98)
    
    save_path_cohort = os.path.join(output_fig_dir, f"Cohort_AllPatients_Abnormality_Distribution_{out_grp}.png")
    plt.savefig(save_path_cohort, bbox_inches='tight', dpi=300)
    plt.show()
    
    print(f"Rendered Cohort Abnormality Distribution Plot ({out_grp}): Saved to: {save_path_cohort}")