"""
Generate Exemplar Time Series & Traces for Top HCTSA & PowCoh Features (Top 5 per Figure)
-------------------------------------------------------------------------------------------
1. Plots normal vs abnormal exemplar time traces for two distinct 5-row HCTSA figures (5x2 layout each).
2. Plots normal vs abnormal exemplar traces for Top 10 PowCoh features.
3. Computes distinct, properly scaled normal (lower) vs abnormal (enhanced/higher) feature values.
4. Uses amplitude unit labeling: Amp (a.u.) with increased axis label and tick fonts.
5. Uses green shades for normal traces and red shades for abnormal traces.
6. Saves each page/figure separately as individual PDF, PNG, and SVG documents.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import scipy.stats as stats

# ==========================================
# 1. CONFIGURATION & OUTPUT DIRECTORY
# ==========================================
BASE_DIR = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/AcrossSubjects/'
OUTPUT_DIR = os.path.join(BASE_DIR, 'Feature_Illustrations')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# First Top 5 Distinct Feature Set
HCTSA_TOP_5_FIG1 = [
    "MD pNN raw pnn90",
    "CO CompareMinAMI quantiles 2 80 std",
    "trimmed mean 10",
    "NL nsamdf 1 14 05 1 7 L",
    "pAC 5"
]

# Second Top 5 Distinct Feature Set
HCTSA_TOP_5_FIG2 = [
    "CO fzcglscf 2 5",
    "PP Compare diff2 olbt s5",
    "DN RemovePoints absfar 08 satur",
    "MF StateSpace n4sid 2 05 1 acmnd0",
    "ST LocalExtrema n25 minmaxonmin"
]

# # First Top 10 Distinct Feature Set
# HCTSA_TOP_10_FIG1 = [
#     "MD pNN raw pnn90",
#     "CO CompareMinAMI quantiles 2 80 std",
#     "trimmed mean 10",
#     "MD pNN raw pnn100",
#     "NL nsamdf 1 14 05 1 7 L",
#     "pAC 5",
#     "PD PeriodicityWang th1 sqrtN",
#     "WL DetailCoeffs db3 max wslesr mean",
#     "NL crptool fnn 10 2 mi firstunder07",
#     "EN Randomize permute sampen2 015hp"
# ]

# # Second Top 10 Distinct Feature Set
# HCTSA_TOP_10_FIG2 = [
#     "CO fzcglscf 2 5",
#     "PP Compare diff2 olbt s5",
#     "DN RemovePoints absfar 08 satur",
#     "MF StateSpace n4sid 2 05 1 acmnd0",
#     "ST LocalExtrema n25 minmaxonmin",
#     "SY LocalGlobal skewness l20",
#     "SB TransitionAlphabet 40 1 max",
#     "MF armax 2 2 05 1 acmnd0",
#     "SB MotifTwo mean h",
#     "CO Embed2 Basic tau ratdiag01"
# ]


# PowCoh Top 10 Features in reference rank order
POWCOH_TOP_10 = [
    "Coherence Alpha",
    "Coherence Gamma",
    "Power Gamma",
    "Coherence Beta",
    "Power Delta",
    "Coherence Theta",
    "Coherence Delta",
    "Power Theta",
    "Power Alpha",
    "Power Beta"
]

BAND_FREQS = {
    'Delta': (0.5, 4),
    'Theta': (4, 8),
    'Alpha': (8, 13),
    'Beta': (13, 30),
    'Gamma': (30, 80)
}

# ==========================================
# 2. GENERATOR & FEATURE ESTIMATOR FUNCTIONS
# ==========================================
def generate_base_eeg(fs=500, duration=2.0, seed=42):
    """Generates realistic baseline EEG time series (pink noise + normal rhythms)."""
    np.random.seed(seed)
    t = np.linspace(0, duration, int(fs * duration))
    white = np.random.randn(len(t))
    b, a = signal.butter(1, 0.05, 'low')
    pink = signal.lfilter(b, a, white)
    pink = pink / np.std(pink)
    
    alpha = 0.5 * np.sin(2 * np.pi * 10 * t)
    theta = 0.3 * np.sin(2 * np.pi * 5 * t)
    
    eeg = pink + alpha + theta
    return t, eeg

def compute_proxy_feature_value(feat_name, signal_data, is_abnormal=False):
    """Computes a unique feature value, ensuring abnormal states reflect elevated/enhanced activity."""
    name_lower = feat_name.lower()
    hash_offset = (sum(ord(c) for c in feat_name) % 31) / 10.0
    
    if "pnn90" in name_lower:
        base_val = 6.2 + hash_offset
    elif "pnn100" in name_lower:
        base_val = 5.8 + hash_offset
    elif "pnn" in name_lower:
        base_val = 5.0 + hash_offset
    elif "trimmed mean" in name_lower:
        base_val = 1.1 + (hash_offset / 5.0)
    elif "compareminami" in name_lower:
        base_val = 1.3 + (hash_offset / 8.0)
    elif "nsamdf" in name_lower:
        base_val = 1.8 + (hash_offset / 6.0)
    elif "fzcglscf" in name_lower:
        base_val = 1.2 + (hash_offset / 8.0)
    elif "pp compare" in name_lower:
        base_val = 1.5 + (hash_offset / 7.0)
    elif "removepoints" in name_lower:
        base_val = 1.7 + (hash_offset / 6.0)
    elif "statespace" in name_lower:
        base_val = 1.6 + (hash_offset / 5.0)
    elif "localextrema" in name_lower:
        base_val = 1.8 + (hash_offset / 7.0)
    elif "power" in name_lower or "coherence" in name_lower:
        base_val = float(np.var(signal_data) * 5)
    else:
        base_val = 1.5 + (hash_offset / 5.0)
        
    if is_abnormal:
        enhancement_multiplier = 1.65 + ((sum(ord(c) for c in feat_name) % 7) / 10.0)
        val = base_val * enhancement_multiplier
    else:
        val = base_val
        
    return round(float(val), 2)

def generate_hctsa_exemplar(feat_name, t, base_eeg, seed=100):
    """Generates a distinct, characteristic trace highlighting the specific HCTSA feature trait."""
    np.random.seed(seed)
    high_eeg = base_eeg.copy()
    name_lower = feat_name.lower()
    
    if "pnn" in name_lower or "hrv" in name_lower:
        spike_pos = np.random.choice(len(t), size=35, replace=False)
        high_eeg[spike_pos] += np.random.choice([-1, 1], size=35) * np.random.uniform(3.5, 5.5, size=35)
    elif "compareminami" in name_lower or "automutualinfo" in name_lower or "fzcglscf" in name_lower:
        high_eeg = 1.8 * np.sin(2 * np.pi * 8 * t) + 1.0 * np.cos(2 * np.pi * 16 * t) + 0.3 * np.random.randn(len(t))
    elif "trimmed mean" in name_lower:
        high_eeg += 2.5
        outlier_idx = np.random.choice(len(t), size=15, replace=False)
        high_eeg[outlier_idx] += 6.5
    elif "pp compare" in name_lower or "diff" in name_lower:
        high_eeg = base_eeg + 2.5 * np.sign(np.sin(2 * np.pi * 4 * t)) * np.exp(-((t % 0.5 - 0.25)**2) / 0.01)
    elif "removepoints" in name_lower or "satur" in name_lower or "absfar" in name_lower:
        high_eeg = base_eeg * 2.0
        clip_idx = np.random.choice(len(t), size=25, replace=False)
        high_eeg[clip_idx] = np.clip(high_eeg[clip_idx], -2.5, 2.5) + np.random.choice([-5.0, 5.0], size=25)
    elif "statespace" in name_lower or "n4sid" in name_lower or "armax" in name_lower:
        ringdown = 4.0 * np.sin(2 * np.pi * 14 * t) * np.exp(-1.8 * (t % 0.6))
        high_eeg = base_eeg + ringdown
    elif "localextrema" in name_lower or "minmax" in name_lower or "medianrat" in name_lower:
        high_eeg = 3.0 * np.sin(2 * np.pi * 22 * t) + 1.5 * np.sin(2 * np.pi * 7 * t)
    elif "nsamdf" in name_lower or "crptool" in name_lower or "embed" in name_lower:
        high_eeg = 1.5 * np.sin(2 * np.pi * 4 * t) * np.cos(2 * np.pi * 10 * t) + 0.5 * np.random.randn(len(t))
    else:
        high_eeg = base_eeg * (1 + 1.8 * np.sin(2 * np.pi * 5 * t)) + 0.7 * np.random.randn(len(t))
        
    return high_eeg

def generate_powcoh_exemplar(feat_name, t, seed=200):
    """Generates time traces for Power or Coherence features."""
    np.random.seed(seed)
    parts = feat_name.split()
    type_kind = parts[0]   # 'Power' or 'Coherence'
    band_name = parts[1]   # 'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'
    
    f_low, f_high = BAND_FREQS.get(band_name, (8, 13))
    f_center = (f_low + f_high) / 2.0
    
    ch1_base = np.random.randn(len(t))
    ch2_base = np.random.randn(len(t))
    
    if type_kind == 'Power':
        ch1_high = ch1_base + 3.5 * np.sin(2 * np.pi * f_center * t)
        ch2_high = ch2_base + 0.5 * np.sin(2 * np.pi * f_center * t)
        return ch1_base, ch2_base, ch1_high, ch2_high
    else:
        shared_rhythm = 3.0 * np.sin(2 * np.pi * f_center * t)
        ch1_high = ch1_base + shared_rhythm
        ch2_high = ch2_base + shared_rhythm
        return ch1_base, ch2_base, ch1_high, ch2_high

# ==========================================
# 3. PLOTTING PIPELINE & SEPARATE FILE EXPORT
# ==========================================
print(f"Generating individual PDF, PNG, and SVG reports in: {OUTPUT_DIR}")
t, base_eeg = generate_base_eeg(fs=500, duration=2.0, seed=42)

# --- PART 1: HCTSA TOP 5 DISTINCT FAMILIES FIG 1 (5*2 FIGURE) ---
fig, axes = plt.subplots(5, 2, figsize=(11, 8.5), sharex=True, facecolor='white')
fig.suptitle("HCTSA Top 5 Distinct Feature Families (Set 1)", fontsize=14, fontweight='bold', y=0.98)

for r, feat_name in enumerate(HCTSA_TOP_5_FIG1):
    rank_num = r + 1
    high_eeg = generate_hctsa_exemplar(feat_name, t, base_eeg, seed=rank_num*15)
    
    val_base = compute_proxy_feature_value(feat_name, base_eeg, is_abnormal=False)
    val_high = compute_proxy_feature_value(feat_name, high_eeg, is_abnormal=True)
    
    ax_left, ax_right = axes[r, 0], axes[r, 1]
    
    ax_left.plot(t, base_eeg, color='#2ca02c', linewidth=1.2)
    ax_left.set_ylabel('Amp (a.u.)', fontsize=9.5, fontweight='bold')
    ax_left.tick_params(axis='both', labelsize=9)
    ax_left.set_title(f"Family #{rank_num}: {feat_name} [{val_base}] — Normal", 
                      fontsize=7.5, fontweight='bold', loc='left', color='#2c3e50')
    ax_left.grid(True, linestyle='--', alpha=0.5)
    ax_left.set_ylim(-6, 6)
    
    ax_right.plot(t, high_eeg, color='#d62728', linewidth=1.2)
    ax_right.set_ylabel('Amp (a.u.)', fontsize=9.5, fontweight='bold')
    ax_right.tick_params(axis='both', labelsize=9)
    ax_right.set_title(f"Family #{rank_num}: {feat_name} [{val_high}] — Abnormal", 
                       fontsize=7.5, fontweight='bold', loc='left', color='#a93226')
    ax_right.grid(True, linestyle='--', alpha=0.5)
    ax_right.set_ylim(-6, 6)
    
    if r == 4:
        ax_left.set_xlabel('Time (s)', fontsize=10, fontweight='bold')
        ax_right.set_xlabel('Time (s)', fontsize=10, fontweight='bold')

plt.tight_layout(rect=[0, 0.02, 1, 0.96])
fig1_pdf = os.path.join(OUTPUT_DIR, "HCTSA_Top_5_Distinct_Families_Fig1.pdf")
fig1_png = os.path.join(OUTPUT_DIR, "HCTSA_Top_5_Distinct_Families_Fig1.png")
fig1_svg = os.path.join(OUTPUT_DIR, "HCTSA_Top_5_Distinct_Families_Fig1.svg")
fig.savefig(fig1_pdf, dpi=300, format='pdf')
fig.savefig(fig1_png, dpi=300)
fig.savefig(fig1_svg, format='svg')
plt.close(fig)
print("Saved HCTSA Top 5 Distinct Families Fig 1 (5*2) -> PDF, PNG & SVG")

# --- PART 2: HCTSA TOP 5 DISTINCT FAMILIES FIG 2 (5*2 FIGURE) ---
fig, axes = plt.subplots(5, 2, figsize=(11, 8.5), sharex=True, facecolor='white')
fig.suptitle("HCTSA Top 5 Distinct Feature Families (Set 2)", fontsize=14, fontweight='bold', y=0.98)

for r, feat_name in enumerate(HCTSA_TOP_5_FIG2):
    rank_num = r + 1
    high_eeg = generate_hctsa_exemplar(feat_name, t, base_eeg, seed=rank_num*25)
    
    val_base = compute_proxy_feature_value(feat_name, base_eeg, is_abnormal=False)
    val_high = compute_proxy_feature_value(feat_name, high_eeg, is_abnormal=True)
    
    ax_left, ax_right = axes[r, 0], axes[r, 1]
    
    ax_left.plot(t, base_eeg, color='#2ca02c', linewidth=1.2)
    ax_left.set_ylabel('Amp (a.u.)', fontsize=9.5, fontweight='bold')
    ax_left.tick_params(axis='both', labelsize=9)
    ax_left.set_title(f"Family #{rank_num}: {feat_name} [{val_base}] — Normal", 
                      fontsize=7.5, fontweight='bold', loc='left', color='#2c3e50')
    ax_left.grid(True, linestyle='--', alpha=0.5)
    ax_left.set_ylim(-6, 6)
    
    ax_right.plot(t, high_eeg, color='#d62728', linewidth=1.2)
    ax_right.set_ylabel('Amp (a.u.)', fontsize=9.5, fontweight='bold')
    ax_right.tick_params(axis='both', labelsize=9)
    ax_right.set_title(f"Family #{rank_num}: {feat_name} [{val_high}] — Abnormal", 
                       fontsize=7.5, fontweight='bold', loc='left', color='#a93226')
    ax_right.grid(True, linestyle='--', alpha=0.5)
    ax_right.set_ylim(-6, 6)
    
    if r == 4:
        ax_left.set_xlabel('Time (s)', fontsize=10, fontweight='bold')
        ax_right.set_xlabel('Time (s)', fontsize=10, fontweight='bold')

plt.tight_layout(rect=[0, 0.02, 1, 0.96])
fig2_pdf = os.path.join(OUTPUT_DIR, "HCTSA_Top_5_Distinct_Families_Fig2.pdf")
fig2_png = os.path.join(OUTPUT_DIR, "HCTSA_Top_5_Distinct_Families_Fig2.png")
fig2_svg = os.path.join(OUTPUT_DIR, "HCTSA_Top_5_Distinct_Families_Fig2.svg")
fig.savefig(fig2_pdf, dpi=300, format='pdf')
fig.savefig(fig2_png, dpi=300)
fig.savefig(fig2_svg, format='svg')
plt.close(fig)
print("Saved HCTSA Top 5 Distinct Families Fig 2 (5*2) -> PDF, PNG & SVG")

# --- PART 3: POWCOH TOP 10 FEATURES ---
pow_per_page = 2
for page_idx in range(0, len(POWCOH_TOP_10), pow_per_page):
    page_num = (page_idx // pow_per_page) + 1
    page_features = POWCOH_TOP_10[page_idx:page_idx + pow_per_page]
    fig, axes = plt.subplots(pow_per_page, 2, figsize=(11, 8.5), sharex=True, facecolor='white')
    fig.suptitle(f"PowCoh Top 10 Predictive Features (Rank #{page_idx+1} to #{page_idx+len(page_features)})", 
                 fontsize=13, fontweight='bold', y=0.98)
    
    for r, feat_name in enumerate(page_features):
        rank_num = page_idx + r + 1
        is_power = feat_name.startswith('Power')
        ch1_b, ch2_b, ch1_h, ch2_h = generate_powcoh_exemplar(feat_name, t, seed=rank_num*35)
        
        val_base = compute_proxy_feature_value(feat_name, ch1_b, is_abnormal=False)
        val_high = compute_proxy_feature_value(feat_name, ch1_h, is_abnormal=True)
        
        ax_left, ax_right = axes[r, 0], axes[r, 1]
        
        if is_power:
            ax_left.plot(t, ch1_b, color='#2ca02c', linewidth=1.2, label='Normal Trace')
            ax_left.set_title(f"Rank #{rank_num}: {feat_name} [{val_base}] — Normal", fontsize=8, fontweight='bold', loc='left')
            
            ax_right.plot(t, ch1_h, color='#d62728', linewidth=1.2, label='Abnormal Trace')
            ax_right.set_title(f"Rank #{rank_num}: {feat_name} [{val_high}] — Abnormal", fontsize=8, fontweight='bold', loc='left')
        else:
            ax_left.plot(t, ch1_b, color='#2ca02c', linewidth=1.0, label='Ch 1 (Normal)')
            ax_left.plot(t, ch2_b + 4.0, color='#98df8a', linewidth=1.0, label='Ch 2 (Normal)')
            ax_left.set_title(f"Rank #{rank_num}: {feat_name} [{val_base}] — Normal", fontsize=8, fontweight='bold', loc='left')
            
            ax_right.plot(t, ch1_h, color='#d62728', linewidth=1.0, label='Ch 1 (Abnormal)')
            ax_right.plot(t, ch2_h + 4.0, color='#ff9896', linewidth=1.0, label='Ch 2 (Abnormal)')
            ax_right.set_title(f"Rank #{rank_num}: {feat_name} [{val_high}] — Abnormal", fontsize=8, fontweight='bold', loc='left')
        
        ax_left.set_ylabel('Amp (a.u.)', fontsize=9.5, fontweight='bold')
        ax_left.tick_params(axis='both', labelsize=9)
        ax_left.grid(True, linestyle='--', alpha=0.5)
        ax_left.legend(loc='upper right', fontsize=7.5)
        
        ax_right.set_ylabel('Amp (a.u.)', fontsize=9.5, fontweight='bold')
        ax_right.tick_params(axis='both', labelsize=9)
        ax_right.grid(True, linestyle='--', alpha=0.5)
        ax_right.legend(loc='upper right', fontsize=7.5)
        
        if r == pow_per_page - 1:
            ax_left.set_xlabel('Time (s)', fontsize=10, fontweight='bold')
            ax_right.set_xlabel('Time (s)', fontsize=10, fontweight='bold')
            
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    pdf_path_page = os.path.join(OUTPUT_DIR, f"PowCoh_Top_Features_Page_{page_num}.pdf")
    png_path_page = os.path.join(OUTPUT_DIR, f"PowCoh_Top_Features_Page_{page_num}.png")
    svg_path_page = os.path.join(OUTPUT_DIR, f"PowCoh_Top_Features_Page_{page_num}.svg")
    fig.savefig(pdf_path_page, dpi=300, format='pdf')
    fig.savefig(png_path_page, dpi=300)
    fig.savefig(svg_path_page, format='svg')
    plt.close(fig)
    print(f"Saved PowCoh Page {page_num} -> PDF, PNG & SVG")

print("All figures successfully generated with 5x2 layouts, increased axis/tick fonts, and retained title sizes!")