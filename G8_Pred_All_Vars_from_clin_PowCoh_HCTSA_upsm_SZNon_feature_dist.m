clc; close all; clearvars;

% =========================================================================
% 1. PATHS & CONFIGURATION
% =========================================================================
base_folder = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';
hctsa_file   = fullfile(base_folder, 'Data_AllContacts_BothDatasets_epochs.mat');
pow_coh_file = fullfile(base_folder, 'FeatureExtraction/Data_AllContacts_BothDatasets_epochs_PowCoh_per_part.mat');
results_file = fullfile(base_folder, 'Classification_Results_LOSO/outcome/Classification_outcome_LOSO_Results.mat');
external_names_file = '/home/uqhkarim/Normative_HCTSA/Scripts/LittsScripts/Feature_names.mat';

% Output folder for figure exports
output_dir = fullfile(base_folder, 'Figures');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Colors: Successful (Green), Failure (Red)
c_succ = [0.13, 0.65, 0.30]; % Green
c_fail = [0.85, 0.20, 0.20]; % Red

freq_bands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};

% =========================================================================
% 2. LOAD DATA & RECONSTRUCT PATIENT-LEVEL FEATURES
% =========================================================================
fprintf('Loading datasets...\n');
load(hctsa_file, 'DS4100', 'DS4100_annotations'); 
pc_data = load(pow_coh_file); 
results_data = load(results_file);
cv_results = results_data.cv_results;

% Extract HCTSA lookup names if available
hctsa_names_table = {};
if exist(external_names_file, 'file')
    ext_data = load(external_names_file);
    if isfield(ext_data, 'Feature_names')
        FT = ext_data.Feature_names;
        col_to_use = ternary(ismember('Name', FT.Properties.VariableNames), 'Name', 'CodeString');
        hctsa_names_table = cellfun(@(x) char(string(x)), FT.(col_to_use), 'UniformOutput', false);
    end
end

% Extract channel status & patient alignment
status_hctsa = lower(string(DS4100_annotations.all_channel_status_description));
status_pc    = lower(string(pc_data.DS4100_ann.all_channel_status_description));

clin_hctsa = DS4100_annotations.all_clinical_info_epoch;
patient_ids_h = erase(string({clin_hctsa.participant_id}), 'sub-')';
patient_ids_p = erase(string({pc_data.DS4100_ann.all_clinical_info_epoch.participant_id}), 'sub-')';

[~, first_idx] = unique(patient_ids_h, 'stable');
T_clin_master = struct2table(clin_hctsa(first_idx));
T_clin_master.participant_id = erase(string(T_clin_master.participant_id), 'sub-');
T_clin_master = sortrows(T_clin_master, 'participant_id'); 

T_clin_master.outcome   = strtrim(upper(string(T_clin_master.outcome)));
T_clin_master.age       = double(T_clin_master.age);
T_clin_master.age_onset = double(T_clin_master.age_onset);
T_clin_master.duration  = T_clin_master.age - T_clin_master.age_onset;

% Compute patient-level aggregated features (SOZ minus NonSOZ)
unique_pats = T_clin_master.participant_id;
X_hctsa_raw_full = DS4100;
X_pc_raw_full    = [pc_data.DS4100_Power, pc_data.DS4100_Coherence_Mean];

X_h_avg_master = nan(height(T_clin_master), size(X_hctsa_raw_full, 2));
X_p_avg_master = nan(height(T_clin_master), size(X_pc_raw_full, 2));

for p = 1:height(T_clin_master)
    curr_pat = unique_pats(p);
    
    % HCTSA
    mask_h = (patient_ids_h == curr_pat);
    is_soz_h = (status_hctsa == "soz" | status_hctsa == "resect,soz");
    Xs_h = X_hctsa_raw_full(mask_h & is_soz_h', :);
    Xn_h = X_hctsa_raw_full(mask_h & ~is_soz_h', :);
    if ~isempty(Xs_h) && ~isempty(Xn_h)
        X_h_avg_master(p, :) = nanmean(Xs_h, 1) - nanmean(Xn_h, 1);
    end
    
    % PowCoh
    mask_p = (patient_ids_p == curr_pat);
    is_soz_p = (status_pc == "soz" | status_pc == "resect,soz");
    Xs_p = X_pc_raw_full(mask_p & is_soz_p', :);
    Xn_p = X_pc_raw_full(mask_p & ~is_soz_p', :);
    if ~isempty(Xs_p) && ~isempty(Xn_p)
        X_p_avg_master(p, :) = nanmean(Xs_p, 1) - nanmean(Xn_p, 1);
    end
end

X_h_avg_master(isnan(X_h_avg_master)) = 0;
X_p_avg_master(isnan(X_p_avg_master)) = 0;

% Filter valid outcome subjects (S vs F)
valid_mask = ismember(T_clin_master.outcome, {'S', 'F'});
T_clin_sub = T_clin_master(valid_mask, :);
X_h_sub    = X_h_avg_master(valid_mask, :);
X_p_sub    = X_p_avg_master(valid_mask, :);

idx_S = (T_clin_sub.outcome == "S");
idx_F = (T_clin_sub.outcome == "F");

% PowCoh Feature Names Mapping
powcoh_names = [arrayfun(@(b) ['Power ' b{1}], freq_bands, 'Un', 0), ...
                arrayfun(@(b) ['Coh ' b{1}], freq_bands, 'Un', 0)];

% =========================================================================
% GROUP 1: CLINICAL VARIABLES (Unnormalized Counts)
% =========================================================================
clin_preds = {'lesion_status', 'sex', 'target', 'therapy', 'hand', 'implant', 'age', 'age_onset', 'duration'};
cont_vars  = {'age', 'age_onset', 'duration'};

fig1 = figure('Name', 'Clinical Variables Distribution (S vs F)', 'Color', 'w', 'Position', [50 50 1400 850]);
sgtitle('Clinical Predictors Distribution Across Outcome Groups', 'FontSize', 14, 'FontWeight', 'bold');

for i = 1:length(clin_preds)
    vn = clin_preds{i};
    subplot(3, 3, i); hold on;
    
    if ismember(vn, cont_vars)
        vals_S = T_clin_sub.(vn)(idx_S);
        vals_F = T_clin_sub.(vn)(idx_F);
        
        vals_S_clean = vals_S(~isnan(vals_S));
        vals_F_clean = vals_F(~isnan(vals_F));
        
        histogram(vals_S_clean, 10, 'FaceColor', c_succ, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Success (S)');
        histogram(vals_F_clean, 10, 'FaceColor', c_fail, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Failure (F)');
        
        mS = nanmean(vals_S_clean); medS = nanmedian(vals_S_clean);
        mF = nanmean(vals_F_clean); medF = nanmedian(vals_F_clean);
        
        xline(mS, '--', 'Color', c_succ, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        xline(medS, ':', 'Color', c_succ, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        xline(mF, '--', 'Color', c_fail, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        xline(medF, ':', 'Color', c_fail, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        
        bf10 = compute_bayes_factor_two_sample(vals_S_clean, vals_F_clean);
        bf_str = format_bf_string(bf10);
        
        xlabel(strrep(vn, '_', ' ')); ylabel('Patient Count');
        title(sprintf('%s (%s)', strrep(vn, '_', ' '), bf_str), 'FontSize', 10, 'FontWeight', 'bold');
    else
        raw_vals = string(T_clin_sub.(vn));
        cats = unique(raw_vals(~ismissing(raw_vals) & raw_vals ~= "" & raw_vals ~= "nan"));
        
        counts_S = arrayfun(@(c) sum(raw_vals(idx_S) == c), cats);
        counts_F = arrayfun(@(c) sum(raw_vals(idx_F) == c), cats);
        
        total_counts = counts_S + counts_F;
        valid_cats = total_counts >= 3;
        
        cats     = cats(valid_cats);
        counts_S = counts_S(valid_cats);
        counts_F = counts_F(valid_cats);
        
        if ~isempty(cats)
            b = bar(1:length(cats), [counts_S, counts_F], 'grouped');
            b(1).FaceColor = c_succ; b(1).DisplayName = 'Success (S)';
            b(2).FaceColor = c_fail; b(2).DisplayName = 'Failure (F)';
            xticks(1:length(cats)); xticklabels(strrep(cats, '_', ' '));
            xtickangle(25);
        end
        
        ylabel('Patient Count');
        title(strrep(vn, '_', ' '), 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    grid on; set(gca, 'FontSize', 8.5);
    if i == 1, legend('Location', 'best'); end
end

% Save Figure 1
exportgraphics(fig1, fullfile(output_dir, 'Clinical_Variables_Distribution.png'), 'Resolution', 300);

% =========================================================================
% GROUP 2: TOP 10 HCTSA FEATURES (Peak-Normalized to 1.0)
% =========================================================================
imp_hctsa = cv_results.HCTSA_Only.feature_importance;
[~, sorted_idx_h] = sort(imp_hctsa, 'descend');
top10_hctsa_idx = sorted_idx_h(1:10);

fig2 = figure('Name', 'Top 10 HCTSA Features Distribution (S vs F)', 'Color', 'w', 'Position', [100 100 1500 800]);
sgtitle('Top 10 HCTSA Features Distribution Across Outcome Groups (Peak Normalized)', 'FontSize', 14, 'FontWeight', 'bold');

for i = 1:10
    f_idx = top10_hctsa_idx(i);
    
    feat_label = sprintf('HCTSA Feat %d', f_idx);
    if f_idx <= length(hctsa_names_table) && ~isempty(hctsa_names_table{f_idx})
        feat_label = hctsa_names_table{f_idx};
    end
    if strlength(feat_label) > 28
        feat_label = extractBefore(feat_label, 26) + "..";
    end
    
    vals_S = X_h_sub(idx_S, f_idx);
    vals_F = X_h_sub(idx_F, f_idx);
    
    subplot(2, 5, i); hold on;
    
    vals_S_clean = vals_S(~isnan(vals_S));
    vals_F_clean = vals_F(~isnan(vals_F));
    
    combined_vals = [vals_S_clean; vals_F_clean];
    [~, bin_edges] = histcounts(combined_vals, 10);
    
    [counts_S, ~] = histcounts(vals_S_clean, bin_edges);
    if max(counts_S) > 0, counts_S = counts_S / max(counts_S); end
    histogram('BinEdges', bin_edges, 'BinCounts', counts_S, 'FaceColor', c_succ, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Success (S)');
    
    [counts_F, ~] = histcounts(vals_F_clean, bin_edges);
    if max(counts_F) > 0, counts_F = counts_F / max(counts_F); end
    histogram('BinEdges', bin_edges, 'BinCounts', counts_F, 'FaceColor', c_fail, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Failure (F)');
    
    mS = nanmean(vals_S); medS = nanmedian(vals_S);
    mF = nanmean(vals_F); medF = nanmedian(vals_F);
    
    xline(mS, '--', 'Color', c_succ, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(medS, ':', 'Color', c_succ, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(mF, '--', 'Color', c_fail, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(medF, ':', 'Color', c_fail, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    bf10 = compute_bayes_factor_two_sample(vals_S_clean, vals_F_clean);
    bf_str = format_bf_string(bf10);
    
    xlabel('Feature Value'); ylabel('Relative Height (Max = 1)');
    title(sprintf('#%d: %s\n(%s)', i, strrep(feat_label, '_', ' '), bf_str), 'FontSize', 8.5, 'FontWeight', 'bold');

    ylim([0 1.1]); grid on; set(gca, 'FontSize', 8);
    if i == 1, legend('Location', 'best'); end
end

% Save Figure 2
exportgraphics(fig2, fullfile(output_dir, 'Top10_HCTSA_Features_Distribution.png'), 'Resolution', 300);

% =========================================================================
% GROUP 3: TOP 10 POWCOH FEATURES (Peak-Normalized to 1.0)
% =========================================================================
imp_powcoh = cv_results.PowCoh_Only.feature_importance;
[~, sorted_idx_p] = sort(imp_powcoh, 'descend');
top10_powcoh_idx = sorted_idx_p(1:min(10, length(imp_powcoh)));

fig3 = figure('Name', 'Top 10 PowCoh Features Distribution (S vs F)', 'Color', 'w', 'Position', [150 150 1500 800]);
sgtitle('Top 10 PowCoh Features Distribution Across Outcome Groups (Peak Normalized)', 'FontSize', 14, 'FontWeight', 'bold');

for i = 1:length(top10_powcoh_idx)
    f_idx = top10_powcoh_idx(i);
    feat_label = powcoh_names{f_idx};
    
    vals_S = X_p_sub(idx_S, f_idx);
    vals_F = X_p_sub(idx_F, f_idx);
    
    subplot(2, 5, i); hold on;
    
    vals_S_clean = vals_S(~isnan(vals_S));
    vals_F_clean = vals_F(~isnan(vals_F));
    
    combined_vals = [vals_S_clean; vals_F_clean];
    [~, bin_edges] = histcounts(combined_vals, 10);
    
    [counts_S, ~] = histcounts(vals_S_clean, bin_edges);
    if max(counts_S) > 0, counts_S = counts_S / max(counts_S); end
    histogram('BinEdges', bin_edges, 'BinCounts', counts_S, 'FaceColor', c_succ, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Success (S)');
    
    [counts_F, ~] = histcounts(vals_F_clean, bin_edges);
    if max(counts_F) > 0, counts_F = counts_F / max(counts_F); end
    histogram('BinEdges', bin_edges, 'BinCounts', counts_F, 'FaceColor', c_fail, 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Failure (F)');
    
    mS = nanmean(vals_S); medS = nanmedian(vals_S);
    mF = nanmean(vals_F); medF = nanmedian(vals_F);
    
    xline(mS, '--', 'Color', c_succ, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(medS, ':', 'Color', c_succ, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(mF, '--', 'Color', c_fail, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(medF, ':', 'Color', c_fail, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    bf10 = compute_bayes_factor_two_sample(vals_S_clean, vals_F_clean);
    bf_str = format_bf_string(bf10);
    
    xlabel('Feature Value'); ylabel('Relative Height (Max = 1)');
    title(sprintf('#%d: %s\n(%s)', i, feat_label, bf_str), 'FontSize', 9, 'FontWeight', 'bold');
    ylim([0 1.1]); grid on; set(gca, 'FontSize', 8);
    if i == 1, legend('Location', 'best'); end
end

% Save Figure 3
exportgraphics(fig3, fullfile(output_dir, 'Top10_PowCoh_Features_Distribution.png'), 'Resolution', 300);

fprintf('All plots successfully exported to: %s\n', output_dir);

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================
function val = ternary(cond, true_val, false_val)
    if cond, val = true_val; else, val = false_val; end
end

function str = format_bf_string(bf10)
    if bf10 > 10000 || bf10 < 0.001
        str = sprintf('BF_{10} = %.2e', bf10);
    else
        str = sprintf('BF_{10} = %.3f', bf10);
    end
end

function bf10 = compute_bayes_factor_two_sample(x1, x2)
    % Computes two-sample Bayes Factor (BF10) using JZS prior framework
    x1 = x1(isfinite(x1));
    x2 = x2(isfinite(x2));
    n1 = length(x1);
    n2 = length(x2);
    
    if n1 < 2 || n2 < 2
        bf10 = 1;
        return;
    end
    
    s1_sq = var(x1);
    s2_sq = var(x2);
    sp_sq = ((n1 - 1) * s1_sq + (n2 - 1) * s2_sq) / (n1 + n2 - 2);
    
    if sp_sq == 0
        bf10 = 1;
        return;
    end
    
    se = sqrt(sp_sq * (1/n1 + 1/n2));
    t_stat = (mean(x1) - mean(x2)) / se;
    df = n1 + n2 - 2;
    r = sqrt(2) / 2; % Medium scale parameter
    
    neff = (n1 * n2) / (n1 + n2);
    integral_fun = @(g) (1 + neff .* g .* (r^2)).^(-0.5) .* ...
                        (1 + (t_stat.^2) ./ ((1 + neff .* g .* (r^2)) .* df)).^(-(df + 1) / 2) .* ...
                        (g.^(-1.5)) .* exp(-1 ./ (2 .* g));
    
    try
        m_null = (1 + (t_stat.^2) / df).^(-(df + 1) / 2);
        m_alt  = integral_fun_eval(integral_fun);
        bf10   = m_alt / m_null;
        if isnan(bf10) || isinf(bf10), bf10 = 1; end
    catch
        bf10 = 1;
    end
end

function val = integral_fun_eval(fun)
    try
        val = integral(fun, 0, Inf);
    catch
        val = integral(fun, 0, 1000);
    end
end