clc;
close all;
clear;

% =========================================================================
% --- CONFIGURATION & PARAMETERS ---
% =========================================================================
N_top_candidates = 100;    % Candidate pool for non-overlapping top features
min_index_distance = 80;   % Minimum x-axis distance to prevent label overlap
PLOT_STANDARDIZED = false; 

% Summary Bar Plot Ns
summary_N.PowCoh = 10;
summary_N.HCTSA  = 40;

base_folder = '/scratch/user/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';
external_names_file = '/scratch/user/uqhkarim/Normative_HCTSA/Scripts/LittsScripts/Feature_names.mat';

% Unified Output Directory
fig_output_dir = fullfile(base_folder, 'AcrossSubjects/Feature_Distributions');
if ~exist(fig_output_dir, 'dir'), mkdir(fig_output_dir); end

feature_modes = {'PowCoh', 'HCTSA'};

% Colors matching Matrix Heatmaps (Orange for PowCoh, Blue for HCTSA)
feature_colors.PowCoh = [0.8500, 0.3250, 0.0980]; 
feature_colors.HCTSA  = [0.0000, 0.4470, 0.7410]; 

% Raw Input Data Files
input_files.PowCoh = fullfile(base_folder, 'FeatureExtraction/DS4100_Power_Coh_Normalized_Final_per_part.mat');
input_files.HCTSA  = fullfile(base_folder, 'DS4100_ZScores_FromNormativeAtlas_Nrmlzd_epochs.mat');

% Classification Output Files
results_files.PowCoh = fullfile(base_folder, 'AcrossSubjects/PowerCoh_Results/PowerCoherence/AcrossSubjects_Final_Results_both_feature_sets_Repeated_fixed_resect.mat');
results_files.HCTSA  = fullfile(base_folder, 'AcrossSubjects/HCTSA_Results/AcrossSubjects_Final_Results_both_feature_sets_Repeated_fixed_resect.mat');

% Histogram Colors (SOZ vs Non-SOZ)
c_soz = [0.8500, 0.1500, 0.1500];
c_non = [0.1500, 0.7000, 0.1500];

% Custom Continuous Colormaps
cmap_hctsa  = [linspace(0, 0.05, 128)', linspace(0, 0.35, 128)', linspace(0, 0.75, 128)'];
cmap_powcoh = [linspace(0, 0.90, 128)', linspace(0, 0.45, 128)', linspace(0, 0.05, 128)'];

combined_data = struct();

% =========================================================================
% --- MAIN EXECUTION PIPELINE ---
% =========================================================================
for m = 1:length(feature_modes)
    mode = feature_modes{m};
    target_file = results_files.(mode);
    
    if ~exist(target_file, 'file')
        fprintf('Skipping %s (Results file not found at %s)\n', mode, target_file);
        continue; 
    end
    
    fprintf('\n=======================================================\n');
    fprintf('PROCESSING FEATURE IMPORTANCES & DISTRIBUTIONS: %s\n', mode);
    fprintf('=======================================================\n');
    
    res_data = load(target_file);
    if ~isfield(res_data, 'all_importances_matrix')
        fprintf('Warning: all_importances_matrix not found in %s file. Skipping...\n', mode);
        continue;
    end
    
    mat = res_data.all_importances_matrix;
    [n_subs, n_feats] = size(mat);
    
    % --- 1. Extract Feature Names ---
    feat_names = cell(n_feats, 1);
    if isfield(res_data, 'feature_names_kept')
        temp_names = cellstr(res_data.feature_names_kept(:));
        feat_names(1:min(n_feats, length(temp_names))) = temp_names(1:min(n_feats, length(temp_names)));
    end
    
    if strcmp(mode, 'HCTSA') && exist(external_names_file, 'file') && ...
       (isempty(feat_names{1}) || startsWith(feat_names{1}, 'Feature_'))
        ext_data = load(external_names_file);
        if isfield(ext_data, 'Feature_names')
            FT = ext_data.Feature_names;
            num_names = height(FT); 
            col_to_use = 'Name';
            if ~ismember('Name', FT.Properties.VariableNames) && ismember('CodeString', FT.Properties.VariableNames)
                col_to_use = 'CodeString';
            end
            for i = 1:min(n_feats, num_names)
                val = FT.(col_to_use){i}; 
                if iscell(val), feat_names{i} = char(val{1});
                else, feat_names{i} = char(val); end
            end
        end
    end
    
    for i = 1:n_feats
        if isempty(feat_names{i})
            feat_names{i} = sprintf('%s_%d', mode, i);
        end
    end

    % --- 2. Matrix Clean-up & Normalization ---
    nan_cols = all(isnan(mat), 1);
    mat_clean = mat(:, ~nan_cols);
    feat_names_clean = feat_names(~nan_cols);
    
    v_min = nanpercentile(mat_clean(:), 1);
    v_max = nanpercentile(mat_clean(:), 98);
    if v_max == v_min, v_max = v_min + 1e-6; end
    mat_norm = (mat_clean - v_min) ./ (v_max - v_min);
    mat_norm(mat_norm > 1) = 1; mat_norm(mat_norm < 0) = 0;
    
    avg_feat = nanmean(mat_clean, 1);
    [~, sorted_idx] = sort(avg_feat, 'descend');
    
    % Select non-overlapping features across rank
    if strcmp(mode, 'PowCoh')
        selected_display_indices = sorted_idx; % Plot all 10 features for PowCoh
    else
        % Candidate selection across rank gradient for HCTSA
        candidate_pool = sorted_idx;
        selected_display_indices = [];
        for c = 1:length(candidate_pool)
            idx = candidate_pool(c);
            if isempty(selected_display_indices) || all(abs(selected_display_indices - idx) >= min_index_distance)
                selected_display_indices(end+1) = idx; %#ok<AGROW>
            end
        end
    end

    % =========================================================================
    % --- SECTION 1: HEATMAP & AVERAGE TRACE PLOT ---
    % =========================================================================
    fig_title = sprintf('%s Feature Importances Across Patient Folds', mode);
    figure('Color', 'w', 'Name', fig_title, 'Position', [50 50 1400 850]);
    
    axMain = subplot('Position', [0.08 0.18 0.85 0.72]);
    imagesc(mat_norm);
    
    if strcmp(mode, 'HCTSA')
        colormap(axMain, cmap_hctsa);
        trace_color = [0.0, 0.4470, 0.7410];
        font_sz = 6.5;
    else
        colormap(axMain, cmap_powcoh);
        trace_color = [0.8500, 0.3250, 0.0980];
        font_sz = 9.5; 
    end
    
    ylabel('Patient #', 'FontSize', 11, 'FontWeight', 'bold');
    set(gca, 'XTick', [], 'TickDir', 'out'); 
    title(fig_title, 'FontSize', 14, 'FontWeight', 'bold');
    
    hold on;
    n_rows = size(mat_norm, 1);
    n_display = length(selected_display_indices);

    % --- Feature Text Placement ---
    for rank_pos = 1:n_display
        idx = selected_display_indices(rank_pos);
        name_str = strrep(feat_names_clean{idx}, '_', ' ');
        
        if n_display > 1
            rel_rank = (rank_pos - 1) / (n_display - 1);
        else
            rel_rank = 0;
        end
        
        if strcmp(mode, 'HCTSA')
            % Truncate long HCTSA names to maximum 25 chars + '...' to keep string height small
            if length(name_str) > 28
                name_str = [name_str(1:25) '...'];
            end
            % Strictly constrain Y center inside [25%, 75%] of row bounds
            y_pos = 0.25 * n_rows + rel_rank * (0.50 * n_rows);
        else
            % PowCoh linear mapping within 15%-85% range
            y_pos = 0.15 * n_rows + rel_rank * (0.70 * n_rows);
        end
        
        txt_color = [1, 1, 1 - rel_rank];
        
        text(idx, y_pos, name_str, 'Rotation', 90, 'FontSize', font_sz, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', txt_color, ...
            'Interpreter', 'none', 'Clipping', 'on');
            
        line([idx idx], [0.5 n_rows+0.5], 'Color', [txt_color 0.15], 'LineWidth', 0.5);
    end
    
    % Bottom Horizontal Average Trace
    axBottom = subplot('Position', [0.08 0.07 0.85 0.08]);
    plot(avg_feat, 'Color', trace_color, 'LineWidth', 1.5);
    xlabel(sprintf('Filtered Feature Index (Total: %d)', size(mat_norm,2)), 'FontSize', 10, 'FontWeight', 'bold'); 
    ylabel('Mean Imp.', 'FontSize', 9, 'FontWeight', 'bold');
    grid on; xlim([0.5, size(mat_norm,2)+0.5]);
    
    cb = colorbar(axMain, 'Position', [0.94 0.18 0.015 0.72]);
    ylabel(cb, 'Normalized Importance Error Delta (a.u.)', 'FontSize', 10, 'FontWeight', 'bold');
    
    linkaxes([axMain, axBottom], 'x');
    
    save_base = fullfile(fig_output_dir, sprintf('Matrix_Feature_Importance_%s', mode));
    exportgraphics(gcf, [save_base '.png'], 'Resolution', 300);
    exportgraphics(gcf, [save_base '.pdf'], 'ContentType', 'vector');
    fprintf('Saved Heatmap (PNG & PDF): %s\n', save_base);
    close(gcf);

    % =========================================================================
    % --- SECTION 2: SUMMARY BOXPLOT WITH DATAPOINTS (TOP N) ---
    % =========================================================================
    curr_n = summary_N.(mode);
    n_sel = min(curr_n, length(sorted_idx));
    
    summary_indices = sorted_idx(1:n_sel);
    mat_sel = mat_clean(:, summary_indices);
    
    x_val_all = [];
    y_cat_all = cell(0, 1);
    
    labels_reversed = feat_names_clean(summary_indices(end:-1:1));
    for k = 1:length(labels_reversed)
        labels_reversed{k} = strrep(labels_reversed{k}, '_', ' ');
    end
    
    for r = 1:n_sel
        feat_col_idx = summary_indices(n_sel - r + 1);
        patient_data = mat_clean(:, feat_col_idx);
        patient_data = patient_data(~isnan(patient_data));
        
        x_val_all = [x_val_all; patient_data(:)]; %#ok<AGROW>
        y_cat_all = [y_cat_all; repmat(labels_reversed(r), length(patient_data), 1)]; %#ok<AGROW>
    end
    
    y_cat_all = categorical(y_cat_all, labels_reversed, 'Ordinal', true);
    
    fig_height = max(550, n_sel * 22);
    fig_box = figure('Color', 'w', 'Name', sprintf('%s_Top_%d_Boxplot', mode, n_sel), ...
                     'Position', [100, 100, 1000, fig_height]);
    
    bc = boxchart(y_cat_all, x_val_all, 'Orientation', 'horizontal', ...
                  'BoxFaceColor', feature_colors.(mode), ...
                  'BoxMedianLineColor', 'k', ...
                  'MarkerColor', [0.2 0.2 0.2], ...
                  'MarkerStyle', 'o', 'MarkerSize', 3.5);
              
    bc.JitterOutliers = 'on';
    bc.MarkerStyle = 'none';
    
    hold on;
    unique_cats = categories(y_cat_all);
    for c = 1:length(unique_cats)
        pts = x_val_all(y_cat_all == unique_cats{c});
        y_jitter = c + (rand(size(pts)) - 0.5) * 0.35;
        plot(pts, y_jitter, '.', 'Color', [0.15 0.15 0.15 0.6], 'MarkerSize', 8);
    end
    
    set(gca, 'FontSize', 8, 'TickLabelInterpreter', 'none');
    xlabel('Importance Delta Error Across Patient Folds', 'FontWeight', 'bold');
    ylabel('Feature Rank', 'FontWeight', 'bold');
    title(sprintf('%s: Top %d Predictive Features (Distribution Across Patients)', mode, n_sel), ...
          'FontSize', 12, 'FontWeight', 'bold');
    grid on; box on;
    
    save_bar_base = fullfile(fig_output_dir, sprintf('Summary_Top_%d_Features_%s', n_sel, mode));
    exportgraphics(fig_box, [save_bar_base '.png'], 'Resolution', 300);
    exportgraphics(fig_box, [save_bar_base '.pdf'], 'ContentType', 'vector');
    fprintf('Saved Boxplot (PNG & PDF): %s\n', save_bar_base);
    close(fig_box);

    % =========================================================================
    % --- SECTION 3: FEATURE DISTRIBUTION HISTOGRAMS ---
    % =========================================================================
    if strcmp(mode, 'PowCoh')
        raw_struct = load(input_files.PowCoh, 'DS4100_Final_Z', 'DS4100_ann', 'feature_names');
        X_raw = abs(raw_struct.DS4100_Final_Z(:, 1:10)); 
        all_channel_status = raw_struct.DS4100_ann.all_channel_status_description;
    else
        raw_struct = load(input_files.HCTSA, 'DS4100_z_scores_from_normative', 'DS4100_annotations');
        X_raw = abs(raw_struct.DS4100_z_scores_from_normative); 
        all_channel_status = raw_struct.DS4100_annotations.all_channel_status_description;
    end
    
    Y_raw = NaN(size(X_raw, 1), 1);
    Y_raw(strcmp(all_channel_status, 'resect,soz')) = 1;
    Y_raw(strcmp(all_channel_status, 'soz') | strcmp(all_channel_status, 'n/a') | strcmp(all_channel_status, 'resect')) = 0;
    
    valid_idx = ~isnan(Y_raw);
    X_full = X_raw(valid_idx, :);
    Y_full = Y_raw(valid_idx);
    
    features_to_keep = true(1, size(X_full, 2));
    if strcmp(mode, 'HCTSA')
        proport = 0.9;
        for j = 1:size(X_full, 2)
            current_col = X_full(:, j);
            if (sum(isnan(current_col) | current_col == 0) > proport * length(current_col)) || all(isinf(current_col))
                features_to_keep(j) = false;
            end
        end
        X_full = X_full(:, features_to_keep);
    end
    
    outls = isoutlier(X_full, 'percentiles', [0.03 99.7]);
    X_full(outls) = NaN;
    for j = 1:size(X_full, 2)
        m_val = nanmean(X_full(:, j));
        if isnan(m_val), m_val = 0; end
        X_full(isnan(X_full(:, j)), j) = m_val;
    end
    
    if PLOT_STANDARDIZED
        mu = mean(X_full, 1); sigma = std(X_full, 0, 1); sigma(sigma == 0) = 1;
        X_plot_data = (X_full - mu) ./ sigma;
        xlabel_str = 'Scaled Classifier Feature Value';
        title_suffix = 'Standardized Input Features';
    else
        X_plot_data = X_full;
        xlabel_str = 'Absolute Normative Z-Score (|Z|)';
        title_suffix = 'Positive Absolute Z-Scores (|Z|)';
    end
    
    mat_imp = res_data.all_importances_matrix;
    if size(mat_imp, 2) == length(features_to_keep)
        mat_imp = mat_imp(:, features_to_keep);
    end
    
    avg_imp = nanmean(mat_imp, 1);
    valid_mask = ~isnan(avg_imp) & ~all(isnan(X_plot_data), 1);
    valid_indices = find(valid_mask);
    
    [~, sort_order] = sort(avg_imp(valid_indices), 'descend');
    sorted_valid_cols = valid_indices(sort_order);
    total_valid = length(sorted_valid_cols);
    
    if strcmp(mode, 'HCTSA')
        top_5_cols   = sorted_valid_cols(1:min(5, total_valid));
        worst_5_cols = sorted_valid_cols(max(1, total_valid-4):total_valid);
        plot_cols   = [top_5_cols, worst_5_cols];
        rank_labels = [1:length(top_5_cols), (total_valid-length(worst_5_cols)+1):total_valid];
        grid_rows   = 2;
        fig_title_str = sprintf('HCTSA: Best & Worst 5 Features - %s', title_suffix);
        
        top_10_cols = sorted_valid_cols(1:min(10, total_valid));
        combined_data.hctsa_best_soz = reshape(X_plot_data(Y_full == 1, top_10_cols), [], 1);
        combined_data.hctsa_best_non = reshape(X_plot_data(Y_full == 0, top_10_cols), [], 1);
    else
        plot_cols   = sorted_valid_cols(1:min(10, total_valid));
        rank_labels = 1:length(plot_cols);
        grid_rows   = 2;
        fig_title_str = sprintf('PowCoh: Ranked Features - %s', title_suffix);
        
        combined_data.powcoh_all_soz = reshape(X_plot_data(Y_full == 1, :), [], 1);
        combined_data.powcoh_all_non = reshape(X_plot_data(Y_full == 0, :), [], 1);
    end

    fig_hist = figure('Name', [mode ' - Classifier Feature Distributions'], 'Color', 'w', 'Position', [50 50 1600 550]);
    sgtitle(fig_title_str, 'FontSize', 14, 'FontWeight', 'bold');

    is_soz = (Y_full == 1);
    is_non = (Y_full == 0);

    for p = 1:length(plot_cols)
        f_col    = plot_cols(p);
        rank_num = rank_labels(p);
        f_name   = char(feat_names_clean{min(f_col, length(feat_names_clean))});
        
        soz_vals = X_plot_data(is_soz, f_col);
        non_vals = X_plot_data(is_non, f_col);
        
        subplot(grid_rows, 5, p); hold on;
        
        if isempty(soz_vals) || isempty(non_vals)
            title(sprintf('Rank #%d: Insufficient Data', rank_num), 'FontSize', 8.5);
            grid on; set(gca, 'FontSize', 8); continue;
        end
        
        median_soz = median(soz_vals);
        median_non = median(non_vals);
        
        all_vals = [soz_vals; non_vals];
        min_v = min(all_vals);
        max_v = prctile(all_vals, 99.0);
        if isnan(min_v) || isnan(max_v) || (min_v >= max_v)
            min_v = min_v - 1e-3; max_v = max_v + 1e-3;
        end
        
        bin_edges = linspace(min_v, max_v, 30);
        
        h_soz = histogram(soz_vals, bin_edges, 'Normalization', 'pdf', ...
            'FaceColor', c_soz, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'SOZ');
        h_non = histogram(non_vals, bin_edges, 'Normalization', 'pdf', ...
            'FaceColor', c_non, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Non-SOZ');
        
        xlim([min_v max_v]);
        max_density = max([h_soz.Values, h_non.Values]);
        if isnan(max_density) || max_density == 0, max_density = 1; end
        ylim([0 max_density * 1.15]);
        
        line([median_soz median_soz], ylim, 'Color', c_soz * 0.8, 'LineWidth', 2.0, 'LineStyle', '--', 'HandleVisibility', 'off');
        line([median_non median_non], ylim, 'Color', c_non * 0.8, 'LineWidth', 2.0, 'LineStyle', '--', 'HandleVisibility', 'off');
        
        clean_title = strrep(f_name, '_', ' ');
        
        title(sprintf('Rank #%d: %s', rank_num, clean_title), 'FontSize', 8.5, 'FontWeight', 'bold');
        xlabel(xlabel_str); ylabel('Probability Density');
        grid on; set(gca, 'FontSize', 8);
        
        if p == 1, legend({'SOZ', 'Non-SOZ'}, 'Location', 'northwest'); end
    end
    
    save_hist_base = fullfile(fig_output_dir, sprintf('%s_Classifier_Input_Feature_Distributions', mode));
    exportgraphics(fig_hist, [save_hist_base '.png'], 'Resolution', 300);
    exportgraphics(fig_hist, [save_hist_base '.pdf'], 'ContentType', 'vector');
    fprintf('Saved Histograms (PNG & PDF): %s\n', save_hist_base);
    close(fig_hist);
end

% =========================================================================
% --- COMBINED COMPARISON FIGURE ---
% =========================================================================
fig_comb = figure('Name', 'Combined Feature Group Comparisons', 'Color', 'w', 'Position', [100 100 1200 500]);
sgtitle('Aggregated Comparison: All PowCoh Features vs. Top 10 HCTSA Features', 'FontSize', 14, 'FontWeight', 'bold');

panel_specs = {
    'Combined All PowCoh Features',     combined_data.powcoh_all_soz,   combined_data.powcoh_all_non,   1;
    'Combined Best 10 HCTSA Features',  combined_data.hctsa_best_soz,  combined_data.hctsa_best_non,   2
};

for p = 1:size(panel_specs, 1)
    p_title  = panel_specs{p, 1};
    soz_vals = panel_specs{p, 2};
    non_vals = panel_specs{p, 3};
    p_idx    = panel_specs{p, 4};
    
    subplot(1, 2, p_idx); hold on;
    
    soz_vals = soz_vals(isfinite(soz_vals));
    non_vals = non_vals(isfinite(non_vals));
    
    median_soz = median(soz_vals);
    median_non = median(non_vals);
    
    all_vals = [soz_vals; non_vals];
    min_v = min(all_vals);
    max_v = prctile(all_vals, 99.0);
    if isnan(min_v) || isnan(max_v) || (min_v >= max_v)
        min_v = min_v - 1e-3; max_v = max_v + 1e-3;
    end
    
    bin_edges = linspace(min_v, max_v, 40);
    
    h_soz = histogram(soz_vals, bin_edges, 'Normalization', 'pdf', ...
        'FaceColor', c_soz, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'SOZ');
    h_non = histogram(non_vals, bin_edges, 'Normalization', 'pdf', ...
        'FaceColor', c_non, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Non-SOZ');
    
    xlim([min_v max_v]);
    max_density = max([h_soz.Values, h_non.Values]);
    if isnan(max_density) || max_density == 0, max_density = 1; end
    ylim([0 max_density * 1.15]);
    
    line([median_soz median_soz], ylim, 'Color', c_soz * 0.8, 'LineWidth', 2.0, 'LineStyle', '--', 'HandleVisibility', 'off');
    line([median_non median_non], ylim, 'Color', c_non * 0.8, 'LineWidth', 2.0, 'LineStyle', '--', 'HandleVisibility', 'off');
    
    title(p_title, 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Absolute Normative Z-Score (|Z|)');
    ylabel('Probability Density');
    grid on; set(gca, 'FontSize', 9);
    
    if p_idx == 1, legend({'SOZ', 'Non-SOZ'}, 'Location', 'northwest'); end
end

save_comb_base = fullfile(fig_output_dir, 'Aggregated_AllPowCoh_vs_TopHCTSA_Comparison');
exportgraphics(fig_comb, [save_comb_base '.png'], 'Resolution', 300);
exportgraphics(fig_comb, [save_comb_base '.pdf'], 'ContentType', 'vector');
fprintf('Saved Aggregated Comparison Figure (PNG & PDF): %s\n\nAll tasks complete.\n', save_comb_base);
close(fig_comb);

% =========================================================================
% --- HELPER FUNCTIONS ---
% =========================================================================
function p = nanpercentile(x, prc)
    x = x(~isnan(x)); if isempty(x), p = 0; return; end
    x = sort(x(:)); n = length(x);
    idx = (prc/100) * (n-1) + 1;
    if idx <= 1, p = x(1); elseif idx >= n, p = x(n);
    else, i = floor(idx); f = idx - i; p = x(i) + f * (x(i+1) - x(i)); end
end