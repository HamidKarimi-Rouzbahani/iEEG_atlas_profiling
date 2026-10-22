clc; clear; close all;

% ==========================================
% --- SECTION 0: PATH CONFIGURATIONS & DATA LOADING ---
% ==========================================
base_dir = '/scratch/user/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';

results_paths.PowCoh = fullfile(base_dir, 'AcrossSubjects/PowerCoh_Results/PowerCoherence/');
results_paths.HCTSA  = fullfile(base_dir, 'AcrossSubjects/HCTSA_Results/');

combined_output_dir = fullfile(base_dir, 'AcrossSubjects/Combined_Plots');
if ~exist(combined_output_dir, 'dir'), mkdir(combined_output_dir); end

feature_modes = {'PowCoh', 'HCTSA'};

% Define consistent color scheme for both feature pipelines
feature_colors.PowCoh = [0.8500, 0.3250, 0.0980]; % Terracotta Red/Orange
feature_colors.HCTSA  = [0.0000, 0.4470, 0.7410]; % Blue

data_loaded = struct();

% --- Load Both Feature Sets ---
for m = 1:length(feature_modes)
    mode = feature_modes{m};
    results_file = fullfile(results_paths.(mode), 'AcrossSubjects_Final_Results_both_feature_sets_Repeated_fixed_resect.mat');
    
    if ~exist(results_file, 'file')
        error('Results file not found for %s at %s.', mode, results_file);
    end
    
    fprintf('Loading results for %s...\n', mode);
    loaded_data = load(results_file);
    
    % Unify field name expectations
    if isfield(loaded_data, 'all_test_labels'), data_loaded.(mode).labels = loaded_data.all_test_labels;
    else, data_loaded.(mode).labels = loaded_data.Y_full; end
    
    if isfield(loaded_data, 'all_test_scores'), data_loaded.(mode).scores = loaded_data.all_test_scores;
    else, error('Scores not found for %s', mode); end
    
    if isfield(loaded_data, 'all_test_patient_ids'), data_loaded.(mode).p_ids = loaded_data.all_test_patient_ids;
    else, data_loaded.(mode).p_ids = loaded_data.patient_ids_full; end
    
    if isfield(loaded_data, 'all_test_clinical_info'), data_loaded.(mode).clinical_info = loaded_data.all_test_clinical_info;
    else, data_loaded.(mode).clinical_info = []; end
    
    % Compute Subject-Level AUCs & Extract Continuous Metadata
    unique_pats = unique(data_loaded.(mode).p_ids);
    sub_aucs    = [];
    sub_ids     = {};
    sub_age     = [];
    sub_onset   = [];
    sub_dur     = [];
    
    clin_tb = data_loaded.(mode).clinical_info;
    has_clin = ~isempty(clin_tb) && istable(clin_tb);
    
    for i = 1:length(unique_pats)
        p_mask = strcmp(data_loaded.(mode).p_ids, unique_pats{i});
        if length(unique(data_loaded.(mode).labels(p_mask))) >= 2
            [~,~,~,pat_auc] = perfcurve(data_loaded.(mode).labels(p_mask), data_loaded.(mode).scores(p_mask), 1);
            sub_aucs = [sub_aucs; pat_auc];
            sub_ids{end+1, 1} = unique_pats{i};
            
            % Extract subject continuous metrics
            if has_clin
                idx = find(p_mask, 1, 'first');
                
                c_age   = NaN; 
                c_onset = NaN; 
                c_dur   = NaN;
                
                if ismember('age', clin_tb.Properties.VariableNames)
                    c_age = double(clin_tb.age(idx));
                end
                if ismember('age_onset', clin_tb.Properties.VariableNames)
                    c_onset = double(clin_tb.age_onset(idx));
                end
                if ismember('duration', clin_tb.Properties.VariableNames)
                    c_dur = double(clin_tb.duration(idx));
                elseif ~isnan(c_age) && ~isnan(c_onset)
                    c_dur = c_age - c_onset;
                end
                
                sub_age   = [sub_age; c_age];
                sub_onset = [sub_onset; c_onset];
                sub_dur   = [sub_dur; c_dur];
            end
        end
    end
    
    data_loaded.(mode).subject_aucs      = sub_aucs;
    data_loaded.(mode).subject_ids       = sub_ids;
    data_loaded.(mode).subject_age       = sub_age;
    data_loaded.(mode).subject_age_onset = sub_onset;
    data_loaded.(mode).subject_duration  = sub_dur;
end

% STATISTICAL COMPARISON OF SUBJECT AUCs (PowCoh vs HCTSA)
fprintf('\n=========================================================\n');
fprintf('  BAYESIAN STATISTICAL COMPARISON: PowCoh vs HCTSA Subject AUCs\n');
fprintf('=========================================================\n');

aucs_pow   = data_loaded.PowCoh.subject_aucs;
aucs_hctsa = data_loaded.HCTSA.subject_aucs;

med_pow   = median(aucs_pow);
med_hctsa = median(aucs_hctsa);

fprintf('PowCoh: Median Subject AUC = %.4f (N = %d)\n', med_pow, length(aucs_pow));
fprintf('HCTSA : Median Subject AUC = %.4f (N = %d)\n', med_hctsa, length(aucs_hctsa));

if length(aucs_pow) == length(aucs_hctsa)
    diff_aucs = aucs_pow - aucs_hctsa;
    [bf_comp, log10_bf_comp] = compute_bayes_factor_one_sample(diff_aucs, 'two-sided');
    fprintf('Test Method: Paired JZS Bayes Factor\n');
    fprintf('BF10 = %10.2e | log10(BF10) = %5.2f [%s]\n', bf_comp, log10_bf_comp, get_bf_category(bf_comp));
else
    diff_aucs = aucs_pow - mean(aucs_hctsa);
    [bf_comp, log10_bf_comp] = compute_bayes_factor_one_sample(diff_aucs, 'two-sided');
    fprintf('Test Method: Two-Sample JZS Bayes Factor\n');
    fprintf('BF10 = %10.2e | log10(BF10) = %5.2f [%s]\n', bf_comp, log10_bf_comp, get_bf_category(bf_comp));
end
fprintf('=========================================================\n\n');

% Expanded Y-Axis limits for visual breathing room and proper whisker displays
common_ylim   = [0, 1.1];
common_yticks = 0.2:0.2:1.0;

% ==========================================
% --- SECTION 1: COMBINED 2x2 GLOBAL PERFORMANCE FIGURE ---
% ==========================================
fprintf('Generating Combined 2x2 Global Performance Figure...\n');
fig_global = figure('Color', 'w', 'Name', 'Combined_Global_Performance', 'Position', [100, 100, 1000, 800]);

% --- Row 1: PowCoh (ROC + Boxplot) ---
subplot(2, 2, 1);
[X_roc1, Y_roc1, ~, overall_auc1] = perfcurve(data_loaded.PowCoh.labels, data_loaded.PowCoh.scores, 1);
plot(X_roc1, Y_roc1, 'LineWidth', 2, 'Color', feature_colors.PowCoh); hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1.5);
xlabel('FPR'); ylabel('TPR');
title(sprintf('PowCoh: Global ROC (AUC = %.4f)', overall_auc1));
grid on; axis square;

subplot(2, 2, 3);
sub_aucs1 = data_loaded.PowCoh.subject_aucs;
h_box1 = boxplot(sub_aucs1, 'Symbol', '', 'Widths', 0.3); hold on;
color_boxplot_gray(h_box1);

jitter1 = (rand(size(sub_aucs1)) - 0.5) * 0.15;
scatter(1 + jitter1, sub_aucs1, 35, feature_colors.PowCoh, 'filled', 'MarkerFaceAlpha', 0.6);
plot([0.5, 1.5], [0.5, 0.5], 'k--', 'LineWidth', 1.5);

[bf_chan1, log10_bf1] = compute_bayes_factor_one_sample(sub_aucs1 - 0.50, 'one-sided');
if bf_chan1 > 100, txt_lbl1 = sprintf('%.1f', log10_bf1); else, txt_lbl1 = sprintf('%.1f', bf_chan1); end
if bf_chan1 > 6, fw1 = 'bold'; fc1 = [0.8 0 0]; else, fw1 = 'normal'; fc1 = [0.2 0.2 0.2]; end

text(1, 1.05, txt_lbl1, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 11, 'FontWeight', fw1, 'Color', fc1);

xlim([0.5, 1.5]); ylim(common_ylim); set(gca, 'YTick', common_yticks);
set(gca, 'XTick', 1, 'XTickLabel', {''}, 'FontSize', 10);
ylabel('Patient AUC');
title(sprintf('PowCoh: Subject-Averaged AUC\n(Median = %.4f)', median(sub_aucs1)));
grid on;

% --- Row 2: HCTSA (ROC + Boxplot) ---
subplot(2, 2, 2);
[X_roc2, Y_roc2, ~, overall_auc2] = perfcurve(data_loaded.HCTSA.labels, data_loaded.HCTSA.scores, 1);
plot(X_roc2, Y_roc2, 'LineWidth', 2, 'Color', feature_colors.HCTSA); hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1.5);
xlabel('FPR'); ylabel('TPR');
title(sprintf('HCTSA: Global ROC (AUC = %.4f)', overall_auc2));
grid on; axis square;

subplot(2, 2, 4);
sub_aucs2 = data_loaded.HCTSA.subject_aucs;
h_box2 = boxplot(sub_aucs2, 'Symbol', '', 'Widths', 0.3); hold on;
color_boxplot_gray(h_box2);

jitter2 = (rand(size(sub_aucs2)) - 0.5) * 0.15;
scatter(1 + jitter2, sub_aucs2, 35, feature_colors.HCTSA, 'filled', 'MarkerFaceAlpha', 0.6);
plot([0.5, 1.5], [0.5, 0.5], 'k--', 'LineWidth', 1.5);

[bf_chan2, log10_bf2] = compute_bayes_factor_one_sample(sub_aucs2 - 0.50, 'one-sided');
if bf_chan2 > 100, txt_lbl2 = sprintf('%.1f', log10_bf2); else, txt_lbl2 = sprintf('%.1f', bf_chan2); end
if bf_chan2 > 6, fw2 = 'bold'; fc2 = [0.8 0 0]; else, fw2 = 'normal'; fc2 = [0.2 0.2 0.2]; end

text(1, 1.05, txt_lbl2, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 11, 'FontWeight', fw2, 'Color', fc2);

xlim([0.5, 1.5]); ylim(common_ylim); set(gca, 'YTick', common_yticks);
set(gca, 'XTick', 1, 'XTickLabel', {''}, 'FontSize', 10);
ylabel('Patient AUC');
title(sprintf('HCTSA: Subject-Averaged AUC\n(Median = %.4f)', median(sub_aucs2)));
grid on;

png_global_file = fullfile(combined_output_dir, 'Combined_2x2_Global_Performance.png');
pdf_global_file = fullfile(combined_output_dir, 'Combined_2x2_Global_Performance.pdf');
saveas(fig_global, png_global_file);
exportgraphics(fig_global, pdf_global_file, 'ContentType', 'vector');
close(fig_global);

% ==========================================
% --- SECTION 2: CLINICAL VARIABLES SIDE-BY-SIDE PLOTS WITH BAYES FACTORS ---
% ==========================================
fprintf('\nGenerating Combined Clinical Breakdown Plots with Bayes Factors...\n');
predictors = {'duration', 'sex', 'hand', 'outcome', 'engel', 'therapy', 'implant', 'target', 'lesion_status'};

for p = 1:length(predictors)
    var_name = predictors{p};
    
    has_pow = ~isempty(data_loaded.PowCoh.clinical_info) && istable(data_loaded.PowCoh.clinical_info);
    has_hctsa = ~isempty(data_loaded.HCTSA.clinical_info) && istable(data_loaded.HCTSA.clinical_info);
    if ~has_pow && ~has_hctsa, continue; end
    
    fig_clin = figure('Color', 'w', 'Name', ['Clinical_', var_name], 'Position', [100, 100, 1200, 550]);
    
    for m = 1:length(feature_modes)
        mode = feature_modes{m};
        clin_info = data_loaded.(mode).clinical_info;
        labels = data_loaded.(mode).labels;
        scores = data_loaded.(mode).scores;
        p_ids = data_loaded.(mode).p_ids;
        current_color = feature_colors.(mode);
        
        if isempty(clin_info) || ~istable(clin_info), continue; end
        
        if strcmp(var_name, 'duration')
            if ~ismember('age', clin_info.Properties.VariableNames) || ...
               ~ismember('age_onset', clin_info.Properties.VariableNames)
                continue; 
            end
            raw_data = clin_info.age - clin_info.age_onset;
        else
            if ~ismember(var_name, clin_info.Properties.VariableNames), continue; end
            raw_data = clin_info.(var_name);
        end
        
        data_for_grouping = cell(size(raw_data));
        
        if strcmp(var_name, 'duration')
            num_vals = raw_data(~isnan(raw_data));
            if isempty(num_vals), continue; end
            bin_edges = floor(min(num_vals)/10)*10 : 10 : ceil(max(num_vals)/10)*10 + 10;
            for k = 1:length(raw_data)
                if isnumeric(raw_data(k)) && ~isnan(raw_data(k))
                    b_idx = discretize(raw_data(k), bin_edges);
                    if ~isnan(b_idx), data_for_grouping{k} = sprintf('%d-%d yrs', bin_edges(b_idx), bin_edges(b_idx+1)-1);
                    else, data_for_grouping{k} = 'unknown'; end
                else, data_for_grouping{k} = 'unknown'; end
            end
        else
            for k = 1:length(raw_data)
                val = raw_data(k);
                if ismissing(val) || isempty(val), data_for_grouping{k} = 'unknown';
                else, data_for_grouping{k} = lower(char(string(val))); end
            end
        end
        
        unique_cats = unique(data_for_grouping);
        unique_cats = unique_cats(~strcmp(unique_cats, 'unknown'));
        
        valid_cats = {};
        cat_aucs_list = {};
        
        for c = 1:length(unique_cats)
            cat = unique_cats{c};
            mask = strcmp(data_for_grouping, cat);
            subs_in_cat = unique(p_ids(mask));
            
            if length(subs_in_cat) < 3, continue; end
            
            per_pat_aucs = [];
            for si = 1:length(subs_in_cat)
                p_mask = mask & strcmp(p_ids, subs_in_cat(si));
                if length(unique(labels(p_mask))) >= 2
                    [~,~,~,sa] = perfcurve(labels(p_mask), scores(p_mask), 1);
                    per_pat_aucs = [per_pat_aucs; sa];
                end
            end
            
            if length(per_pat_aucs) >= 3
                valid_cats{end+1} = cat;
                cat_aucs_list{end+1} = per_pat_aucs;
            end
        end
        
        subplot(1, 2, m);
        if ~isempty(valid_cats)
            all_plot_vals = [];
            group_ids = [];
            for gc = 1:length(cat_aucs_list)
                all_plot_vals = [all_plot_vals; cat_aucs_list{gc}];
                group_ids = [group_ids; repmat(gc, length(cat_aucs_list{gc}), 1)];
            end
            
            h_box = boxplot(all_plot_vals, group_ids, 'Labels', valid_cats, 'Symbol', '', 'Widths', 0.4); hold on;
            color_boxplot_gray(h_box);
            
            for gc = 1:length(cat_aucs_list)
                y_vals = cat_aucs_list{gc};
                x_jit = gc + (rand(size(y_vals)) - 0.5) * 0.18;
                scatter(x_jit, y_vals, 30, current_color, 'filled', 'MarkerFaceAlpha', 0.7);
            end
            
            plot(xlim, [0.5 0.5], 'k--', 'LineWidth', 1.5);
            set(gca, 'XTickLabel', valid_cats, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
            
            for gc = 1:length(valid_cats)
                curr_aucs = cat_aucs_list{gc};
                max_val   = max(curr_aucs);
                
                [bf_gc, log10_bf_gc] = compute_bayes_factor_one_sample(curr_aucs - 0.50, 'one-sided');
                
                if bf_gc > 100, txt = sprintf('%.1f', log10_bf_gc);
                else, txt = sprintf('%.1f', bf_gc); end
                
                if bf_gc > 6, f_weight = 'bold'; f_color  = [0.8 0 0];
                else, f_weight = 'normal'; f_color  = [0.2 0.2 0.2]; end
                
                text(gc, min(1.05, max_val + 0.03), txt, 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'bottom', 'FontSize', 9, 'FontWeight', f_weight, 'Color', f_color);
            end
            
            n_groups = length(valid_cats);
            sig_count = 0;
            if n_groups >= 2
                base_y = min(1.08, max(all_plot_vals) + 0.05);
                step_y = 0.05;
                
                for g1 = 1:(n_groups - 1)
                    for g2 = (g1 + 1):n_groups
                        v1 = cat_aucs_list{g1};
                        v2 = cat_aucs_list{g2};
                        
                        if length(v1) == length(v2)
                            [bf_pair, log10_bf_pair] = compute_bayes_factor_one_sample(v1 - v2, 'two-sided');
                        else
                            [bf_pair, log10_bf_pair] = compute_bayes_factor_one_sample(v1 - mean(v2), 'two-sided');
                        end
                        
                        if bf_pair > 3
                            sig_count = sig_count + 1;
                            y_pos = base_y + (sig_count * step_y);
                            
                            if bf_pair > 100, l_width = 2.2;
                            elseif bf_pair > 10, l_width = 1.6;
                            elseif bf_pair > 6, l_width = 1.2;
                            else, l_width = 0.8; end
                            
                            plot([g1, g2], [y_pos, y_pos], 'k-', 'LineWidth', l_width);
                            
                            if bf_pair > 100, txt_pair = sprintf('%.1f', log10_bf_pair);
                            else, txt_pair = sprintf('%.1f', bf_pair); end
                            
                            fw_pair = ternary(bf_pair > 6, 'bold', 'normal');
                            
                            text((g1 + g2)/2, y_pos + 0.01, txt_pair, ...
                                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                                'FontSize', 9, 'FontWeight', fw_pair);
                        end
                    end
                end
            end
        else
            text(0.5, 0.5, 'Insufficient Multi-Class Data', 'HorizontalAlignment', 'center');
        end
        
        title(sprintf('%s: %s Breakdown', mode, strrep(var_name, '_', ' ')));
        ylabel('Patient AUC'); 
        ylim(common_ylim); set(gca, 'YTick', common_yticks);
        grid on;
    end
    
    png_clin_file = fullfile(combined_output_dir, ['Combined_AUC_', var_name, '.png']);
    pdf_clin_file = fullfile(combined_output_dir, ['Combined_AUC_', var_name, '.pdf']);
    saveas(fig_clin, png_clin_file);
    exportgraphics(fig_clin, pdf_clin_file, 'ContentType', 'vector');
    close(fig_clin);
end

% =========================================================================
% --- SECTION 3: EPILEPSY DURATION CORRELATION (1x2 Plot) ---
% =========================================================================
fprintf('\nGenerating Epilepsy Duration Scatter Plots with Spearman Correlation...\n');

fig_corr = figure('Color', 'w', 'Name', 'Epilepsy_Duration_Correlation', 'Position', [100, 100, 800, 400]);

for m = 1:length(feature_modes)
    mode = feature_modes{m};
    cur_color = feature_colors.(mode);
    sub_aucs  = data_loaded.(mode).subject_aucs;
    x_vals    = data_loaded.(mode).subject_duration;
    
    subplot(1, 2, m); hold on;
    
    valid_mask = ~isnan(x_vals) & ~isnan(sub_aucs);
    x_clean = x_vals(valid_mask);
    y_clean = sub_aucs(valid_mask);
    
    if length(x_clean) >= 3
        scatter(x_clean, y_clean, 38, cur_color, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'w');
        
        p_fit = polyfit(x_clean, y_clean, 1);
        x_grid = linspace(min(x_clean), max(x_clean), 100);
        y_grid = polyval(p_fit, x_grid);
        plot(x_grid, y_grid, 'Color', cur_color * 0.65, 'LineWidth', 2);
        
        [r_val, p_val] = corr(x_clean, y_clean, 'Type', 'Spearman');
        
        title(sprintf('%s: Disease Duration\nr_s = %.3f (p = %.3f)', mode, r_val, p_val), 'FontSize', 10, 'FontWeight', 'bold');
    else
        text(0.5, 0.5, 'Insufficient Data', 'HorizontalAlignment', 'center');
        title(sprintf('%s: Disease Duration', mode), 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    xlabel('Disease Duration (years)');
    ylabel('Patient AUC');
    ylim(common_ylim); set(gca, 'YTick', common_yticks);
    grid on; set(gca, 'FontSize', 9);
end

fig_corr_file = fullfile(combined_output_dir, 'Combined_Duration_Correlation.png');
fig_corr_pdf  = fullfile(combined_output_dir, 'Combined_Duration_Correlation.pdf');
saveas(fig_corr, fig_corr_file);
exportgraphics(fig_corr, fig_corr_pdf, 'ContentType', 'vector');
close(fig_corr);
fprintf('Saved Spearman correlation plot: %s\n', fig_corr_file);

fprintf('\nDone! Script completed.\n');

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================
function color_boxplot_gray(h_box)
    gray_color = [0.4, 0.4, 0.4];
    for i = 1:size(h_box, 1)
        for j = 1:size(h_box, 2)
            set(h_box(i, j), 'Color', gray_color, 'LineWidth', 1.2);
        end
    end
end

function [bf10, log10_bf10] = compute_bayes_factor_one_sample(x, test_type)
    x = x(~isnan(x)); N = length(x);
    if N < 2 || std(x) == 0, bf10 = 1; log10_bf10 = 0; return; end
    t_stat = mean(x) / (std(x) / sqrt(N)); r = sqrt(2) / 2;
    integral_fun = @(g) (1 + N.*g.*(r^2)).^(-0.5) .* (1 + (t_stat.^2) ./ ((1 + N.*g.*(r^2)).*(N - 1))).^(-N/2) .* (g.^(-1.5)) .* exp(-1 ./ (2.*g));
    try, marg_lik_h1 = integral(integral_fun, 0, Inf, 'ArrayValued', true); catch, marg_lik_h1 = 1e-300; end
    marg_lik_h0 = (1 + (t_stat^2)/(N - 1))^(-N/2);
    bf10_twosided = marg_lik_h1 / marg_lik_h0;
    bf10 = ternary(strcmpi(test_type, 'one-sided'), ternary(t_stat > 0, bf10_twosided * 2, bf10_twosided / 2), bf10_twosided);
    bf10 = max(bf10, 1e-300); log10_bf10 = log10(bf10);
end

function cat_str = get_bf_category(bf)
    if bf >= 100, cat_str = 'Decisive Evidence';
    elseif bf >= 30, cat_str = 'Very Strong Evidence';
    elseif bf >= 10, cat_str = 'Strong Evidence';
    elseif bf >= 3, cat_str = 'Moderate Evidence';
    elseif bf >= 1, cat_str = 'Anecdotal / Weak Evidence';
    else, cat_str = 'Evidence for Null (H0)'; end
end

function val = ternary(cond, true_val, false_val)
    if cond, val = true_val; else, val = false_val; end
end