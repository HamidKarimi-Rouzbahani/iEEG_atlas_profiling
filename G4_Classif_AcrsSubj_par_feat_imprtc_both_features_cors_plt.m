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

data_loaded = struct();

% Shared Color Schemes
color_black    = [0.0, 0.0, 0.0];
color_darkgray = [0.2, 0.2, 0.2];
color_midgray  = [0.5, 0.5, 0.5];

% Common Y-Axis limits for consistency
common_ylim   = [0.15, 1.25];
common_yticks = 0.2:0.2:1.0;

% --- Load Both Feature Sets ---
for m = 1:length(feature_modes)
    mode = feature_modes{m};
    results_file = fullfile(results_paths.(mode), 'AcrossSubjects_Final_Results_both_feature_sets_Repeated_fixed_resect.mat');
    
    if ~exist(results_file, 'file')
        error('Results file not found for %s at %s.', mode, results_file);
    end
    
    fprintf('Loading results for %s...\n', mode);
    loaded_data = load(results_file);
    
    if isfield(loaded_data, 'all_test_labels'), data_loaded.(mode).labels = loaded_data.all_test_labels;
    else, data_loaded.(mode).labels = loaded_data.Y_full; end
    
    if isfield(loaded_data, 'all_test_scores'), data_loaded.(mode).scores = loaded_data.all_test_scores;
    else, error('Scores not found for %s', mode); end
    
    if isfield(loaded_data, 'all_test_patient_ids'), data_loaded.(mode).p_ids = loaded_data.all_test_patient_ids;
    else, data_loaded.(mode).p_ids = loaded_data.patient_ids_full; end
    
    if isfield(loaded_data, 'all_test_clinical_info'), data_loaded.(mode).clinical_info = loaded_data.all_test_clinical_info;
    else, data_loaded.(mode).clinical_info = []; end
    
    unique_pats = unique(data_loaded.(mode).p_ids);
    sub_aucs    = [];
    sub_ids     = {};
    
    for i = 1:length(unique_pats)
        p_mask = strcmp(data_loaded.(mode).p_ids, unique_pats{i});
        if length(unique(data_loaded.(mode).labels(p_mask))) >= 2
            [~,~,~,pat_auc] = perfcurve(data_loaded.(mode).labels(p_mask), data_loaded.(mode).scores(p_mask), 1);
            sub_aucs = [sub_aucs; pat_auc];
            sub_ids{end+1, 1} = unique_pats{i};
        end
    end
    
    data_loaded.(mode).subject_aucs = sub_aucs;
    data_loaded.(mode).subject_ids  = sub_ids;
end

% Match subject IDs between PowCoh and HCTSA pipelines
[common_pts, idx_pow, idx_hctsa] = intersect(data_loaded.PowCoh.subject_ids, data_loaded.HCTSA.subject_ids, 'stable');

% =========================================================================
% --- SECTION 1: POWCOH vs. HCTSA SUBJECT AUC CORRELATION ---
% =========================================================================
fprintf('\nGenerating PowCoh vs HCTSA Patient AUC Correlation Scatter Plot...\n');

if ~isempty(common_pts) && length(common_pts) >= 3
    x_auc = data_loaded.PowCoh.subject_aucs(idx_pow);
    y_auc = data_loaded.HCTSA.subject_aucs(idx_hctsa);
    
    [r_spearman, p_spearman] = corr(x_auc, y_auc, 'Type', 'Spearman');
    
    fig_auc_corr = figure('Color', 'w', 'Name', 'PowCoh_vs_HCTSA_Patient_AUC_Correlation', 'Position', [150, 150, 650, 550]);
    hold on;
    
    plot([0.15, 1.25], [0.15, 1.25], 'k--', 'LineWidth', 1.2, 'DisplayName', 'Identity Line (1:1)');
    scatter(x_auc, y_auc, 45, color_black, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'k', 'DisplayName', 'Patients');
    
    p_fit = polyfit(x_auc, y_auc, 1);
    x_grid = linspace(min(x_auc), max(x_auc), 100);
    y_grid = polyval(p_fit, x_grid);
    plot(x_grid, y_grid, 'Color', color_midgray, 'LineWidth', 2.2, 'DisplayName', 'Linear Fit');
    
    xlabel('PowCoh Patient AUC', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('HCTSA Patient AUC', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('Patient AUC Concordance: PowCoh vs. HCTSA\nr_s = %.3f (p_s = %.3e)', ...
          r_spearman, p_spearman), 'FontSize', 11, 'FontWeight', 'bold');
        
    xlim(common_ylim); set(gca, 'XTick', common_yticks);
    ylim(common_ylim); set(gca, 'YTick', common_yticks);
    grid on; legend('Location', 'northwest'); set(gca, 'FontSize', 10);
    
    stats_str = sprintf('Spearman r_s: %.3f (p = %.3e)', ...
                        r_spearman, p_spearman);
    text(0.95, 0.08, stats_str, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'bottom', 'FontSize', 9.5, 'BackgroundColor', [1 1 1 0.85], 'EdgeColor', [0.7 0.7 0.7]);
       
    fig_auc_corr_png = fullfile(combined_output_dir, 'PowCoh_vs_HCTSA_Subject_AUC_Scatter.png');
    fig_auc_corr_pdf = fullfile(combined_output_dir, 'PowCoh_vs_HCTSA_Subject_AUC_Scatter.pdf');
    saveas(fig_auc_corr, fig_auc_corr_png);
    exportgraphics(fig_auc_corr, fig_auc_corr_pdf, 'ContentType', 'vector');
    close(fig_auc_corr);
    
    fprintf('=========================================================\n');
    fprintf('  PATIENT-LEVEL CONCORDANCE: PowCoh vs HCTSA AUCs\n');
    fprintf('=========================================================\n');
    fprintf('Matched Cohort N = %d Patients\n', length(common_pts));
    fprintf('Spearman : r_s = %.4f, p = %.4e\n', r_spearman, p_spearman);
    fprintf('Saved Correlation Figures: %s (.png & .pdf)\n', fig_auc_corr_png);
    fprintf('=========================================================\n');
else
    fprintf('Warning: Could not match patient IDs between PowCoh and HCTSA to generate correlation scatter.\n');
end

% =========================================================================
% --- SECTION 2: CLINICAL GROUP BREAKDOWN OF POWCOH VS HCTSA AUC CORRELATIONS ---
% =========================================================================
fprintf('\nGenerating PowCoh vs HCTSA AUC Correlations by Clinical Subgroups (Boxplots & Continuous Scatter)...\n');

if ~isempty(data_loaded.PowCoh.clinical_info) && istable(data_loaded.PowCoh.clinical_info) && ~isempty(common_pts)
    clin_info_master = data_loaded.PowCoh.clinical_info;
    p_ids_master     = data_loaded.PowCoh.p_ids;
    
    eval_predictors = {'duration', 'sex', 'hand', 'outcome', 'engel', 'therapy', 'implant', 'target', 'lesion_status'};
    
    % Part 1: Categorical Subgroup Correlation Boxplots
    for p = 1:length(eval_predictors)
        var_name = eval_predictors{p};
        
        if strcmp(var_name, 'duration')
            if ~ismember('age', clin_info_master.Properties.VariableNames) || ...
               ~ismember('age_onset', clin_info_master.Properties.VariableNames)
                continue; 
            end
            raw_data = clin_info_master.age - clin_info_master.age_onset;
        else
            if ~ismember(var_name, clin_info_master.Properties.VariableNames), continue; end
            raw_data = clin_info_master.(var_name);
        end
        
        pt_subgroups = cell(length(common_pts), 1);
        
        if strcmp(var_name, 'duration')
            num_vals = raw_data(~isnan(raw_data));
            if isempty(num_vals), continue; end
            bin_edges = floor(min(num_vals)/10)*10 : 10 : ceil(max(num_vals)/10)*10 + 10;
            for i = 1:length(common_pts)
                p_id = common_pts{i};
                idx = find(strcmp(p_ids_master, p_id), 1, 'first');
                if ~isempty(idx) && isnumeric(raw_data(idx)) && ~isnan(raw_data(idx))
                    b_idx = discretize(raw_data(idx), bin_edges);
                    if ~isnan(b_idx), pt_subgroups{i} = sprintf('%d-%d yrs', bin_edges(b_idx), bin_edges(b_idx+1)-1);
                    else, pt_subgroups{i} = 'unknown'; end
                else, pt_subgroups{i} = 'unknown'; end
            end
        else
            for i = 1:length(common_pts)
                p_id = common_pts{i};
                idx = find(strcmp(p_ids_master, p_id), 1, 'first');
                if ~isempty(idx)
                    val = raw_data(idx);
                    if ismissing(val) || isempty(val), pt_subgroups{i} = 'unknown';
                    else, pt_subgroups{i} = lower(char(string(val))); end
                else, pt_subgroups{i} = 'unknown'; end
            end
        end
        
        unique_groups = unique(pt_subgroups);
        unique_groups = unique_groups(~strcmp(unique_groups, 'unknown'));
        
        valid_groups = {};
        group_corrs_list = {};
        
        for g = 1:length(unique_groups)
            grp_name = unique_groups{g};
            grp_mask = strcmp(pt_subgroups, grp_name);
            
            x_grp = x_auc(grp_mask);
            y_grp = y_auc(grp_mask);
            
            if length(x_grp) >= 3
                n_pts_grp = length(x_grp);
                sub_corrs = zeros(n_pts_grp, 1);
                for jk = 1:n_pts_grp
                    jk_idx = true(n_pts_grp, 1); jk_idx(jk) = false;
                    sub_corrs(jk) = corr(x_grp(jk_idx), y_grp(jk_idx), 'Type', 'Spearman');
                end
                
                valid_groups{end+1} = grp_name;
                group_corrs_list{end+1} = sub_corrs;
            end
        end
        
        if isempty(valid_groups), continue; end
        
        fig_grp_corr = figure('Color', 'w', 'Name', ['Correlation_Boxplot_', var_name], 'Position', [100, 100, 700, 550]);
        hold on;
        
        all_plot_vals = [];
        group_ids = [];
        for gc = 1:length(group_corrs_list)
            all_plot_vals = [all_plot_vals; group_corrs_list{gc}];
            group_ids = [group_ids; repmat(gc, length(group_corrs_list{gc}), 1)];
        end
        
        h_box = boxplot(all_plot_vals, group_ids, 'Labels', valid_groups, 'Symbol', '', 'Widths', 0.4);
        color_boxplot_darkgray(h_box, color_darkgray);
        
        for gc = 1:length(group_corrs_list)
            y_vals = group_corrs_list{gc};
            x_jit = gc + (rand(size(y_vals)) - 0.5) * 0.18;
            scatter(x_jit, y_vals, 35, color_black, 'filled', 'MarkerFaceAlpha', 0.7);
        end
        
        % 0-Dashed reference line changed to BLACK
        plot(xlim, [0 0], 'k--', 'LineWidth', 1.5);
        
        % Center X-ticks cleanly relative to the data range
        xlim([0.5, length(valid_groups) + 0.5]);
        set(gca, 'XTick', 1:length(valid_groups), 'XTickLabel', valid_groups, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
        
        for gc = 1:length(valid_groups)
            curr_corrs = group_corrs_list{gc};
            max_val    = max(curr_corrs);
            
            [bf_gc, log10_bf_gc] = compute_bayes_factor_one_sample(curr_corrs, 'one-sided');
            
            if bf_gc > 100, txt = sprintf('%.1f', log10_bf_gc); else, txt = sprintf('%.1f', bf_gc); end
            if bf_gc > 6, f_weight = 'bold'; f_color = [0.8 0 0]; else, f_color = [0.2 0.2 0.2]; f_weight = 'normal'; end
            
            text(gc, min(1.20, max_val + 0.05), txt, 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', 'FontSize', 9, 'FontWeight', f_weight, 'Color', f_color);
        end
        
        n_groups = length(valid_groups);
        sig_count = 0;
        if n_groups >= 2
            base_y = max(all_plot_vals) + 0.15;
            step_y = 0.12;
            
            for g1 = 1:(n_groups - 1)
                for g2 = (g1 + 1):n_groups
                    v1 = group_corrs_list{g1};
                    v2 = group_corrs_list{g2};
                    
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
                        
                        text((g1 + g2)/2, y_pos + 0.02, txt_pair, ...
                            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                            'FontSize', 9, 'FontWeight', fw_pair);
                    end
                end
            end
        end
        
        title(sprintf('PowCoh vs HCTSA AUC Correlation by %s', strrep(var_name, '_', ' ')), 'FontWeight', 'bold');
        ylabel('Spearman Correlation (r_s)');
        ylim([-0.5, max(1.25, 0.25 + max(all_plot_vals) + (sig_count + 1) * 0.12)]);
        grid on;
        
        fig_grp_png = fullfile(combined_output_dir, ['Combined_Correlation_', var_name, '.png']);
        fig_grp_pdf = fullfile(combined_output_dir, ['Combined_Correlation_', var_name, '.pdf']);
        saveas(fig_grp_corr, fig_grp_png);
        exportgraphics(fig_grp_corr, fig_grp_pdf, 'ContentType', 'vector');
        close(fig_grp_corr);
        fprintf('Saved correlation breakdown plot: %s (.png & .pdf)\n', fig_grp_png);
    end

    % Part 2: Continuous Disease Duration Scatter Plot
    n_pts = length(common_pts);
    fig_cont_scatter = figure('Color', 'w', 'Name', 'Combined_Duration_Correlation_Scatter', 'Position', [100, 100, 600, 420]);
    
    x_vals = zeros(n_pts, 1);
    for i = 1:n_pts
        p_id = common_pts{i};
        idx = find(strcmp(p_ids_master, p_id), 1, 'first');
        if ~isempty(idx) && ismember('age', clin_info_master.Properties.VariableNames) && ismember('age_onset', clin_info_master.Properties.VariableNames)
            x_vals(i) = double(clin_info_master.age(idx) - clin_info_master.age_onset(idx));
        else, x_vals(i) = NaN; end
    end
    
    valid_mask = ~isnan(x_vals) & ~isnan(x_auc) & ~isnan(y_auc);
    x_clean = x_vals(valid_mask);
    x_a     = x_auc(valid_mask);
    y_a     = y_auc(valid_mask);
    
    [x_clean, sort_idx] = sort(x_clean);
    x_a = x_a(sort_idx);
    y_a = y_a(sort_idx);
    
    win_size = min(12, max(5, floor(length(x_clean)/3)));
    y_local_rs = zeros(length(x_clean), 1);
    
    for k = 1:length(x_clean)
        [~, dist_idx] = sort(abs(x_clean - x_clean(k)));
        w_idx = dist_idx(1:win_size);
        
        w_idx_loo = w_idx(w_idx ~= k);
        y_local_rs(k) = corr(x_a(w_idx_loo), y_a(w_idx_loo), 'Type', 'Spearman');
    end
    
    hold on;
    if length(x_clean) >= 3
        scatter(x_clean, y_local_rs, 40, color_black, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'w');
        
        p_fit = polyfit(x_clean, y_local_rs, 1);
        x_grid = linspace(min(x_clean), max(x_clean), 100);
        y_grid = polyval(p_fit, x_grid);
        plot(x_grid, y_grid, 'Color', color_midgray, 'LineWidth', 2.2);
        
        [r_val, p_val] = corr(x_clean, y_local_rs, 'Type', 'Spearman');
        
        title(sprintf('Disease Duration Concordance\nr_s = %.3f (p = %.3f)', r_val, p_val), 'FontSize', 11, 'FontWeight', 'bold');
    else
        text(0.5, 0.5, 'Insufficient Data', 'HorizontalAlignment', 'center');
        title('Disease Duration Concordance', 'FontSize', 11, 'FontWeight', 'bold');
    end
    
    xlabel('Disease Duration (years)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('PowCoh-HCTSA AUC Correlation (r_s)', 'FontSize', 10, 'FontWeight', 'bold');
    ylim([0.4, 1.1]); % Set ylim for duration plot from 0.4 to 1.1
    grid on; set(gca, 'FontSize', 10);
    
    fig_cont_png = fullfile(combined_output_dir, 'Combined_Duration_Correlation_Scatter.png');
    fig_cont_pdf = fullfile(combined_output_dir, 'Combined_Duration_Correlation_Scatter.pdf');
    saveas(fig_cont_scatter, fig_cont_png);
    exportgraphics(fig_cont_scatter, fig_cont_pdf, 'ContentType', 'vector');
    close(fig_cont_scatter);
    fprintf('Saved continuous disease duration correlation scatter figure: %s (.png & .pdf)\n', fig_cont_png);
end

fprintf('\nDone! Script finished. All plots saved in: %s\n', combined_output_dir);

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================
function color_boxplot_darkgray(h_box, box_color)
    for i = 1:size(h_box, 1)
        for j = 1:size(h_box, 2)
            set(h_box(i, j), 'Color', box_color, 'LineWidth', 1.2);
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

function val = ternary(cond, true_val, false_val)
    if cond, val = true_val; else, val = false_val; end
end