clc; clear; close all;

% ==========================================
% --- SECTION 0: PATH CONFIGURATIONS ---
% ==========================================
% base_dir = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';
base_dir = '/scratch/user/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';

% Results input directory for HFO Spikes
results_dir = fullfile(base_dir, 'AcrossSubjects/HFO_Spike_Results/');
results_file = fullfile(results_dir, 'AcrossSubjects_Final_Results_HFO_Spikes.mat');

if ~exist(results_file, 'file')
    error('Results file not found at %s. Please run the classification script first.', results_file);
end

% Plots output directory
output_dir = fullfile(base_dir, 'AcrossSubjects/HFO_Spike_Plots');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Define primary color for HFO Spikes (Purple/Indigo theme)
hfo_color = [0.4940, 0.1840, 0.5560]; 

fprintf('Loading HFO-Spike results from: %s...\n', results_file);
loaded_data = load(results_file);

% --- Ingest & Structure Data ---
data = struct();
data.labels = loaded_data.all_test_labels;
data.scores = loaded_data.all_test_scores;
data.p_ids = loaded_data.all_test_patient_ids;
data.clinical_info = loaded_data.all_test_clinical_info;

% Extract Subject-Level Metrics
unique_pats = unique(data.p_ids);
sub_aucs    = [];
sub_ids     = {};
sub_age     = [];
sub_onset   = [];
sub_dur     = [];

clin_tb = data.clinical_info;
has_clin = ~isempty(clin_tb) && istable(clin_tb);

for i = 1:length(unique_pats)
    p_mask = strcmp(data.p_ids, unique_pats{i});
    if length(unique(data.labels(p_mask))) >= 2
        [~,~,~,pat_auc] = perfcurve(data.labels(p_mask), data.scores(p_mask), 1);
        sub_aucs = [sub_aucs; pat_auc];
        sub_ids{end+1, 1} = unique_pats{i};
        
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

data.subject_aucs      = sub_aucs;
data.subject_ids       = sub_ids;
data.subject_age       = sub_age;
data.subject_age_onset = sub_onset;
data.subject_duration  = sub_dur;


% ==========================================
% --- SECTION 1: GLOBAL PERFORMANCE FIGURE ---
% ==========================================
fprintf('Generating Global Performance Figure...\n');
figure('Color', 'w', 'Name', 'HFO_Spike_Global_Performance', 'Position', [100, 100, 850, 380]);

% --- Plot 1: Global ROC Curve ---
subplot(1, 2, 1);
[X_roc, Y_roc, ~, overall_auc] = perfcurve(data.labels, data.scores, 1);
plot(X_roc, Y_roc, 'LineWidth', 2.5, 'Color', hfo_color); hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1.8);
xlabel('FPR', 'FontSize', 11); ylabel('TPR', 'FontSize', 11);
title(sprintf('HFO Spikes: Global ROC (AUC = %.4f)', overall_auc), 'FontSize', 12, 'FontWeight', 'bold');
grid on; axis square; set(gca, 'FontSize', 10, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
xlim([0 1]); ylim([0 1]);

% --- Plot 2: Subject-Averaged AUC Boxplot & Scatter ---
subplot(1, 2, 2);
boxplot(sub_aucs, 'Symbol', '', 'Widths', 0.35, 'Colors', [0.2 0.2 0.2]); hold on;
jitter = (rand(size(sub_aucs)) - 0.5) * 0.15;
scatter(1 + jitter, sub_aucs, 40, hfo_color, 'filled', 'MarkerFaceAlpha', 0.65);

% Black Dashed Chance Line
plot([0.5, 1.5], [0.5, 0.5], 'k--', 'LineWidth', 1.8);

% Above-Chance Bayes Factor (vs 0.50 Baseline)
[bf_chan, log10_bf] = compute_bayes_factor_one_sample(sub_aucs - 0.50, 'one-sided');
if bf_chan > 100, txt_lbl = sprintf('%.1f', log10_bf); else, txt_lbl = sprintf('%.1f', bf_chan); end
if bf_chan > 6, fw = 'bold'; fc = [0.8 0 0]; else, fw = 'normal'; fc = [0.2 0.2 0.2]; end

max_y = max(sub_aucs);
text(1, max_y + 0.05, txt_lbl, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 11, 'FontWeight', fw, 'Color', fc);

xlim([0.5, 1.5]); ylim([0.15, 1.2]);
set(gca, 'XTick', 1, 'XTickLabel', {'Cohort Patients'}, 'YTick', 0.2:0.2:1.0, 'FontSize', 10);
ylabel('Patient AUC', 'FontSize', 11);
title(sprintf('HFO Spikes: Subject-Averaged AUC\n(Median = %.4f)', median(sub_aucs)), ...
    'FontSize', 11, 'FontWeight', 'bold');
grid on;

fig1_file = fullfile(output_dir, 'HFO_Spikes_Global_Performance.png');
saveas(gcf, fig1_file);
close(gcf);
fprintf('Saved Global Performance Plot: %s\n', fig1_file);


% ==========================================
% --- SECTION 2: CATEGORICAL CLINICAL BREAKDOWN PLOTS ---
% ==========================================
fprintf('\nGenerating Categorical Clinical Breakdown Plots with Bayes Factors...\n');
predictors = {'duration', 'sex', 'hand', 'outcome', 'engel', 'therapy', 'implant', 'target', 'lesion_status'};

if has_clin
    for p = 1:length(predictors)
        var_name = predictors{p};
        
        % Dynamic Calculation for Duration
        if strcmp(var_name, 'duration')
            if ~ismember('age', clin_tb.Properties.VariableNames) || ...
               ~ismember('age_onset', clin_tb.Properties.VariableNames)
                continue; 
            end
            raw_data = clin_tb.age - clin_tb.age_onset;
        else
            if ~ismember(var_name, clin_tb.Properties.VariableNames), continue; end
            raw_data = clin_tb.(var_name);
        end
        
        data_for_grouping = cell(size(raw_data));
        
        % Binning / Categorization Logic
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
            subs_in_cat = unique(data.p_ids(mask));
            
            if length(subs_in_cat) < 3, continue; end
            
            per_pat_aucs = [];
            for si = 1:length(subs_in_cat)
                p_mask = mask & strcmp(data.p_ids, subs_in_cat(si));
                if length(unique(data.labels(p_mask))) >= 2
                    [~,~,~,sa] = perfcurve(data.labels(p_mask), data.scores(p_mask), 1);
                    per_pat_aucs = [per_pat_aucs; sa];
                end
            end
            
            if length(per_pat_aucs) >= 3
                valid_cats{end+1} = cat;
                cat_aucs_list{end+1} = per_pat_aucs;
            end
        end
        
        if isempty(valid_cats), continue; end
        
        fig_clin = figure('Color', 'w', 'Name', ['Clinical_', var_name], 'Position', [100, 100, 650, 500]);
        hold on;
        
        all_plot_vals = [];
        group_ids = [];
        for gc = 1:length(cat_aucs_list)
            all_plot_vals = [all_plot_vals; cat_aucs_list{gc}];
            group_ids = [group_ids; repmat(gc, length(cat_aucs_list{gc}), 1)];
        end
        
        boxplot(all_plot_vals, group_ids, 'Labels', valid_cats, 'Symbol', '', 'Widths', 0.4, 'Colors', [0.2 0.2 0.2]);
        
        % Scatter overlay
        for gc = 1:length(cat_aucs_list)
            y_vals = cat_aucs_list{gc};
            x_jit = gc + (rand(size(y_vals)) - 0.5) * 0.18;
            scatter(x_jit, y_vals, 35, hfo_color, 'filled', 'MarkerFaceAlpha', 0.7);
        end
        
        % Black Dashed Chance Line
        plot([0, length(valid_cats) + 1], [0.5 0.5], 'k--', 'LineWidth', 1.8);
        set(gca, 'XTickLabel', valid_cats, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
        
        % --- Above-Chance Bayes Factor (vs 0.50) ---
        for gc = 1:length(valid_cats)
            curr_aucs = cat_aucs_list{gc};
            max_val   = max(curr_aucs);
            [bf_gc, log10_bf_gc] = compute_bayes_factor_one_sample(curr_aucs - 0.50, 'one-sided');
            
            if bf_gc > 100, txt = sprintf('%.1f', log10_bf_gc); else, txt = sprintf('%.1f', bf_gc); end
            if bf_gc > 6, f_weight = 'bold'; f_color = [0.8 0 0]; else, f_weight = 'normal'; f_color = [0.2 0.2 0.2]; end
            
            text(gc, max_val + 0.04, txt, 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', f_weight, 'Color', f_color);
        end
        
        % --- Pairwise Category Bayes Factor Comparisons ---
        n_groups = length(valid_cats);
        sig_count = 0;
        if n_groups >= 2
            base_y = max(all_plot_vals) + 0.12;
            step_y = 0.08;
            
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
        
        title(sprintf('HFO Spikes: %s Breakdown', strrep(var_name, '_', ' ')), 'FontWeight', 'bold');
        ylabel('Patient AUC', 'FontSize', 11); 
        xlim([0.5, length(valid_cats) + 0.5]);
        ylim([0.15, max(1.20, 0.22 + max(all_plot_vals) + (sig_count + 1) * 0.08)]); 
        set(gca, 'YTick', 0.2:0.2:1.0, 'FontSize', 10);
        grid on;
        
        fig_clin_file = fullfile(output_dir, ['HFO_Spikes_AUC_', var_name, '.png']);
        saveas(fig_clin, fig_clin_file);
        close(fig_clin);
    end
end


% =========================================================================
% --- SECTION 3: CONTINUOUS CLINICAL CORRELATIONS (Disease Duration Only) ---
% =========================================================================
fprintf('\nGenerating Continuous Clinical Scatter Plot (Disease Duration Only)...\n');

fig_corr = figure('Color', 'w', 'Name', 'Continuous_Clinical_Correlation_Duration', 'Position', [100, 100, 500, 400]);

x_vals = data.subject_duration;
valid_mask = ~isnan(x_vals) & ~isnan(sub_aucs);
x_clean = x_vals(valid_mask);
y_clean = sub_aucs(valid_mask);

hold on;
if length(x_clean) >= 3
    scatter(x_clean, y_clean, 45, hfo_color, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'w');
    
    p_fit = polyfit(x_clean, y_clean, 1);
    x_grid = linspace(min(x_clean), max(x_clean), 100);
    y_grid = polyval(p_fit, x_grid);
    plot(x_grid, y_grid, 'Color', hfo_color * 0.65, 'LineWidth', 2);
    
    [r_val, p_val] = corr(x_clean, y_clean, 'Type', 'Spearman');
    
    title(sprintf('HFO Spikes: Disease Duration\nr_s = %.3f (p = %.3f)', r_val, p_val), 'FontSize', 11, 'FontWeight', 'bold');
else
    text(0.5, 0.5, 'Insufficient Data', 'HorizontalAlignment', 'center');
    title('HFO Spikes: Disease Duration', 'FontSize', 11, 'FontWeight', 'bold');
end

xlabel('Disease Duration (years)', 'FontSize', 10);
ylabel('Patient AUC', 'FontSize', 10);
ylim([0.15, 1.05]);
set(gca, 'YTick', 0.2:0.2:1.0, 'FontSize', 9);
grid on;

fig_corr_file = fullfile(output_dir, 'HFO_Spikes_Continuous_Clinical_Correlations.png');
saveas(fig_corr, fig_corr_file);
close(fig_corr);
fprintf('Saved Continuous Clinical Correlation Plot: %s\n', fig_corr_file);

fprintf('\nDone! All plots for HFO-Spikes saved in: %s\n', output_dir);


% =========================================================================
% HELPER FUNCTIONS FOR BAYESIAN STATISTICAL CALCULATIONS
% =========================================================================

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