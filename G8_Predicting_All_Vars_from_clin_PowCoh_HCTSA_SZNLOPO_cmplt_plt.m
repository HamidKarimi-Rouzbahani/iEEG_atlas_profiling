clc; clearvars; close all;

% =========================================================================
% CONFIGURATION: SELECT TARGETS & PATHS
% =========================================================================
base_folder = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';
output_dir  = fullfile(base_folder, 'Classification_Results_LOSO_cmplt/'); 
external_names_file = '/home/uqhkarim/Normative_HCTSA/Scripts/LittsScripts/Feature_names.mat';

freq_bands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
selected_targets = 'all'; 
chance_level = 50; % Baseline (%)

% Number of features to display in feature importance plots
k_show_default = 10;
k_show_hctsa   = 50; % Set to 15, 20, etc. for more HCTSA features

% Define Base Colors
c_clin   = [0.55, 0.25, 0.65]; % Purple
c_hctsa  = [0.4660, 0.6740, 0.1880]; % Green
c_powcoh = [0.8500, 0.3250, 0.0980]; % Orange

model_colors = containers.Map();
model_colors('Clin_Only')    = c_clin;
model_colors('HCTSA_Only')   = c_hctsa;
model_colors('PowCoh_Only')  = c_powcoh;
model_colors('HCTSA_PowCoh') = (c_hctsa + c_powcoh) / 2;
model_colors('Clin_HCTSA')   = (c_clin + c_hctsa) / 2;
model_colors('Clin_PowCoh')  = (c_clin + c_powcoh) / 2;
model_colors('All_Features') = (c_clin + c_hctsa + c_powcoh) / 3;

% Single-Domain Models selected for the Feature Importance plot
single_domain_models = {'Clin_Only', 'HCTSA_Only', 'PowCoh_Only'};

% Clinical Variable Order & Label Mapping
clin_vars_order = {'lesion_status', 'sex', 'target', 'therapy', 'hand', 'implant', 'age', 'age_onset', 'duration'};
clin_labels_map = containers.Map(...
    {'lesion_status', 'sex', 'target', 'therapy', 'hand', 'implant', 'age', 'age_onset', 'duration', 'engel'}, ...
    {'Lesion Status', 'Sex', 'Target', 'Therapy', 'Handedness', 'Implant Type', 'Age', 'Age Onset', 'Duration', 'Engel Score'} ...
);

% Load External HCTSA Feature Table for Lookup
hctsa_names_table = {};
if exist(external_names_file, 'file')
    fprintf('Loading HCTSA feature names lookup: %s\n', external_names_file);
    ext_data = load(external_names_file);
    if isfield(ext_data, 'Feature_names')
        FT = ext_data.Feature_names;
        col_to_use = 'Name';
        if ~ismember('Name', FT.Properties.VariableNames) && ismember('CodeString', FT.Properties.VariableNames)
            col_to_use = 'CodeString';
        end
        
        num_names = height(FT);
        hctsa_names_table = cell(num_names, 1);
        for i = 1:num_names
            val = FT.(col_to_use){i};
            if iscell(val), hctsa_names_table{i} = char(val{1});
            else, hctsa_names_table{i} = char(val); end
        end
    end
end

% =========================================================================
% 1. IDENTIFY TARGET FOLDERS TO PROCESS
% =========================================================================
if ischar(selected_targets) && strcmpi(selected_targets, 'all')
    dir_info = dir(output_dir);
    dir_info = dir_info([dir_info.isdir] & ~startsWith({dir_info.name}, '.'));
    target_folders = {dir_info.name};
else
    target_folders = selected_targets;
end

if isempty(target_folders)
    error('No output directories found matching selection in: %s', output_dir);
end

% Fixed Axis Scaling for Performance Plot
GLOBAL_Y_MIN  = 30;   
GLOBAL_Y_MAX  = 210;  
GLOBAL_TICKS  = 30:10:100; 
BASE_Y_BAR    = 105;  
STEP_Y_BAR    = 4.5;  

% =========================================================================
% 2. ITERATE, COMPUTE STATS, AND PLOT
% =========================================================================
for t = 1:length(target_folders)
    curr_target = target_folders{t};
    target_path = fullfile(output_dir, curr_target);
    
    if ~exist(target_path, 'dir'), continue; end

    mat_files = dir(fullfile(target_path, 'Classification_*_LOSO_Results_complete.mat'));
    if isempty(mat_files)
        mat_files = dir(fullfile(target_path, 'Classification_*_Results.mat'));
    end
    
    for m = 1:length(mat_files)
        mat_filepath = fullfile(target_path, mat_files(m).name);
        data = load(mat_filepath);
        if ~isfield(data, 'cv_results'), continue; end

        cv_results = data.cv_results;
        comp_name  = 'target';
        if isfield(data, 'comp_name'), comp_name = data.comp_name; end
        
        model_names = fieldnames(cv_results);
        num_models  = length(model_names);
        names_clean = strrep(model_names, '_', ' ');
        num_repeats = length(cv_results.(model_names{1}).raw_repeat_accs);
        
        raw_acc_matrix = zeros(num_repeats, num_models);
        mean_accs      = zeros(1, num_models);
        se_accs        = zeros(1, num_models);
        
        for i = 1:num_models
            accs = cv_results.(model_names{i}).raw_repeat_accs;
            raw_acc_matrix(:, i) = accs;
            mean_accs(i)         = mean(accs);
            se_accs(i)           = std(accs) / sqrt(length(accs));
        end
        
        % Bayesian Calculations
        bf_chance     = zeros(1, num_models);
        log10_bf_chan = zeros(1, num_models);
        for i = 1:num_models
            data_diff = raw_acc_matrix(:, i) - chance_level;
            [bf_chance(i), log10_bf_chan(i)] = compute_bayes_factor_one_sample(data_diff, 'one-sided');
        end
        
        bf_pair_mat  = zeros(num_models, num_models);
        log_pair_mat = zeros(num_models, num_models);
        for i = 1:(num_models - 1)
            for j = (i + 1):num_models
                diff_ij = raw_acc_matrix(:, i) - raw_acc_matrix(:, j);
                [bf_ij, log10_bf_ij] = compute_bayes_factor_one_sample(diff_ij, 'two-sided');
                bf_pair_mat(i, j)  = bf_ij;        bf_pair_mat(j, i)  = bf_ij;
                log_pair_mat(i, j) = log10_bf_ij; log_pair_mat(j, i) = log10_bf_ij;
            end
        end
        
        % -----------------------------------------------------------------
        % FIGURE 1: ROC & ACCURACY (SWAPPED PANELS)
        % -----------------------------------------------------------------
        fig_perf = figure('Name', [comp_name ' - Performance'], 'Color', 'w', ...
            'Position', [50 100 1300 520], 'Visible', 'on');
        
        % PANEL 1 (LEFT): ROC Curves
        subplot(1, 2, 1); hold on;
        legend_labels = cell(1, num_models);
        for i = 1:num_models
            m_key = model_names{i};
            m_color = ternary(isKey(model_colors, m_key), model_colors(m_key), [0.4 0.4 0.4]);
            m_struct = cv_results.(m_key);
            auc_val = NaN;
            if isfield(m_struct, 'roc') && ~isempty(m_struct.roc)
                plot(m_struct.roc(:, 1), m_struct.roc(:, 2), 'LineWidth', 2.2, 'Color', m_color);
                auc_val = ternary(isfield(m_struct, 'auc') && ~isempty(m_struct.auc), mean(m_struct.auc), trapz(m_struct.roc(:, 1), m_struct.roc(:, 2)));
            end
            legend_labels{i} = ternary(~isnan(auc_val), sprintf('%s (AUC = %.2f)', names_clean{i}, auc_val), names_clean{i});
        end
        plot([0 1], [0 1], 'k--', 'LineWidth', 1.2);
        title(sprintf('Aggregated ROC Curves\n[%s]', strrep(comp_name, '_', ' ')), 'Interpreter', 'none');
        xlabel('False Positive Rate'); ylabel('True Positive Rate');
        legend(legend_labels, 'Location', 'southeast', 'FontSize', 9);
        grid on; set(gca, 'FontSize', 10);
        
        % PANEL 2 (RIGHT): Accuracy Boxplots & Statistics
        subplot(1, 2, 2); hold on;
        plot(1:num_models, raw_acc_matrix', '-', 'Color', [0.82 0.82 0.82 0.5], 'LineWidth', 0.6);
        bh = boxplot(raw_acc_matrix, 'Colors', 'k', 'Symbol', '', 'Widths', 0.35);
        set(bh, 'LineWidth', 1.2);
        
        for i = 1:num_models
            m_key = model_names{i};
            m_color = ternary(isKey(model_colors, m_key), model_colors(m_key), [0.4 0.4 0.4]);
            x_jitter = i + (rand(num_repeats, 1) - 0.5) * 0.12;
            scatter(x_jitter, raw_acc_matrix(:, i), 28, m_color, 'filled', ...
                'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
        end
        
        yline(chance_level, '--r', 'Chance (50%)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left', 'FontSize', 9);
        
        max_y_per_model = max(raw_acc_matrix, [], 1);
        for i = 1:num_models
            bf_val  = bf_chance(i); log_val = log10_bf_chan(i);
            txt_label = sprintf('%.1f', ternary(bf_val > 100, log_val, bf_val));
            f_weight  = ternary(bf_val > 6, 'bold', 'normal');
            f_color   = ternary(bf_val > 6, [0.8 0 0], [0.2 0.2 0.2]);
            text(i, max_y_per_model(i) + 2.5, txt_label, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', 9, 'FontWeight', f_weight, 'Color', f_color);
        end
        
        sig_bar_count = 0;
        for i = 1:(num_models - 1)
            for j = (i + 1):num_models
                bf_val = bf_pair_mat(i, j); log_val = log_pair_mat(i, j);
                if bf_val > 3
                    sig_bar_count = sig_bar_count + 1;
                    y_pos = BASE_Y_BAR + (sig_bar_count * STEP_Y_BAR);
                    l_width = ternary(bf_val > 100, 2.2, ternary(bf_val > 10, 1.6, ternary(bf_val > 6, 1.2, 0.8)));
                    plot([i, j], [y_pos, y_pos], 'k-', 'LineWidth', l_width);
                    txt_label = sprintf('%.1f', ternary(bf_val > 100, log_val, bf_val));
                    text((i + j) / 2, y_pos + 0.5, txt_label, 'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'bottom', 'FontSize', 8.5, 'FontWeight', ternary(bf_val > 6, 'bold', 'normal'));
                end
            end
        end
        
        title(sprintf('Classification Accuracy Across Repeats\n[%s]', strrep(comp_name, '_', ' ')), 'Interpreter', 'none');
        ylabel('Accuracy (%)'); grid on;
        ylim([GLOBAL_Y_MIN, GLOBAL_Y_MAX]);
        set(gca, 'YTick', GLOBAL_TICKS, 'YTickLabel', arrayfun(@(v) sprintf('%d%%', v), GLOBAL_TICKS, 'UniformOutput', false));
        set(gca, 'XTick', 1:num_models, 'XTickLabel', names_clean, 'XTickLabelRotation', 30, 'FontSize', 10);
        
        % Save PNG & Vector PDF for Figure 1
        fig_perf_base = fullfile(target_path, sprintf('Classification_%s_LOSO_Accuracy_ROC', comp_name));
        saveas(fig_perf, [fig_perf_base '.png']);
        save_figure_as_pdf(fig_perf, [fig_perf_base '.pdf']);
        fprintf('Saved performance plot (PNG & PDF): %s\n', fig_perf_base);
        
        % -----------------------------------------------------------------
        % FIGURE 2: SINGLE DOMAIN FEATURE IMPORTANCES (3 PANELS ONLY)
        % -----------------------------------------------------------------
        num_sd_models = length(single_domain_models);
        fig_feats = figure('Name', [comp_name ' - Single Domain Features'], 'Color', 'w', ...
            'Position', [100 150 1250 560], 'Visible', 'on');
        
        sgtitle(sprintf('Top Predictive Features (Single Domain Models) [%s]', strrep(comp_name, '_', ' ')), ...
            'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
            
        fields_to_check = {'feature_importance', 'feature_importances', 'all_importances_matrix', ...
                           'importances', 'importance', 'weights', 'feat_imp', 'w'};
        
        for i = 1:num_sd_models
            m_key = single_domain_models{i};
            if ~isfield(cv_results, m_key), continue; end
            
            m_color = ternary(isKey(model_colors, m_key), model_colors(m_key), [0.4 0.4 0.4]);
            m_struct = cv_results.(m_key);
            
            subplot(1, num_sd_models, i); hold on;
            
            imp_matrix = [];
            for f_idx = 1:length(fields_to_check)
                fn = fields_to_check{f_idx};
                if isfield(m_struct, fn) && ~isempty(m_struct.(fn))
                    imp_matrix = m_struct.(fn); break;
                end
            end
            
            if isempty(imp_matrix)
                text(0.5, 0.5, sprintf('No Importances\nFound'), 'HorizontalAlignment', 'center', 'FontSize', 9);
                title(strrep(m_key, '_', ' '), 'FontSize', 10, 'Interpreter', 'none');
                axis off; continue;
            end
            
            imp_vec = ternary(isvector(imp_matrix), imp_matrix(:)', nanmean(imp_matrix, 1));
            
            raw_names = {};
            if isfield(m_struct, 'feature_names')
                raw_names = cellstr(m_struct.feature_names);
            elseif isfield(data, 'feature_names')
                raw_names = cellstr(data.feature_names);
            end
            
            n_f = length(imp_vec);
            mapped_names = cell(n_f, 1);
            
            for f_i = 1:n_f
                curr_name = '';
                if f_i <= length(raw_names) && ~isempty(raw_names{f_i})
                    curr_name = char(raw_names{f_i});
                end
                
                % Clinical Variable Name Mapping (Dynamic Lookup by Feature Name)
                if contains(m_key, 'Clin')
                    v_key = lower(strtrim(curr_name));
                    if isKey(clin_labels_map, v_key)
                        curr_name = clin_labels_map(v_key);
                    end
                end
                
                % PowCoh Band Renaming
                for b_idx = 1:5
                    str_num = num2str(b_idx); b_name = freq_bands{b_idx};
                    curr_name = strrep(curr_name, ['Power Band ' str_num], ['Power Band ' b_name]);
                    curr_name = strrep(curr_name, ['Coherence Pair ' str_num], ['Coherence Pair ' b_name]);
                    curr_name = strrep(curr_name, ['Band ' str_num], ['Band ' b_name]);
                    curr_name = strrep(curr_name, ['Pair ' str_num], ['Pair ' b_name]);
                end
                
                if contains(m_key, 'PowCoh') && ~startsWith(curr_name, 'Power') && ~startsWith(curr_name, 'Coh')
                    if f_i <= 5, curr_name = ['Power ' freq_bands{f_i}];
                    elseif f_i <= 10, curr_name = ['Coh ' freq_bands{f_i-5}]; end
                end
                
                % HCTSA Name Mapping
                if contains(m_key, 'HCTSA') && ~isempty(hctsa_names_table)
                    if startsWith(curr_name, 'HCTSA_Feat_')
                        idx_num = str2double(erase(curr_name, 'HCTSA_Feat_'));
                        if ~isnan(idx_num) && idx_num >= 1 && idx_num <= length(hctsa_names_table)
                            curr_name = hctsa_names_table{idx_num};
                        end
                    elseif ~isnan(str2double(curr_name))
                        idx_num = str2double(curr_name);
                        if idx_num >= 1 && idx_num <= length(hctsa_names_table)
                            curr_name = hctsa_names_table{idx_num};
                        end
                    end
                end
                
                mapped_names{f_i} = ternary(isempty(curr_name), sprintf('Feat_%d', f_i), curr_name);
            end
            
            T_m = table(string(mapped_names), imp_vec(:), 'VariableNames', {'Name', 'Score'});
            T_m = sortrows(T_m(~isnan(T_m.Score), :), 'Score', 'descend');
            
            % Set top-k dynamically per domain
            k_target = ternary(contains(m_key, 'HCTSA'), k_show_hctsa, k_show_default);
            k_show   = min(k_target, height(T_m));
            
            T_top  = T_m(1:k_show, :);
            scores_plot = T_top.Score(k_show:-1:1);
            labels_plot = cell(k_show, 1);
            
            % Max character limit (longer string display for HCTSA)
            max_char_len = ternary(contains(m_key, 'HCTSA'), 34, 26);
            
            for f_idx = 1:k_show
                str_n = strrep(T_top.Name(k_show - f_idx + 1), '_', ' ');
                if strlength(str_n) > max_char_len
                    labels_plot{f_idx} = char(extractBefore(str_n, max_char_len - 2) + "..");
                else
                    labels_plot{f_idx} = char(str_n);
                end
            end
            
            y_ticks_pos = 1:k_show;
            barh(y_ticks_pos, scores_plot, 0.65, 'FaceColor', m_color, 'EdgeColor', 'none');
            yticks(y_ticks_pos); yticklabels(labels_plot);
            xlabel('Importance Score');
            title(strrep(m_key, '_', ' '), 'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'none');
            grid on; 
            
            % Adjust font sizes dynamically based on panel density
            tick_font_size = ternary(k_show > 12, 7.5, 8.5);
            text_font_size = ternary(k_show > 12, 7.0, 7.5);
            set(gca, 'FontSize', tick_font_size, 'TickLabelInterpreter', 'none');
            
            max_sc = max(scores_plot);
            if max_sc > 0, xlim([0, max_sc * 1.30]); end
            
            % Adaptive decimal formatting based on magnitude of values
            for b = 1:k_show
                val = scores_plot(b);
                if max_sc < 0.001
                    str_val = sprintf('%.2e', val);
                elseif max_sc < 0.05
                    str_val = sprintf('%.4f', val);
                else
                    str_val = sprintf('%.3f', val);
                end
                
                text(val + (max_sc * 0.02), y_ticks_pos(b), str_val, ...
                    'VerticalAlignment', 'middle', 'FontSize', text_font_size, 'Color', [0.2 0.2 0.2]);
            end
        end
        
        % Save PNG & Vector PDF for Figure 2
        fig_feats_base = fullfile(target_path, sprintf('Classification_%s_LOSO_TopFeatures_SingleDomain', comp_name));
        saveas(fig_feats, [fig_feats_base '.png']);
        save_figure_as_pdf(fig_feats, [fig_feats_base '.pdf']);
        fprintf('Saved single-domain feature plot (PNG & PDF): %s\n', fig_feats_base);
    end
end

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function save_figure_as_pdf(fig_handle, pdf_filepath)
    try
        exportgraphics(fig_handle, pdf_filepath, 'ContentType', 'vector');
    catch
        set(fig_handle, 'PaperPositionMode', 'auto');
        print(fig_handle, pdf_filepath, '-dpdf', '-bestfit');
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