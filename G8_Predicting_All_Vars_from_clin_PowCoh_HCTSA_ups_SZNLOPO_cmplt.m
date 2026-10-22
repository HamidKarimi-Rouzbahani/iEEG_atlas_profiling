clc; close all; clearvars;

% --- Start Parallel Pool if not already active ---
if isempty(gcp('nocreate'))
    parpool(); % Automatically starts pool with available CPU cores
end

% --- 1. Setup Paths & Output Directory ---
base_folder = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';
hctsa_file   = fullfile(base_folder, 'Data_AllContacts_BothDatasets_epochs.mat');
pow_coh_file = fullfile(base_folder, 'FeatureExtraction/Data_AllContacts_BothDatasets_epochs_PowCoh_per_part.mat');

output_dir = fullfile(base_folder, 'Classification_Results_LOSO_cmplt/');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

fprintf('Loading Raw Data...\n');
load(hctsa_file, 'DS4100', 'DS4100_annotations'); 
pc_data = load(pow_coh_file); 

% Extract HCTSA Feature Names / IDs if available
if isfield(DS4100_annotations, 'Operations') && isfield(DS4100_annotations.Operations, 'Name')
    hctsa_feat_names = string(DS4100_annotations.Operations.Name);
else
    hctsa_feat_names = arrayfun(@(x) sprintf('HCTSA_Feat_%d', x), 1:size(DS4100, 2), 'UniformOutput', false)';
end

% Channel Status Extraction
status_hctsa = lower(string(DS4100_annotations.all_channel_status_description));
status_pc    = lower(string(pc_data.DS4100_ann.all_channel_status_description));

% --- 2. Master Table & Alignment ---
clin_hctsa = DS4100_annotations.all_clinical_info_epoch;
patient_ids_h = erase(string({clin_hctsa.participant_id}), 'sub-')';
patient_ids_p = erase(string({pc_data.DS4100_ann.all_clinical_info_epoch.participant_id}), 'sub-')';

[~, first_idx] = unique(patient_ids_h, 'stable');
T_clin_master = struct2table(clin_hctsa(first_idx));
T_clin_master.participant_id = erase(string(T_clin_master.participant_id), 'sub-');
T_clin_master = sortrows(T_clin_master, 'participant_id'); 

% Standardize labels
T_clin_master.outcome = strtrim(upper(string(T_clin_master.outcome)));
T_clin_master.duration = double(T_clin_master.age) - double(T_clin_master.age_onset);

% Clinical Predictors (Excluding 'engel' to avoid data leakage)
clin_preds = {'lesion_status', 'sex', 'target', 'therapy', 'hand', 'implant', 'age', 'age_onset', 'duration'};
for i = 1:length(clin_preds)
    vn = clin_preds{i};
    if ~ismember(vn, {'age', 'age_onset', 'duration'})
        T_clin_master.(vn) = categorical(string(T_clin_master.(vn))); 
    else
        T_clin_master.(vn) = double(T_clin_master.(vn));
    end
end

% Raw Features
X_hctsa_raw_full = DS4100;
X_pc_raw_full    = [pc_data.DS4100_Power, pc_data.DS4100_Coherence_Mean];

% Define PowCoh Feature Names (Mapped to Delta, Theta, Alpha, Beta, Gamma)
num_pow = size(pc_data.DS4100_Power, 2);
num_coh = size(pc_data.DS4100_Coherence_Mean, 2);

standard_pow_bands = {'Power_Delta', 'Power_Theta', 'Power_Alpha', 'Power_Beta', 'Power_Gamma'};
if num_pow == 5
    pow_names = standard_pow_bands';
else
    pow_names = arrayfun(@(x) sprintf('Power_Band_%d', x), 1:num_pow, 'UniformOutput', false)';
end

standard_coh_bands = {'Coh_Delta', 'Coh_Theta', 'Coh_Alpha', 'Coh_Beta', 'Coh_Gamma'};
if num_coh == 5
    coh_names = standard_coh_bands';
else
    coh_names = arrayfun(@(x) sprintf('Coherence_Pair_%d', x), 1:num_coh, 'UniformOutput', false)';
end

powcoh_feat_names = [pow_names; coh_names];

% Single mode for high efficiency
current_mode = 'SOZ_minus_NonSOZ';

% Target Variables to Run: ONLY outcome, lesion_status, target
target_variables = {'outcome', 'lesion_status', 'target'};

% Analysis Configuration
num_repeats = 100; 
min_samples_per_cat = 2; 

% --- 3. Synchronized Feature Aggregation ---
fprintf('\nAggregating Features for Mode: %s...\n', current_mode);
unique_pats = T_clin_master.participant_id;

X_h_avg_master = nan(height(T_clin_master), size(X_hctsa_raw_full, 2));
X_p_avg_master = nan(height(T_clin_master), size(X_pc_raw_full, 2));

for p = 1:height(T_clin_master)
    curr_pat = unique_pats(p);
    
    % HCTSA Masking
    mask_h = (patient_ids_h == curr_pat);
    is_soz_h = (status_hctsa == "soz" | status_hctsa == "resect,soz");
    is_non_h = ~is_soz_h;
    
    Xs_h = X_hctsa_raw_full(mask_h & is_soz_h', :);
    Xn_h = X_hctsa_raw_full(mask_h & is_non_h', :);
    if ~isempty(Xs_h) && ~isempty(Xn_h)
        X_h_avg_master(p, :) = nanmean(Xs_h, 1) - nanmean(Xn_h, 1);
    end
    
    % PowCoh Masking
    mask_p = (patient_ids_p == curr_pat);
    is_soz_p = (status_pc == "soz" | status_pc == "resect,soz");
    is_non_p = ~is_soz_p;
    
    Xs_p = X_pc_raw_full(mask_p & is_soz_p', :);
    Xn_p = X_pc_raw_full(mask_p & is_non_p', :);
    if ~isempty(Xs_p) && ~isempty(Xn_p)
        X_p_avg_master(p, :) = nanmean(Xs_p, 1) - nanmean(Xn_p, 1);
    end
end

% Impute missing features with zero
X_h_avg_master(isnan(X_h_avg_master)) = 0; 
X_p_avg_master(isnan(X_p_avg_master)) = 0;

num_hctsa_feats  = size(X_h_avg_master, 2);
num_powcoh_feats = size(X_p_avg_master, 2);

% --- 4. Iterate Over Target Variables ---
for tv = 1:length(target_variables)
    target_var = target_variables{tv};
    target_out_dir = fullfile(output_dir, target_var);
    
    % Resume check
    if exist(target_out_dir, 'dir')
        existing_mat_files = dir(fullfile(target_out_dir, 'Classification_*.mat'));
        if ~isempty(existing_mat_files)
            fprintf('\n[SKIP] Target [%s]: Results already exist in %s.\n', target_var, target_out_dir);
            continue;
        end
    else
        mkdir(target_out_dir);
    end
    
    raw_labels = string(T_clin_master.(target_var));
    valid_mask = ~ismember(upper(raw_labels), ["", "NAN", "N/A", "UNKNOWN", "<UNDEFINED>", "<MISSING>"]) & ~ismissing(raw_labels);
    labels_sub = raw_labels(valid_mask);
    
    % Filter categories with < min_samples_per_cat
    [cat_counts, cat_names] = groupcounts(labels_sub);
    valid_categories = cat_names(cat_counts >= min_samples_per_cat);
    sample_keep_mask = ismember(labels_sub, valid_categories);
    
    T_clin_sub = T_clin_master(valid_mask, :);
    T_clin_sub = T_clin_sub(sample_keep_mask, :);
    
    X_h_sub = X_h_avg_master(valid_mask, :);
    X_h_sub = X_h_sub(sample_keep_mask, :);
    
    X_p_sub = X_p_avg_master(valid_mask, :);
    X_p_sub = X_p_sub(sample_keep_mask, :);
    
    labels_sub = labels_sub(sample_keep_mask);
    unique_cats = unique(labels_sub);
    
    if length(unique_cats) < 2
        fprintf('[SKIP] Target [%s]: Fewer than 2 categories have >= %d samples.\n', target_var, min_samples_per_cat);
        continue;
    end
    
    if length(unique_cats) == 2
        pair_list = {unique_cats'};
    else
        pairs = nchoosek(unique_cats, 2);
        pair_list = cell(size(pairs, 1), 1);
        for k = 1:size(pairs, 1)
            pair_list{k} = pairs(k, :);
        end
    end
    
    fprintf('\n=======================================================\n');
    fprintf('  RUNNING TARGET: %s (%d Pairwise Comparison(s))\n', target_var, length(pair_list));
    fprintf('=======================================================\n');
    
    % Clinical Predictors logic according to target
    if strcmp(target_var, 'outcome')
        curr_preds = setdiff(clin_preds, 'outcome');
    else
        % For lesion_status and target: drop outcome, engel, therapy
        curr_preds = setdiff(clin_preds, {target_var, 'outcome', 'therapy', 'engel'});
    end

    model_names = {'Clin_Only', 'HCTSA_Only', 'PowCoh_Only', 'HCTSA_PowCoh', 'Clin_HCTSA', 'Clin_PowCoh', 'All_Features'};
    
    macro_repeat_accs = cell(length(model_names), 1);
    macro_aucs        = cell(length(model_names), 1);
    macro_rocs        = cell(length(model_names), 1);
    macro_importances = cell(length(model_names), 1);
    
    for m = 1:length(model_names)
        macro_repeat_accs{m} = zeros(num_repeats, 0);
        macro_aucs{m}        = [];
        macro_rocs{m}        = [];
        macro_importances{m} = [];
    end

    pair_results = struct();

    for pair_idx = 1:length(pair_list)
        curr_pair = pair_list{pair_idx};
        pair_str = sprintf('%s_vs_%s', matlab.lang.makeValidName(curr_pair(1)), matlab.lang.makeValidName(curr_pair(2)));
        
        pair_mask = ismember(labels_sub, curr_pair);
        
        T_clin = T_clin_sub(pair_mask, :);
        
        % Clean up unused categorical levels for this pairwise subset
        for col = 1:width(T_clin)
            if iscategorical(T_clin{:, col})
                T_clin{:, col} = removecats(T_clin{:, col});
            end
        end
        
        X_h_avg = X_h_sub(pair_mask, :);
        X_p_avg = X_p_sub(pair_mask, :);
        
        % Strict string conversion to ensure consistent class ordering
        Y_str = string(labels_sub(pair_mask));
        unique_labels = unique(Y_str, 'stable');
        pos_class_str = unique_labels{1}; % First label set as explicit positive class
        
        empty_T = table('Size', [height(T_clin) 0]);
        
        configs = {
            'Clin_Only',    T_clin(:, curr_preds);
            'HCTSA_Only',   empty_T; 
            'PowCoh_Only',  empty_T;
            'HCTSA_PowCoh', empty_T;
            'Clin_HCTSA',   T_clin(:, curr_preds); 
            'Clin_PowCoh',  T_clin(:, curr_preds);
            'All_Features', T_clin(:, curr_preds)
        };

        % Set up Leave-One-Subject-Out CV partition
        num_subjects = height(T_clin);
        cv = cvpartition(num_subjects, 'LeaveOut');
        num_folds = cv.NumTestSets; % Equal to total number of subjects (N)
        
        for c = 1:size(configs, 1)
            m_name = configs{c,1};
            T_cfg  = configs{c,2};
            
            repeat_accs     = zeros(num_repeats, 1);
            repeat_scores   = cell(num_repeats, 1);
            repeat_true_bin = cell(num_repeats, 1);
            
            has_hctsa  = contains(m_name, 'HCTSA')  || strcmp(m_name, 'All_Features');
            has_powcoh = contains(m_name, 'PowCoh') || strcmp(m_name, 'All_Features');
            has_clin   = contains(m_name, 'Clin')   || strcmp(m_name, 'All_Features');
            
            % Build composite feature names for model logging
            model_feat_names = string.empty;
            if has_clin
                model_feat_names = [model_feat_names; string(curr_preds')];
            end
            if has_hctsa
                model_feat_names = [model_feat_names; hctsa_feat_names(:)];
            end
            if has_powcoh
                model_feat_names = [model_feat_names; powcoh_feat_names(:)];
            end
            
            total_model_feats = length(model_feat_names);
            repeat_imp = cell(num_repeats, 1);
            
            % --- REPEAT LOOP ---
            for r = 1:num_repeats
                
                fold_accs      = zeros(num_folds, 1);
                fold_scores    = zeros(num_folds, 1);
                fold_true_bin  = zeros(num_folds, 1);
                fold_feat_imps = zeros(num_folds, total_model_feats);
                
                parfor f = 1:num_folds
                    trIdx = training(cv, f); teIdx = test(cv, f);
                    T_tr = T_cfg(trIdx, :); T_te = T_cfg(teIdx, :);
                    coeff_h = [];
                    
                    % Add HCTSA
                    if has_hctsa
                        mu = mean(X_h_avg(trIdx,:), 1); sd = std(X_h_avg(trIdx,:), 0, 1); sd(sd==0)=1;
                        tr_z = (X_h_avg(trIdx,:) - mu)./sd; tr_z(isnan(tr_z)) = 0;
                        n_comp = min(10, max(1, size(tr_z, 1) - 1));
                        [coeff_h, score] = pca(tr_z, 'NumComponents', n_comp);
                        te_pca = ((X_h_avg(teIdx,:) - mu)./sd) * coeff_h;
                        h_names = arrayfun(@(x) sprintf('H_PC%d',x), 1:n_comp, 'Un', 0);
                        T_tr = [T_tr, array2table(score, 'VariableNames', h_names)];
                        T_te = [T_te, array2table(te_pca, 'VariableNames', h_names)];
                    end
                    
                    % Add PowCoh
                    if has_powcoh
                        mu = mean(X_p_avg(trIdx,:), 1); sd = std(X_p_avg(trIdx,:), 0, 1); sd(sd==0)=1;
                        tr_z = (X_p_avg(trIdx,:) - mu)./sd; tr_z(isnan(tr_z)) = 0;
                        te_z = (X_p_avg(teIdx,:) - mu)./sd; te_z(isnan(te_z)) = 0;
                        p_names = arrayfun(@(x) sprintf('PC_%d',x), 1:size(X_p_avg,2), 'Un', 0);
                        T_tr = [T_tr, array2table(tr_z, 'VariableNames', p_names)];
                        T_te = [T_te, array2table(te_z, 'VariableNames', p_names)];
                    end
                    
                    % Class Balancing (Upsampling on training fold)
                    Y_tr_str = Y_str(trIdx); 
                    uL = unique(Y_tr_str);
                    counts = arrayfun(@(x) sum(Y_tr_str == x), uL);
                    n_max = max(counts);
                    bal_idx = [];
                    for k = 1:length(uL)
                        idx = find(Y_tr_str == uL(k));
                        if ~isempty(idx)
                            bal_idx = [bal_idx; idx(randi(length(idx), n_max, 1))];
                        end
                    end
                    
                    % Train TreeBagger with OOB Feature Importance enabled
                    Mdl = TreeBagger(150, T_tr(bal_idx,:), Y_tr_str(bal_idx), ...
                        'Method', 'classification', 'OOBPredictorImportance', 'on');
                    
                    [f_preds, f_scores] = predict(Mdl, T_te);
                    
                    % --- Feature Importances ---
                    raw_imp = Mdl.OOBPermutedPredictorDeltaError;
                    var_names = string(Mdl.PredictorNames);
                    
                    fold_imp = zeros(1, total_model_feats);
                    curr_offset = 0;
                    
                    % 1) Clinical features
                    if has_clin
                        n_clin = length(curr_preds);
                        for cl_i = 1:n_clin
                            c_name = string(curr_preds{cl_i});
                            match_idx = find(var_names == c_name);
                            if ~isempty(match_idx)
                                fold_imp(cl_i) = raw_imp(match_idx);
                            end
                        end
                        curr_offset = curr_offset + n_clin;
                    end
                    
                    % 2) HCTSA features (PCA Re-projection Fix applied)
                    if has_hctsa
                        if ~isempty(coeff_h)
                            h_pc_mask = startsWith(var_names, "H_PC");
                            if any(h_pc_mask)
                                pc_imp = raw_imp(h_pc_mask);
                                % Normalize PCA loadings weight so projected feature importance 
                                % accurately reflects the total PCA signal contribution relative to uncompressed features
                                load_weights = coeff_h.^2;
                                norm_weights = load_weights ./ max(eps, sum(load_weights, 1));
                                fold_imp(curr_offset + (1:num_hctsa_feats)) = pc_imp * norm_weights'; 
                            end
                        end
                        curr_offset = curr_offset + num_hctsa_feats;
                    end
                    
                    % 3) PowCoh features
                    if has_powcoh
                        pc_mask = startsWith(var_names, "PC_");
                        if any(pc_mask)
                            fold_imp(curr_offset + (1:num_powcoh_feats)) = raw_imp(pc_mask);
                        end
                    end
                    
                    fold_feat_imps(f, :) = fold_imp;
                    
                    % --- STRICT PROBABILITY MATCHING FIX ---
                    class_names_str = string(Mdl.ClassNames);
                    colIdx = find(class_names_str == pos_class_str);
                    
                    if ~isempty(colIdx)
                        fold_scores(f) = double(f_scores(1, colIdx));
                    else
                        fold_scores(f) = 0;
                    end
                    
                    fold_true_bin(f) = double(Y_str(teIdx) == pos_class_str);
                    fold_accs(f)     = double(string(f_preds{1}) == Y_str(teIdx));
                end
                
                repeat_accs(r) = mean(fold_accs);
                repeat_scores{r} = fold_scores;
                repeat_true_bin{r} = fold_true_bin;
                repeat_imp{r} = mean(fold_feat_imps, 1);
            end
            
            % Consolidate outputs across repeats
            all_scores   = vertcat(repeat_scores{:});
            all_true_bin = vertcat(repeat_true_bin{:});
            
            % Compute ROC with explicit positive class = 1
            [Xroc, Yroc, ~, AUC] = perfcurve(all_true_bin, all_scores, 1);
            
            % Upstream Check: If scores flip, reflect coordinates
            if AUC < 0.5
                [Xroc, Yroc, ~, AUC] = perfcurve(all_true_bin, -all_scores, 1);
            end
            
            [unique_Xroc, unique_idx] = unique(Xroc, 'stable');
            unique_Yroc = Yroc(unique_idx);
            
            fpr_grid = linspace(0, 1, 100)';
            tpr_interp = interp1(unique_Xroc, unique_Yroc, fpr_grid, 'linear', 'extrap');
            tpr_interp = max(0, min(1, tpr_interp)); % Clamp bounds to [0, 1]
            
            avg_pair_imp = mean(vertcat(repeat_imp{:}), 1);
            
            % Store individual pair metrics
            pair_results.(pair_str).(m_name).mean_acc = mean(repeat_accs * 100);
            pair_results.(pair_str).(m_name).std_acc  = std(repeat_accs * 100);
            pair_results.(pair_str).(m_name).var_acc  = var(repeat_accs * 100);
            pair_results.(pair_str).(m_name).raw_repeat_accs = repeat_accs * 100;
            pair_results.(pair_str).(m_name).auc = AUC;
            pair_results.(pair_str).(m_name).roc = [fpr_grid, tpr_interp];
            pair_results.(pair_str).(m_name).feature_importance = avg_pair_imp;
            pair_results.(pair_str).(m_name).feature_names = model_feat_names;
            
            % Accumulate macro-averages
            macro_repeat_accs{c} = [macro_repeat_accs{c}, repeat_accs * 100];
            macro_aucs{c}        = [macro_aucs{c}, AUC];
            macro_rocs{c}        = [macro_rocs{c}, tpr_interp];
            macro_importances{c} = [macro_importances{c}; avg_pair_imp];
        end
    end
    
    % --- 5. Consolidate Macro-Averaged Results Across Pairs ---
    cv_results = struct();
    p_fields = fieldnames(pair_results);
    
    for c = 1:length(model_names)
        m_name = model_names{c};
        
        avg_repeats = mean(macro_repeat_accs{c}, 2);
        cv_results.(m_name).mean_acc = mean(avg_repeats);
        cv_results.(m_name).std_acc  = std(avg_repeats);
        cv_results.(m_name).var_acc  = var(avg_repeats);
        cv_results.(m_name).raw_repeat_accs = avg_repeats;
        
        cv_results.(m_name).auc = mean(macro_aucs{c});
        
        fpr_grid = linspace(0, 1, 100)';
        mean_tpr = mean(macro_rocs{c}, 2);
        cv_results.(m_name).roc = [fpr_grid, mean_tpr];
        
        % Macro feature importances
        cv_results.(m_name).feature_importance = mean(macro_importances{c}, 1);
        
        % Retain exact model feature names
        cv_results.(m_name).feature_names = pair_results.(p_fields{1}).(m_name).feature_names;
    end
    
    % --- 6. Plotting & MAT Export ---
    comp_name = target_var;
    fig = figure('Color','w','Position',[50 50 1400 450], 'Visible', 'off');
    names = fieldnames(cv_results);
    names_clean = strrep(names, '_', ' ');

    % Plot 1: Mean Accuracy (%)
    subplot(1,3,1); 
    bar(cellfun(@(n) cv_results.(n).mean_acc, names), 'FaceColor', [0.2 0.4 0.6]); hold on;
    errorbar(1:length(names), cellfun(@(n) cv_results.(n).mean_acc, names), ...
        cellfun(@(n) cv_results.(n).std_acc, names), 'k.', 'LineWidth', 1.2);
    title(sprintf('Mean Accuracy (%%)\n[%s]', strrep(comp_name,'_',' ')), 'Interpreter', 'none'); 
    grid on; set(gca, 'XTick', 1:length(names), 'XTickLabel', names_clean, 'XTickLabelRotation', 30);
    ylim([0 100]); yline(50, '--r', 'Chance');

    % Plot 2: Accuracy Variance
    subplot(1,3,2); 
    bar(cellfun(@(n) cv_results.(n).var_acc, names), 'FaceColor', [0.8 0.4 0.2]);
    title(sprintf('Variance across Runs\n[%s]', strrep(comp_name,'_',' ')), 'Interpreter', 'none'); 
    ylabel('Variance (%^2)'); grid on; set(gca, 'XTick', 1:length(names), 'XTickLabel', names_clean, 'XTickLabelRotation', 30);

    % Plot 3: Aggregated ROC Curves
    subplot(1,3,3); 
    hold on; colors = lines(length(names));
    for i = 1:length(names)
        plot(cv_results.(names{i}).roc(:,1), cv_results.(names{i}).roc(:,2), 'LineWidth', 2, 'Color', colors(i,:));
    end
    plot([0 1], [0 1], 'k--'); 
    title(sprintf('Macro-Averaged ROC\n[%s]', strrep(comp_name,'_',' ')), 'Interpreter', 'none'); 
    legend(names_clean, 'Location', 'southeast'); grid on;

    % Save MAT file
    mat_filename = fullfile(target_out_dir, sprintf('Classification_%s_LOSO_Results_complete.mat', comp_name));
    save(mat_filename, 'cv_results', 'pair_results', 'comp_name', 'num_folds', 'num_repeats');

    % Save Figure
    fig_filename = fullfile(target_out_dir, sprintf('Classification_%s_LOSO.png', comp_name));
    saveas(fig, fig_filename); close(fig);

    fprintf('Completed target [%s]. LOSO results saved to MAT file.\n', target_var);
end

fprintf('\nPipeline processing complete!\n');