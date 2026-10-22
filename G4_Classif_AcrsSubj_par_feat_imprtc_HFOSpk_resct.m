clc; clear; close all;

% ==========================================
% --- SECTION 0: PATH CONFIGURATIONS ---
% ==========================================
base_dir = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';

% Define precise input file and output directory targets
input_file = fullfile(base_dir, 'FeatureExtraction/DS4100_HFO_Spikes_Features.mat');
output_dir = fullfile(base_dir, 'AcrossSubjects/HFO_Spike_Results/');

% Ensure output directory exists
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Open parallel processing pool if not already active
if isempty(gcp('nocreate')), fprintf('Starting parallel pool...\n'); parpool; end

fprintf('\n=======================================================\n');
fprintf('STARTING UNIFIED PIPELINE FOR FEATURE MODE: HFO_Spikes\n');
fprintf('=======================================================\n');

% ==========================================
% --- SECTION 1: UNIFIED DATA INGESTION ---
% ==========================================
load(input_file, 'DS4100_Num_Spikes', 'DS4100_ann');

% Extract the 1-dimensional feature matrix (Channel count x 1)
X_raw = double(DS4100_Num_Spikes); 
if size(X_raw, 2) ~= 1
    X_raw = X_raw'; % Ensure column vector representation
end
feature_names_kept = {'HFO_Spike_Count'};

% Extract channel status and clinical metadata
all_channel_status = DS4100_ann.all_channel_status_description;
if isfield(DS4100_ann, 'all_clinical_info_epoch')
    clin_info_raw = DS4100_ann.all_clinical_info_epoch;
else
    clin_info_raw = DS4100_ann.all_clinical_info;
end

% Unify patient identifiers and build valid classification masks
temp_pids = {clin_info_raw.participant_id}';
patient_ids_raw = erase(string(temp_pids), 'sub-');

Y_raw = NaN(size(X_raw, 1), 1);
Y_raw(strcmp(all_channel_status, 'resect,soz')) = 1;
Y_raw(strcmp(all_channel_status, 'soz') | strcmp(all_channel_status, 'n/a') | strcmp(all_channel_status, 'resect')) = 0;

valid_idx = ~isnan(Y_raw);
X_full = X_raw(valid_idx, :);
Y_full = Y_raw(valid_idx);
patient_ids_full = patient_ids_raw(valid_idx);
clin_info_full = clin_info_raw(valid_idx);

% Unify sparse structure definitions into a regularized table representation
all_field_names = {};
for i = 1:length(clin_info_full)
    all_field_names = union(all_field_names, fieldnames(clin_info_full(i)));
end
temp_struct = clin_info_full;
for i = 1:length(temp_struct)
    missing_fields = setdiff(all_field_names, fieldnames(temp_struct(i)));
    for f = 1:length(missing_fields)
        if ismember(missing_fields{f}, {'age', 'age_onset', 'outcome', 'engel'})
            temp_struct(i).(missing_fields{f}) = NaN;
        else
            temp_struct(i).(missing_fields{f}) = '';
        end
    end
end
all_test_clinical_info = struct2table(temp_struct, 'AsArray', true);

% ==========================================
% --- SECTION 2: CROSS-VALIDATION SETUP ---
% ==========================================
n_times = 10; % 10 repetitions to eliminate initialization variance
rng('default'); 

unique_patients = unique(patient_ids_full);
num_patients = length(unique_patients);

iter_overall_auc = zeros(n_times, 1);
iter_all_test_scores = cell(n_times, 1);
iter_feat_importance_matrices = cell(n_times, 1);

% ==========================================
% --- SECTION 3: CORE MODEL TRAINING LOOP ---
% ==========================================
for t = 1:n_times
    fprintf('>>> Repetition %d/%d <<<\n', t, n_times);
    
    fold_auc = zeros(num_patients, 1);
    fold_scores = cell(num_patients, 1);
    fold_feat_imp = cell(num_patients, 1);
    
    parfor fold = 1:num_patients
        test_patient = unique_patients(fold);
        test_mask = (patient_ids_full == test_patient);
        train_mask = ~test_mask;
        
        X_train = X_full(train_mask, :);
        Y_train = Y_full(train_mask);
        X_test  = X_full(test_mask, :);
        Y_test  = Y_full(test_mask);
        
        % Dynamic feature tracking mask (1-dimensional)
        features_to_keep = true(1, size(X_train, 2));
        
        % --- Outlier Handling and Mean Imputation ---
        outls = isoutlier(X_train, 'percentiles', [0.03 99.7]);
        X_train(outls) = NaN;
        
        for j = 1:size(X_train, 2)
            m_val = nanmean(X_train(:, j));
            if isnan(m_val), m_val = 0; end % Fallback safety
            X_train(isnan(X_train(:, j)), j) = m_val;
            X_test(isnan(X_test(:, j)), j)   = m_val;
        end
        
        % --- Unified Dynamic Class Upsampling Mechanics ---
        pos_idx = find(Y_train == 1); 
        neg_idx = find(Y_train == 0);
        
        if isempty(pos_idx) || isempty(neg_idx)
            X_train_bal = X_train; 
            Y_train_bal = Y_train;
        else
            num_pos = length(pos_idx);
            num_neg = length(neg_idx);
            
            if num_pos < num_neg
                shortfall = num_neg - num_pos;
                resample_idx = randsample(pos_idx, shortfall, true);
                X_train_bal = [X_train; X_train(resample_idx, :)];
                Y_train_bal = [Y_train; Y_train(resample_idx)];
            elseif num_neg < num_pos
                shortfall = num_pos - num_neg;
                resample_idx = randsample(neg_idx, shortfall, true);
                X_train_bal = [X_train; X_train(resample_idx, :)];
                Y_train_bal = [Y_train; Y_train(resample_idx)];
            else
                X_train_bal = X_train; 
                Y_train_bal = Y_train;
            end
            
            % Scramble training rows to break up sampling sequential bias
            shuf_idx = randperm(length(Y_train_bal));
            X_train_bal = X_train_bal(shuf_idx, :); 
            Y_train_bal = Y_train_bal(shuf_idx);
        end
        
        % --- Z-Score Normalization ---
        mu = mean(X_train_bal, 1); 
        sigma = std(X_train_bal, 0, 1); 
        sigma(sigma == 0) = 1; % Shield against divide-by-zero errors
        
        X_train_sc = (X_train_bal - mu) ./ sigma;
        X_test_sc  = (X_test - mu) ./ sigma;
        
        % --- Random Forest Classification Ensemble ---
        Mdl = TreeBagger(100, X_train_sc, Y_train_bal, 'Method', 'classification', ...
            'OOBPredictorImportance', 'on', 'MinLeafSize', 5);
        
        % Extract and map feature importance scores back to global index array
        imp_full = NaN(1, size(X_full, 2));
        imp_full(features_to_keep) = Mdl.OOBPermutedPredictorDeltaError;
        fold_feat_imp{fold} = imp_full;
        
        [~, scores] = predict(Mdl, X_test_sc);
        pos_class_col = find(str2double(Mdl.ClassNames) == 1);
        if ~isempty(pos_class_col)
            fold_scores{fold} = scores(:, pos_class_col);
        else
            fold_scores{fold} = zeros(size(X_test, 1), 1);
        end
        
        if length(unique(Y_test)) >= 2
            [~,~,~,auc_val] = perfcurve(Y_test, fold_scores{fold}, 1);
            fold_auc(fold) = auc_val;
        else
            fold_auc(fold) = NaN;
        end
    end
    
    iter_overall_auc(t) = nanmean(fold_auc);
    
    % --- Rank-normalize concatenated LOPO scores per repetition ---
    r_scores = vertcat(fold_scores{:});
    iter_all_test_scores{t} = tiedrank(r_scores) / length(r_scores);
    
    iter_feat_importance_matrices{t} = vertcat(fold_feat_imp{:});
end

% ==========================================
% --- SECTION 4: UNIFIED RESULT MATRIX OUTPUT ---
% ==========================================
final_mean_auc = mean(iter_overall_auc);
overall_auc_scores = iter_overall_auc;

% Compute average vectors/matrices across iterations
all_test_scores = mean(cat(3, iter_all_test_scores{:}), 3);
all_importances_matrix = mean(cat(3, iter_feat_importance_matrices{:}), 3);

all_test_labels = Y_full; 
all_test_patient_ids = patient_ids_full;

% Generate structural data file for downstream execution
output_filename = fullfile(output_dir, 'AcrossSubjects_Final_Results_HFO_Spikes.mat');
save(output_filename, 'overall_auc_scores', 'final_mean_auc', ...
     'all_test_labels', 'all_test_scores', 'all_test_patient_ids', ...
     'all_test_clinical_info', 'all_importances_matrix', 'feature_names_kept', '-v7.3');
 
fprintf('SUCCESS: Finished processing and saved results for HFO_Spikes feature mode to %s\n', output_filename);