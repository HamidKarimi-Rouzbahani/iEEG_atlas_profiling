clc;
close all;
clear all;

% Define output folder for saving results
output_folder = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';
% Ensure the output folder exists
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% --- Load the necessary data from the previously saved file ---
% input_data_filename = fullfile(output_folder, 'Data_AllContacts_BothDatasets.mat');
input_data_filename = fullfile(output_folder, 'Data_AllContacts_BothDatasets_epochs.mat');

if exist(input_data_filename, 'file')
    fprintf('Loading pre-processed data from %s...\n', input_data_filename);
    load(input_data_filename, ...
         'DS4100', ...
         'Frauscher', ...
         'DS4100_atlases_processed', ...
         'Frauscher_atlases_processed', ...
         'all_rois', ...
         'DS4100_annotations');
    fprintf('Data loaded successfully.\n');
else
    error('Pre-processed data file not found at %s. Please run the previous data preparation script first.', input_data_filename);
end


% --- Normalise Frauscher data to [0, 1] for every column separately ---
fprintf('\n--- Normalizing Frauscher data features to the [0, 1] range and saving min/max values ---\n');
num_frauscher_contacts = size(Frauscher, 1);
num_features = size(Frauscher, 2);

frauscher_min_vals = NaN(1, num_features);
frauscher_max_vals = NaN(1, num_features);
frauscher_normalized_between_0_1 = NaN(num_frauscher_contacts, num_features);

for f = 1:num_features
    current_col = Frauscher(:, f);
    min_val = min(current_col, [], 'omitnan');
    max_val = max(current_col, [], 'omitnan');
    
    frauscher_min_vals(f) = min_val; % Store min value for this feature
    frauscher_max_vals(f) = max_val; % Store max value for this feature
    
    range = max_val - min_val;
    
    if range == 0
        frauscher_normalized_between_0_1(:, f) = zeros(num_frauscher_contacts, 1);
    else
        frauscher_normalized_between_0_1(:, f) = (current_col - min_val) / range;
    end
end
Frauscher = frauscher_normalized_between_0_1;
fprintf('Frauscher data normalized successfully.\n');

% -------------------------------------------------------------------------
% --- REVISED: Normalise DS4100 data to Frauscher min/max, for each subject separately ---
% This new block iterates through each unique patient in DS4100 and applies the normalization.
fprintf('\n--- REVISING: Normalizing DS4100 data features for each patient separately, using Frauscher min/max values ---\n');

if isfield(DS4100_annotations, 'all_clinical_info_epoch') && ...
   ~isempty(DS4100_annotations.all_clinical_info_epoch) && ...
   isfield(DS4100_annotations.all_clinical_info_epoch(1), 'participant_id')
    
    patient_ids_per_channel = {DS4100_annotations.all_clinical_info_epoch.participant_id}';
    unique_patients = unique(patient_ids_per_channel);
    
    DS4100_normalized_per_patient = NaN(size(DS4100));

    for p_id = 1:length(unique_patients)
        current_patient_id = unique_patients{p_id};
        
        % Find all contacts belonging to the current patient
        patient_contact_indices = find(strcmp(patient_ids_per_channel, current_patient_id));
        
        if isempty(patient_contact_indices)
            continue;
        end
        
        patient_data = DS4100(patient_contact_indices, :);
        
        % Normalize each feature (column) for this patient's data
        normalized_patient_data = NaN(size(patient_data));
        for f = 1:size(patient_data, 2)
            current_col = patient_data(:, f);
            
            % Use the min/max values from the Frauscher dataset for this feature
            min_val = frauscher_min_vals(f);
            max_val = frauscher_max_vals(f);
            range = max_val - min_val;
            
            if range == 0
                normalized_patient_data(:, f) = (current_col == min_val);
                normalized_patient_data(~(current_col == min_val), f) = NaN;
            else
                normalized_patient_data(:, f) = (current_col - min_val) / range;
            end
        end
        
        % Place the normalized data back into the main matrix
        DS4100_normalized_per_patient(patient_contact_indices, :) = normalized_patient_data;
    end
    
    % Overwrite the old DS4100 data with the new normalized data
    DS4100 = DS4100_normalized_per_patient;
    fprintf('DS4100 data normalized successfully for %d patients.\n', length(unique_patients));

else
    error('Participant IDs not found in DS4100_annotations.all_clinical_info_epoch. Cannot perform patient-level normalization.');
end
% -------------------------------------------------------------------------

% Get dimensions (now using the newly normalized matrices)
num_ds4100_contacts = size(DS4100, 1);
num_features = size(DS4100, 2);

% Initialize matrix to store Z-scores for DS4100 contacts
DS4100_z_scores_from_normative = NaN(num_ds4100_contacts, num_features);

fprintf('\n--- Performing ROI-based Z-score normalization for DS4100 contacts against Frauscher normative data ---\n');

% Loop through each DS4100 contact
for i = 1:num_ds4100_contacts
    current_ds4100_contact_data = DS4100(i, :);
    current_ds4100_contact_roi = DS4100_atlases_processed{i};

    if strcmp(current_ds4100_contact_roi, 'No ROI') || isempty(current_ds4100_contact_roi)
        fprintf('  Contact %d (DS4100) has no assigned ROI. Skipping Z-score calculation for this contact.\n', i);
        continue; % Skip contacts without a defined ROI
    end

    % Find Frauscher contacts within the same ROI
    frauscher_roi_indices = find(strcmp(Frauscher_atlases_processed, current_ds4100_contact_roi));

    if isempty(frauscher_roi_indices)
        fprintf('  No Frauscher contacts found in ROI "%s" (for DS4100 contact %d). Skipping Z-score calculation for this contact.\n', current_ds4100_contact_roi, i);
        continue; % Cannot calculate Z-score if no normative data
    end

    % Get the normative data (Frauscher) for the current ROI
    normative_roi_data = Frauscher(frauscher_roi_indices, :);

    % Calculate mean and standard deviation for each feature in the normative ROI data
    mu_normative = mean(normative_roi_data, 1, 'omitnan'); % Mean across contacts for each feature
    sigma_normative = std(normative_roi_data, 0, 1, 'omitnan'); % Standard deviation across contacts for each feature

    zero_std_features = (sigma_normative == 0);

    current_contact_z_scores = NaN(1, num_features);

    for f = 1:num_features
        if zero_std_features(f)
            if current_ds4100_contact_data(f) == mu_normative(f)
                current_contact_z_scores(f) = 0;
            else
                current_contact_z_scores(f) = NaN;
                fprintf('Warning: Feature %d in ROI "%s" has zero standard deviation in normative data, and DS4100 contact %d value is different. Z-score set to NaN.\n', f, current_ds4100_contact_roi, i);
            end
        else
            current_contact_z_scores(f) = (current_ds4100_contact_data(f) - mu_normative(f)) / sigma_normative(f);
        end
    end

    DS4100_z_scores_from_normative(i, :) = current_contact_z_scores;

    if mod(i, 100) == 0 || i == num_ds4100_contacts
        fprintf('  Processed %d/%d DS4100 contacts...\n', i, num_ds4100_contacts);
    end
end

fprintf('\nROI-based Z-score calculation complete for DS4100 contacts.\n');

% --- Save Results ---
% output_zscore_filename = fullfile(output_folder, 'DS4100_ZScores_FromNormativeAtlas_Nrmlzd.mat');
output_zscore_filename = fullfile(output_folder, 'DS4100_ZScores_FromNormativeAtlas_Nrmlzd_epochs.mat');
fprintf('Saving DS4100 Z-scores (relative to Frauscher normative data) to %s...\n', output_zscore_filename);

save(output_zscore_filename, ...
     'DS4100_z_scores_from_normative', ...
     'DS4100', ...
     'DS4100_atlases_processed', ...
     'Frauscher_atlases_processed', ...
     'all_rois', ...
     'DS4100_annotations', ...
     '-v7.3');

fprintf('Z-score results saved successfully.\n');