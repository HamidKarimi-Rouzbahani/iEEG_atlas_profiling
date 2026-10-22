clc; clear; close all;

% --- Setup Folders ---
output_folder = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/FeatureExtraction/';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

% --- Parameters ---
fs = 200; % Ensure this matches your data's sampling rate
freq_bands = [0.5 4; 4 8; 8 12; 12 30; 30 80]; 
band_names = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
num_bands = size(freq_bands, 1);

%% 1. Process Frauscher Dataset
fprintf('Loading Frauscher data...\n');
Frauscher_data = load('/home/uqhkarim/Normative_HCTSA/Scripts/Frauscher backup/HCTSA_W.mat');
Frauscher_annotations = load('/home/uqhkarim/Normative_HCTSA/Dataset/MNI_Frauscher/MatlabFile/MatlabFile.mat');

num_channels_f = length(Frauscher_data.TimeSeries.Data);
Frauscher_Power = zeros(num_channels_f, num_bands);
Frauscher_Coherence_Mean = zeros(num_channels_f, num_bands);

% Extract Patient IDs for Frauscher (1 to 110)
frauscher_pids = Frauscher_annotations.Patient; % Assuming this is a vector of same length as channels
unique_frauscher_pats = unique(frauscher_pids);

fprintf('Extracting features for Frauscher by patient (%d patients, %d channels)...\n', ...
    length(unique_frauscher_pats), num_channels_f);

for p = 1:length(unique_frauscher_pats)
    current_pid = unique_frauscher_pats(p);
    p_idx = find(frauscher_pids == current_pid);
    
    fprintf('  Processing Frauscher Patient: %d (%d channels)\n', current_pid, length(p_idx));
    
    % 1a. Power Calculation (Per Channel for this Patient)
    for i = 1:length(p_idx)
        n = p_idx(i);
        sig = Frauscher_data.TimeSeries.Data{n, 1};
        if isempty(sig), continue; end
        [pxx, f] = periodogram(sig, rectwin(length(sig)), length(sig), fs);
        for b = 1:num_bands
            idx = f >= freq_bands(b,1) & f <= freq_bands(b,2);
            Frauscher_Power(n, b) = mean(pxx(idx));
        end
    end
    
    % 1b. Coherence Calculation (Within Patient Only)
    if length(p_idx) > 1
        for b = 1:num_bands
            adj_matrix = eye(length(p_idx));
            for i = 1:length(p_idx)
                for j = i+1:length(p_idx)
                    sig1 = Frauscher_data.TimeSeries.Data{p_idx(i), 1};
                    sig2 = Frauscher_data.TimeSeries.Data{p_idx(j), 1};
                    if isempty(sig1) || isempty(sig2), continue; end
                    [mscoh, f_coh] = mscohere(sig1, sig2, hamming(fs), fs/2, fs, fs);
                    f_idx = f_coh >= freq_bands(b,1) & f_coh <= freq_bands(b,2);
                    adj_matrix(i,j) = mean(mscoh(f_idx));
                    adj_matrix(j,i) = adj_matrix(i,j);
                end
            end
            
            % Calculate mean coherence per contact for this patient
            for i = 1:length(p_idx)
                row_vals = adj_matrix(i, :);
                row_vals(i) = []; % Exclude self-loop
                Frauscher_Coherence_Mean(p_idx(i), b) = mean(row_vals, 'omitnan');
            end
        end
    end
end

save(fullfile(output_folder, 'Frauscher_Features_New_per_part.mat'), 'Frauscher_Power', 'Frauscher_Coherence_Mean', 'Frauscher_annotations', 'Frauscher_data', '-v7.3');

%% 2. Process DS4100 Dataset
fprintf('\nLoading DS4100 data...\n');
DS4100_data = load('/home/uqhkarim/Normative_HCTSA/Scripts/LittsScripts/HCTSA.mat');
DS4100_ann = load('/home/uqhkarim/Normative_HCTSA/Scripts/LittsScripts/ConcatenatedChannelMetadata_more_epochs.mat');

num_channels_ds = length(DS4100_data.TimeSeries.Data);
DS4100_Power = zeros(num_channels_ds, num_bands);
DS4100_Coherence_Mean = zeros(num_channels_ds, num_bands);

% Extract Patient IDs
patient_ids = {DS4100_ann.all_clinical_info_epoch.participant_id}';
unique_patients = unique(patient_ids);

fprintf('Extracting features for DS4100 by patient...\n');
for p = 1:length(unique_patients)
    pid = unique_patients{p};
    p_idx = find(strcmp(patient_ids, pid));
    
    % 2a. Power for patient channels
    for i = 1:length(p_idx)
        n = p_idx(i);
        sig = DS4100_data.TimeSeries.Data{n, 1};
        if isempty(sig), continue; end
        [pxx, f] = periodogram(sig, rectwin(length(sig)), length(sig), fs);
        for b = 1:num_bands
            idx = f >= freq_bands(b,1) & f <= freq_bands(b,2);
            DS4100_Power(n, b) = mean(pxx(idx));
        end
    end
    
    % 2b. Coherence for patient channels
    if length(p_idx) > 1
        for b = 1:num_bands
            adj_matrix = eye(length(p_idx));
            for i = 1:length(p_idx)
                for j = i+1:length(p_idx)
                    sig1 = DS4100_data.TimeSeries.Data{p_idx(i), 1};
                    sig2 = DS4100_data.TimeSeries.Data{p_idx(j), 1};
                    if isempty(sig1) || isempty(sig2), continue; end
                    
                    [mscoh, f_coh] = mscohere(sig1, sig2, hamming(fs), fs/2, fs, fs);
                    f_idx = f_coh >= freq_bands(b,1) & f_coh <= freq_bands(b,2);
                    adj_matrix(i,j) = mean(mscoh(f_idx));
                    adj_matrix(j,i) = adj_matrix(i,j);
                end
            end
            % Calculate mean coherence per contact
            for i = 1:length(p_idx)
                row_vals = adj_matrix(i, :);
                row_vals(i) = [];
                DS4100_Coherence_Mean(p_idx(i), b) = mean(row_vals, 'omitnan');
            end
        end
    end
end

% --- Final Save ---
save(fullfile(output_folder, 'DS4100_Features_New_per_part.mat'), ...
    'DS4100_Power', 'DS4100_Coherence_Mean', 'DS4100_ann', 'DS4100_data', 'band_names', '-v7.3');

fprintf('\nSuccess: Feature extraction complete for both datasets.\n');