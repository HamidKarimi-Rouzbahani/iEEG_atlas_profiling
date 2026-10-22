clc; clear; close all;

% --- Setup Folders ---
output_folder = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/FeatureExtraction/';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

% --- Toolbox Path Configuration ---
addpath(genpath('/home/uqhkarim/Normative_HCTSA/Scripts/Spike-Ripple-Detector-Method-master'));

% --- Parameters ---
fs = 200; % Ensure this matches your data's sampling rate
Event_type = 1; % 1 = spike&ripple (HFO-spikes)

%% 2. Process DS4100 Dataset (HFO-Spike Extraction)
fprintf('\nLoading DS4100 data for HFO-Spike extraction...\n');
DS4100_data = load('/scratch/user/uqhkarim/Normative_HCTSA/Scripts/LittsScripts/HCTSA.mat');
DS4100_ann = load('/home/uqhkarim/Normative_HCTSA/Scripts/LittsScripts/ConcatenatedChannelMetadata_more_epochs.mat');

num_channels_ds = length(DS4100_data.TimeSeries.Data);
DS4100_Num_Spikes = zeros(num_channels_ds, 1);
DS4100_Spike_Start = cell(num_channels_ds, 1);

% Extract Patient IDs
patient_ids = {DS4100_ann.all_clinical_info_epoch.participant_id}';
unique_patients = unique(patient_ids);

fprintf('Extracting HFO-spike features for DS4100 by patient (%d patients, %d channels)...\n', ...
    length(unique_patients), num_channels_ds);

for p = 1:length(unique_patients)
    pid = unique_patients{p};
    p_idx = find(strcmp(patient_ids, pid));
    
    fprintf('  Processing DS4100 Patient: %s (%d channels)\n', string(pid), length(p_idx));
    
    for i = 1:length(p_idx)
        n = p_idx(i);
        sig = DS4100_data.TimeSeries.Data{n, 1};
        if isempty(sig), continue; end
        
        EEG = sig(:); % Column vector
        t = (0:length(EEG)-1)' / fs;
        
        [res, ~] = spike_ripple_detector(EEG, t, Event_type);
        
        if ~isempty(res.INPOS)
            raw_spikes = res.INPOS;
            if length(raw_spikes) > 1
                diffs = diff(raw_spikes);
                keep_idx = [true; diffs > 0.5];
                DS4100_Spike_Start{n} = raw_spikes(keep_idx);
            else
                DS4100_Spike_Start{n} = raw_spikes;
            end
            DS4100_Num_Spikes(n) = length(DS4100_Spike_Start{n});
        else
            DS4100_Spike_Start{n} = [];
            DS4100_Num_Spikes(n) = 0;
        end
    end
end

% --- Final Save ---
save_filename = fullfile(output_folder, 'DS4100_HFO_Spikes_Features.mat');
save(save_filename, 'DS4100_Num_Spikes', 'DS4100_Spike_Start', 'DS4100_ann', 'DS4100_data', '-v7.3');

fprintf('\nSuccess: HFO-Spike extraction complete for DS4100. Saved to %s\n', save_filename);