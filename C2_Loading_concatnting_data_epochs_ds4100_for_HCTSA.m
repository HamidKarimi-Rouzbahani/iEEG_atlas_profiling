clc;
clear all;
close all;

% Adding FieldTrip for consistency, though not strictly needed for HCTSA part
addpath(genpath('/home/uqhkarim/Normative_HCTSA/fieldtrip-20250523'))
ft_defaults % Initialize FieldTrip defaults

%% --- Configuration Constants ---
% Define directories for input (preprocessed) and output (HCTSA files)
PREPROCESSED_DATA_DIR = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset'; % Assuming preprocessed .mat files are in the current directory
HCTSA_OUTPUT_DIR = './';    % Directory to save HCTSA input file and HCTSA.mat
MAX_EPOCHS_TO_CHECK = 5;    % Maximum number of epochs expected per trial (5 based on typical Litt dataset split)

% --- Analysis Type (fixed to interictal as per request) ---
ictal_or_inter = 'interictal'; % Matches the naming convention of saved files

%% HCTSA Data Preparation Loop
% Initialize cell arrays to store all time series, labels, and keywords
all_timeSeriesData = {};
all_labels = {};
all_keywords = {};
series_counter = 0; % Global counter for all time series across all patients, trials, and epochs

fprintf('Starting HCTSA data preparation and concatenating all epochs...\n');

% The patient and trial loop must match exactly how your preprocessed files were generated
for Patient = [1:56]
    % Patient-specific data definitions (must match preprocessing script)
    if Patient==1
        Patient_initials='060'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==2
        Patient_initials='064'; modality='ecog'; inter_ictal_trials=[2];
    elseif Patient==3
        Patient_initials='065'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==4
        Patient_initials='070'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==5
        Patient_initials='074'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==6
        Patient_initials='075'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==7
        Patient_initials='080'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==8
        Patient_initials='082'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==9
        Patient_initials='086'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==10
        Patient_initials='087'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==11
        Patient_initials='088'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==12
        Patient_initials='089'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==13
        Patient_initials='094'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==14
        Patient_initials='097'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==15
        Patient_initials='105'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==16
        Patient_initials='106'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==17
        Patient_initials='107'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==18
        Patient_initials='111'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==19
        Patient_initials='112'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==20
        Patient_initials='114'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==21
        Patient_initials='116'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==22
        Patient_initials='117'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==23
        Patient_initials='123'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==24
        Patient_initials='126'; modality='ecog'; inter_ictal_trials=[1:2];
    elseif Patient==25
        Patient_initials='130'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==26
        Patient_initials='133'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==27
        Patient_initials='134'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==28
        Patient_initials='135'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==29
        Patient_initials='138'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==30
        Patient_initials='139'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==31
        Patient_initials='140'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==32
        Patient_initials='141'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==33
        Patient_initials='142'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==34
        Patient_initials='144'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==35
        Patient_initials='146'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==36
        Patient_initials='148'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==37
        Patient_initials='150'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==38
        Patient_initials='151'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==39
        Patient_initials='157'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==40
        Patient_initials='158'; modality='seeg'; inter_ictal_trials=[];
    elseif Patient==41
        Patient_initials='160'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==42
        Patient_initials='162'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==43
        Patient_initials='163'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==44
        Patient_initials='164'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==45
        Patient_initials='166'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==46
        Patient_initials='171'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==47
        Patient_initials='172'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==48
        Patient_initials='173'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==49
        Patient_initials='177'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==50
        Patient_initials='179'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==51
        Patient_initials='180'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==52
        Patient_initials='181'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==53
        Patient_initials='185'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==54
        Patient_initials='187'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==55
        Patient_initials='188'; modality='seeg'; inter_ictal_trials=[1:2];
    elseif Patient==56
        Patient_initials='190'; modality='seeg'; inter_ictal_trials=[1:2];
    end

    trials_to_process = inter_ictal_trials; % Always use the interictal trials list

    % Skip the current patient if no interictal trials are defined
    if isempty(trials_to_process)
        fprintf('Patient %d (HUP%s): No interictal trials defined or found for HCTSA. Skipping.\n', Patient, Patient_initials);
        continue;
    end
    
    % Loop through each specified interictal trial for the current patient
    for trial = trials_to_process
        fprintf('  Processing Patient %d (HUP%s), Trial %d...\n', Patient, Patient_initials, trial);
        
        epochs_found_count = 0; % Counter for epochs in this specific trial

        % --- MODIFIED: Loop through each expected epoch (1 to MAX_EPOCHS_TO_CHECK) ---
        for epoch_idx = 1:MAX_EPOCHS_TO_CHECK
            
            % Construct the filename for the epoched .mat file
            filename = sprintf('HUP%s_SEEG_InterIctal_Processed_Trial%02d_Epoch%02d_Bipolar_GM_only.mat', Patient_initials, trial, epoch_idx);
            filepath = fullfile(PREPROCESSED_DATA_DIR, filename);

            % Check if the epoched file exists. If not, break the epoch loop for this trial.
            if ~exist(filepath, 'file')
                if epochs_found_count == 0 % Only skip trial if the first epoch is missing
                     % fprintf('    WARNING: Epoch 1 not found for Trial %d (e.g. %s). Skipping trial.\n', trial, filename);
                else
                     fprintf('    Stop: Epoch %d not found for Trial %d. Finished epochs for this trial.\n', epoch_idx, trial);
                end
                break; % Exit epoch loop for current trial
            end
            
            epochs_found_count = epochs_found_count + 1;
            
            % Load the preprocessed data (only 'signal_epoch' and 'channels_id' needed)
            try
                loaded_data = load(filepath, 'signal_epoch', 'channels_id');
            catch ME
                fprintf(2, '    ERROR loading data from %s: %s\n', filepath, ME.message);
                continue; % Skip to next epoch
            end
            
            signal_epoch = loaded_data.signal_epoch;
            channels_id = loaded_data.channels_id;

            % Iterate through each channel in the loaded signal_epoch
            num_channels = size(signal_epoch, 1);
            for channel_idx = 1:num_channels
                series_counter = series_counter + 1;
                
                % Extract time series (ensure it's a row vector)
                tmp_data = signal_epoch(channel_idx, :); 
                all_timeSeriesData{series_counter} = tmp_data; 
                
                % --- Electrode Name Refinement ---
                original_channel_label = channels_id{channel_idx};
                
                % 1. Remove non-informative strings: 'EEG ' and '-Ref'
                cleaned_label = strrep(original_channel_label, 'EEG ', ''); % Remove 'EEG '
                cleaned_label = strrep(cleaned_label, '-Ref', ''); % Remove '-Ref'
                
                % 2. Add patient initials, trial_number, and EPOCH NUMBER to make unique
                unique_electrode_label = sprintf('HUP%s-%s-Trl%02d-Ep%02d', Patient_initials, cleaned_label, trial, epoch_idx);
                
                all_labels{series_counter} = unique_electrode_label; 
                
                % Create keywords for each time series. Include patient, trial, epoch, modality, unique electrode label.
                all_keywords{series_counter} = sprintf('%s,Patient%s,Trial%02d,Epoch%02d,Channel_%s', ...
                                                        modality, Patient_initials, trial, epoch_idx, unique_electrode_label);
            end % End of channel_idx loop
        end % End of epoch_idx loop

        if epochs_found_count > 0
            fprintf('  Added %d channels (%d epochs) from Patient %d (HUP%s), Trial %d. Total series collected: %d\n', ...
                num_channels * epochs_found_count, epochs_found_count, Patient, Patient_initials, trial, series_counter);
        end
        % Clear signal_epoch and loaded data for the next trial
        clear signal_epoch channels_id loaded_data; 
    end % End of trial loop
end % End of patient loop

fprintf('\nFinished preparing data for HCTSA. Total time series collected: %d\n', series_counter);

% Check if any time series were collected
if series_counter == 0
    fprintf('No time series were collected for HCTSA. Please check your preprocessed data files and paths.\n');
    return;
end

% Set the filename for the HCTSA input .mat file
hctsa_input_filename = fullfile(HCTSA_OUTPUT_DIR, 'HCTSA_AllPreprocessedSignals_epochs.mat');

% Save these variables out to the HCTSA input .mat file:
fprintf('Saving HCTSA input file to %s...\n', hctsa_input_filename);
% Rename variables to match HCTSA's expected input variable names
timeSeriesData = all_timeSeriesData;
labels = all_labels;
keywords = all_keywords;
save(hctsa_input_filename, 'timeSeriesData', 'labels', 'keywords', '-v7.3'); % Use -v7.3 for potentially large files

fprintf('HCTSA input file saved successfully.\n');