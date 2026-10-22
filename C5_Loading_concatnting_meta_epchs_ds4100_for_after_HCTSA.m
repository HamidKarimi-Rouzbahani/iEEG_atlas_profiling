clc;
clear all;
close all;

% Adding FieldTrip for consistency, though not strictly needed for this modified script
addpath(genpath('/home/uqhkarim/Normative_HCTSA/fieldtrip-20250523'))
ft_defaults % Initialize FieldTrip defaults

%% --- Configuration Constants ---
PREPROCESSED_DATA_DIR = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset'; % Directory where preprocessed .mat files are located
OUTPUT_DIR = './'; % Directory to save the concatenated metadata file
MAX_EPOCHS_TO_CHECK = 5; % Maximum number of epochs expected per trial

% --- Analysis Type (fixed to interictal as per request) ---
ictal_or_inter = 'interictal'; % Matches the naming convention of saved files

%% Data Concatenation Loop
% Initialize variables to store all concatenated metadata
% For struct arrays, initialize as empty struct arrays to ensure correct concatenation
all_channel_info_epoch = struct([]); 
all_electrode_loc_info_epoch = struct([]); 
all_clinical_info_epoch = struct([]); % Initialize for clinical info
% For cell arrays, initialize as empty cell arrays
all_channel_status_description = {}; 
all_channels_id = {}; 

fprintf('Starting metadata concatenation across all epochs...\n');

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
        Patient_initials='158'; modality='seeg'; inter_ictal_trials=[]; % Patient 158 seems to have no interictal trials listed
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

    trials_to_process = inter_ictal_trials;    

    % --- Skip patient if no interictal trials ---
    if isempty(trials_to_process)
        fprintf('Patient %d (HUP%s): No interictal trials defined or found. Skipping.\n', Patient, Patient_initials);
        continue;
    end

    for trial = trials_to_process
        fprintf('  Processing Patient %d (HUP%s), Trial %d...\n', Patient, Patient_initials, trial);
        
        % --- New Inner Loop: Iterate through each expected epoch (1 to 5) ---
        epochs_found_count = 0;
        for epoch_idx = 1:MAX_EPOCHS_TO_CHECK
            
            % Construct filename for the epoched .mat file
            filename = sprintf('HUP%s_SEEG_InterIctal_Processed_Trial%02d_Epoch%02d_Bipolar_GM_only.mat', Patient_initials, trial, epoch_idx);
            filepath = fullfile(PREPROCESSED_DATA_DIR, filename);

            % Check if the file exists (to stop when all epochs for a trial are processed)
            if ~exist(filepath, 'file')
                if epoch_idx == 1
                    % Only issue a warning if the first epoch file is missing
                    fprintf('    WARNING: No epoched files found for Trial %d (e.g. %s). Skipping trial.\n', trial, filename);
                else
                    fprintf('    Stop: Epoch %d not found. Finished epochs for Trial %d.\n', epoch_idx, trial);
                end
                break; % Exit epoch loop for current trial
            end
            
            % Load processed data including all metadata variables
            try
                % Load metadata from the single epoch file
                loaded_data = load(filepath, 'channels_id', 'channel_info_epoch', ...
                                   'electrode_loc_info_epoch', 'channel_status_description', ...
                                   'clinical_info_epoch');
            catch ME
                fprintf(2, '    ERROR loading metadata from %s: %s\n', filepath, ME.message);
                continue; % Skip to next epoch
            end

            % Concatenate loaded data to the respective master variables
            % CRITICAL: All concatenation must happen inside this epoch loop!
            if isfield(loaded_data, 'channel_info_epoch') && ~isempty(loaded_data.channel_info_epoch)
                all_channel_info_epoch = [all_channel_info_epoch, loaded_data.channel_info_epoch];
            end
            if isfield(loaded_data, 'electrode_loc_info_epoch') && ~isempty(loaded_data.electrode_loc_info_epoch)
                all_electrode_loc_info_epoch = [all_electrode_loc_info_epoch, loaded_data.electrode_loc_info_epoch];
            end
            if isfield(loaded_data, 'channel_status_description') && ~isempty(loaded_data.channel_status_description)
                % Ensure correct concatenation for cell array
                all_channel_status_description = [all_channel_status_description, loaded_data.channel_status_description];
            end
            if isfield(loaded_data, 'channels_id') && ~isempty(loaded_data.channels_id)
                % Ensure correct concatenation for cell array
                all_channels_id = [all_channels_id, loaded_data.channels_id];
            end
            
            % Process and concatenate clinical info
            if isfield(loaded_data, 'clinical_info_epoch') && ~isempty(loaded_data.clinical_info_epoch)
                current_clinical_info = loaded_data.clinical_info_epoch;
                % Ensure it's a row vector for horizontal concatenation (assuming it was saved as column)
                if iscolumn(current_clinical_info)
                    current_clinical_info = current_clinical_info'; 
                end
                all_clinical_info_epoch = [all_clinical_info_epoch, current_clinical_info];
            end
            
            epochs_found_count = epochs_found_count + 1;
            fprintf('    Concatenated Epoch %02d.\n', epoch_idx);

            % Clear temporary variables for the next epoch iteration
            clear loaded_data;
        end % End of epoch_idx loop

        if epochs_found_count > 0
             fprintf('  Finished processing Patient %d (HUP%s), Trial %d. Total epochs concatenated: %d.\n', Patient, Patient_initials, trial, epochs_found_count);
        end
    end % End of trial loop
end % End of patient loop

fprintf('\nFinished concatenating metadata.\n');

% Set the filename for the output .mat file containing concatenated metadata
concatenated_output_filename = fullfile(OUTPUT_DIR, 'ConcatenatedChannelMetadata_more_epochs.mat');

% Save the concatenated metadata variables
fprintf('Saving concatenated metadata to %s...\n', concatenated_output_filename);
save(concatenated_output_filename, ...
     'all_channel_info_epoch', ...
     'all_electrode_loc_info_epoch', ...
     'all_channel_status_description', ...
     'all_channels_id', ...
     'all_clinical_info_epoch', ... 
     '-v7.3'); % Use -v7.3 for large files

fprintf('Metadata concatenation and saving complete.\n');