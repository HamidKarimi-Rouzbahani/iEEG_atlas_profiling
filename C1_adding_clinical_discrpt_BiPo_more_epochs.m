clc;
clear all;
close all;

% Adding toolboxes (Ensure FieldTrip is correctly installed and accessible at this path)
addpath(genpath('/home/uqhkarim/Normative_HCTSA/fieldtrip-20250523'))
addpath(genpath('/home/uqhkarim/Normative_HCTSA/Coordi_to_atlases'))
addpath(genpath('/home/uqhkarim/Normative_HCTSA/NIfTI_20140122'))
output_folder='/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';
ft_defaults % Initialize FieldField defaults, clearing previous paths if needed

% --- Configuration Constants for Pre-processing ---
TARGET_FS = 200; % Target sampling rate in Hz for resampling
LP_FILTER_FREQ = 80; % Low-pass filter frequency for anti-aliasing in Hz
HP_FILTER_FREQ = 0.5; % High-pass filter frequency for line noise pre-processing in Hz
HP_FILTER_ORDER = 1000; % FIR filter order for high-pass
LINE_FREQ = 60; % Fundamental power line frequency in Hz
EPOCH_DURATION_SEC = 60; % Desired duration of each extracted epoch in seconds
MAX_BIPOLAR_DISTANCE_MM = 10; % Maximum allowed distance for bipolar pairing in mm
MAX_EPOCHS_TO_CHECK = 5; % The maximum number of epochs expected to be processed per trial

% --- Analysis Type (fixed to interictal as per request) ---
ictal_or_inter = 'interictal';    

%% --- Load Participant Clinical Data (No change in this section) ---
TSV_FILE_PATH = '/QRISdata/Q8553/Normative_project/participants.tsv';
fprintf('Loading participant clinical data from: %s\n', TSV_FILE_PATH);

% Initialize clinical_table as empty in case file doesn't exist or is empty
clinical_table = table(); 
empty_clinical_struct = struct('participant_id', '', 'age', NaN, 'sex', '', 'hand', '', ...
                               'outcome', '', 'engel', '', 'therapy', '', 'implant', '', ...
                               'target', '', 'lesion_status', '', 'age_onset', NaN, 'PatientID_Num', '');

if exist(TSV_FILE_PATH, 'file')
    try
        clinical_table = readtable(TSV_FILE_PATH, 'FileType','text','Delimiter', '\t', 'TreatAsEmpty', {'n/a', 'N/A'});
        clinical_table.PatientID_Num = cellfun(@(x) x(8:end), clinical_table.participant_id, 'UniformOutput', false);
        
        if ~isempty(clinical_table)
            temp_struct_for_fields = table2struct(clinical_table(1,:));
            fields = fieldnames(temp_struct_for_fields);
            for k = 1:length(fields)
                if ischar(temp_struct_for_fields.(fields{k}))
                    empty_clinical_struct.(fields{k}) = '';
                elseif isnumeric(temp_struct_for_fields.(fields{k}))
                    empty_clinical_struct.(fields{k}) = NaN;
                else
                    empty_clinical_struct.(fields{k}) = [];
                end
            end
        end
        fprintf('Participant clinical data loaded successfully.\n');
    catch ME
        warning('Error loading participant clinical data from %s: %s. Clinical information will not be added.', TSV_FILE_PATH, ME.message);
        clinical_table = table();
    end
else
    warning('Participant clinical data file not found at: %s. Clinical information will not be added.', TSV_FILE_PATH);
    clinical_table = table();
end

%% Data Processing Loop (Modified to include Epoch loop)
for Patient = [1:56]
    % Patient-specific data definitions (Your definitions remain here)
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
    
    % --- Get clinical data for current patient ---
    patient_clinical_struct = empty_clinical_struct;
    if ~isempty(clinical_table)
        current_patient_clinical_idx = find(strcmp(clinical_table.PatientID_Num, Patient_initials));
        if ~isempty(current_patient_clinical_idx)
            patient_clinical_struct = table2struct(clinical_table(current_patient_clinical_idx, :));
        else
            warning('Clinical data not found in TSV for Patient HUP%s.', Patient_initials);
        end
    end

    % --- Skip patient if no interictal trials are defined ---
    if isempty(trials_to_process)
        fprintf('Patient %d (HUP%s): No interictal trials defined. Skipping.\n', Patient, Patient_initials);
        continue;
    end

    for trial = trials_to_process
        fprintf('  Processing Patient %d (HUP%s), Trial %d...\n', Patient, Patient_initials, trial);
        
        % --- New Inner Loop: Iterate through each expected epoch ---
        epochs_processed_count = 0;
        
        for epoch_idx = 1:MAX_EPOCHS_TO_CHECK % Loop 1 to 5 times
            
            % Construct the full filename (for loading and saving)
            filename = sprintf('HUP%s_SEEG_InterIctal_Processed_Trial%02d_Epoch%02d_Bipolar_GM_only.mat', Patient_initials, trial, epoch_idx);
            filepath = fullfile(output_folder, filename);

            % Check if the epoched file exists. If not, break the epoch loop for this trial.
            if ~exist(filepath, 'file')
                if epoch_idx == 1
                     fprintf('    WARNING: No epoched files found for Trial %d (e.g. %s). Skipping trial.\n', trial, filename);
                else
                     fprintf('    Stop: Epoch %d not found for Trial %d. Finished epochs for this trial.\n', epoch_idx, trial);
                end
                break; % Exit epoch loop for current trial
            end
            
            % Load processed data
            try
                % Load signal_epoch (now expected to be 2D: Channels x Time)
                loaded_data = load(filepath, 'signal_epoch', 'channels_id', 'channel_info_epoch', 'electrode_loc_info_epoch');
            catch ME
                fprintf(2, '    ERROR loading data from %s: %s\n', filepath, ME.message);
                continue; % Skip to next epoch
            end

            % Assign loaded variables
            signal_epoch = loaded_data.signal_epoch; % Now 2D
            channels_id = loaded_data.channels_id;
            channel_info_epoch = loaded_data.channel_info_epoch;
            electrode_loc_info_epoch = loaded_data.electrode_loc_info_epoch;

            % --- Add Clinical and Status Info (identical for all epochs of the same trial) ---
            num_channels = size(signal_epoch, 1);
            
            % Prepare channel_status_description
            if ~isempty(channel_info_epoch)
                channel_status_description = {channel_info_epoch.status};
            else
                channel_status_description = {};
            end

            % Create clinical_info_epoch (duplicate patient info for each electrode)
            clinical_info_epoch = repmat(patient_clinical_struct, num_channels, 1);

            % --- Saving the updated file (using the same filename) ---
            output_filename = filename; % The filename is already correctly constructed

            fprintf('    Saving updated Epoch %02d to: %s\n', epoch_idx, output_filename);
            
            % Save all original variables plus the new 'clinical_info_epoch'
            save([output_folder, output_filename], 'signal_epoch', 'channels_id', ...
                 'channel_info_epoch', 'electrode_loc_info_epoch', ...
                 'channel_status_description', 'clinical_info_epoch', ... 
                 'TARGET_FS', '-v7.3', '-append'); % Use -append to only add/overwrite new variables
            
            epochs_processed_count = epochs_processed_count + 1;
            
            % Clear temporary variables for the next epoch iteration
            clear signal_epoch channels_id channel_info_epoch electrode_loc_info_epoch loaded_data;
        end % End of epoch_idx loop

        fprintf('  Finished Trial %d (HUP%s). Total Epochs Updated: %d.\n\n', trial, Patient_initials, epochs_processed_count);
        
        % Clear variables specific to this trial/patient before the next trial/patient
        clearvars -except output_folder Patient Patient_initials inter_ictal_trials ictal_or_inter modality ...
                               TARGET_FS LP_FILTER_FREQ HP_FILTER_FREQ HP_FILTER_ORDER LINE_FREQ EPOCH_DURATION_SEC ...
                               MAX_BIPOLAR_DISTANCE_MM TSV_FILE_PATH clinical_table empty_clinical_struct patient_clinical_struct MAX_EPOCHS_TO_CHECK; 
    end % End of trial loop
end % End of patient loop

fprintf('Script completed.\n');