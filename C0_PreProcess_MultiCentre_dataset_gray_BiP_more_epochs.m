clc;
clear all;
close all;

% Adding toolboxes (Ensure FieldTrip is correctly installed and accessible at this path)
addpath(genpath('/home/uqhkarim/Normative_HCTSA/fieldtrip-20250523'))
addpath(genpath('/home/uqhkarim/Normative_HCTSA/Coordi_to_atlases'))
addpath(genpath('/home/uqhkarim/Normative_HCTSA/NIfTI_20140122'))
output_folder='/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';
ft_defaults % Initialize FieldTrip defaults, clearing previous paths if needed

% --- Configuration Constants for Pre-processing ---
TARGET_FS = 200; % Target sampling rate in Hz for resampling
LP_FILTER_FREQ = 80; % Low-pass filter frequency for anti-aliasing in Hz
HP_FILTER_FREQ = 0.5; % High-pass filter frequency for line noise pre-processing in Hz
HP_FILTER_ORDER = 1000; % FIR filter order for high-pass
LINE_FREQ = 60; % Fundamental power line frequency in Hz
EPOCH_DURATION_SEC = 60; % Desired duration of each extracted epoch in seconds
MAX_BIPOLAR_DISTANCE_MM = 10; % NEW: Maximum allowed distance for bipolar pairing in mm

% --- Analysis Type (fixed to interictal as per request) ---
ictal_or_inter = 'interictal';    

%% Data Processing Loop
for Patient = [1]
    % Patient-specific data definitions
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
        Patient_initials='180'; modality='seeg'; inter_ictal_trials=[1:2];        
    elseif Patient==51
        Patient_initials='181'; modality='seeg'; inter_ictal_trials=[1:2];        
    elseif Patient==52
        Patient_initials='185'; modality='seeg'; inter_ictal_trials=[1:2];        
    elseif Patient==53
        Patient_initials='187'; modality='seeg'; inter_ictal_trials=[1:2];        
    elseif Patient==54
        Patient_initials='188'; modality='seeg'; inter_ictal_trials=[1:2];        
    elseif Patient==55
        Patient_initials='190'; modality='seeg'; inter_ictal_trials=[1:2];        
    end

    patient_data_address = ['/QRISdata/Q8553/Normative_project/sub-HUP',Patient_initials,'/ses-presurgery/ieeg/'];

    trials = inter_ictal_trials; % Always use the interictal trials list

    % Skip the current patient if no interictal trials are defined
    if isempty(trials)
        fprintf('Patient %d (HUP%s): No interictal trials defined or found. Skipping.\n', Patient, Patient_initials);
        continue;    
    end
    
    % Loop through each specified interictal trial for the current patient
    for trial = trials
        fprintf('Processing Patient %d (HUP%s), Trial %d...\n', Patient, Patient_initials, trial);

        % Construct file paths
        folder_address = [patient_data_address, 'sub-HUP', Patient_initials, '_ses-presurgery_task-', ictal_or_inter, '_acq-', modality, '_run-', sprintf('%02d',trial)];
        signal_address = [folder_address, '_ieeg.edf'];
        channel_tsv_path = [folder_address, '_channels.tsv'];
        electrode_tsv_path = [patient_data_address, 'sub-HUP', Patient_initials, '_ses-presurgery_acq-', modality, '_space-fsaverage_electrodes.tsv'];
        
        % --- File Existence Checks ---
        if ~exist(signal_address, 'file')
            fprintf('    WARNING: Signal file not found: %s. Skipping this trial.\n', signal_address);
            continue;
        end
        if ~exist(channel_tsv_path, 'file')
            fprintf('    WARNING: Channel info file not found: %s. Skipping this trial.\n', channel_tsv_path);
            continue;
        end
        if ~exist(electrode_tsv_path, 'file')
            fprintf('    WARNING: Electrode location info file not found: %s. Skipping this trial.\n', electrode_tsv_path);
            continue;
        end

        % Load header and TSV files
        header = ft_read_header(signal_address);
        channel_info = tdfread(channel_tsv_path);
        electrode_loc_info = tdfread(electrode_tsv_path);
        
        % --- Initial Channel Selection (excluding bad, ECG/EKG, and n/a channels) ---
        % Store info for valid individual contacts before bipolarization
        individual_channels_info = struct('label', {}, 'original_idx', {}, 'status', {}, 'group', {}, ...
                                          'contact_num', {}, ... 
                                          'x', {}, 'y', {}, 'z', {}, 'is_gray_matter', {});
        
        % Iterate through all channels in the header
        for channel_num = 1:size(header.label, 1)
            current_label = header.label{channel_num};
            
            % Ensure the channel_info and electrode_loc_info have data for this channel_num
            if channel_num > size(channel_info.status_description, 1) || channel_num > size(electrode_loc_info.x, 1) % Using status_description here
                fprintf('    SKIPPING Channel %s: Channel info or electrode location info incomplete for this channel number.\n', current_label);
                continue; % Skip if info is missing for this channel
            end

            % Convert to char arrays for string comparison (tdfread loads strings as char arrays padded with spaces)
            % IMPORTANT: Using channel_info.status_description as per user's instruction
            current_status_good_bad = strtrim(channel_info.status(channel_num, :)); 
            current_status = strtrim(channel_info.status_description(channel_num, :)); 
            current_electrode_x_str = strtrim(electrode_loc_info.x(channel_num, :)); 
            
            % Apply exclusion criteria for individual contacts
            is_bad_channel = strcmp(current_status_good_bad, 'bad');
            is_ecg_channel = contains(current_label, 'EKG', 'IgnoreCase', true) || contains(current_label, 'ECG', 'IgnoreCase', true);
            is_na_electrode_loc = strcmp(current_electrode_x_str, 'n/a'); % Assuming 'n/a' indicates bad location

            if ~is_bad_channel && ~is_ecg_channel && ~is_na_electrode_loc
                % Extract MNI coordinates
                x_coord = str2double(current_electrode_x_str);
                y_coord = str2double(strtrim(electrode_loc_info.y(channel_num, :)));
                z_coord = str2double(strtrim(electrode_loc_info.z(channel_num, :)));

                % Remove unwanted strings from label for grouping
                cleaned_label = strrep(strrep(current_label, 'EEG ', ''), '-Ref', '');
                
                % Determine electrode group (e.g., 'LEF', 'RIN')
                % IMPROVED REGEX: Captures anything before the last digits as group, last digits as contact_num
                token = regexp(cleaned_label, '^(.*?)(\d+)$', 'tokens', 'once'); 
                electrode_group = '';
                contact_number = NaN;
                
                % Robust checks for the output of regexp 
                if iscell(token) && ~isempty(token) && iscell(token{1}) && numel(token{1}) == 2
                    electrode_group = token{1}{1};
                    contact_number = str2double(token{1}{2});
                else
                    % Fallback: try to find just a letter-based group if the primary regex fails
                    token_letters_only = regexp(cleaned_label, '([A-Za-z]+)', 'tokens', 'once');
                    
                    if iscell(token_letters_only) && ~isempty(token_letters_only) && iscell(token_letters_only{1}) && numel(token_letters_only{1}) >= 1
                        electrode_group = token_letters_only{1}{1};
                        contact_number = 0; % Indicate no specific contact number found for this type
                    else
                        % If no identifiable letter group or number, use the full cleaned label as the group
                        electrode_group = cleaned_label; 
                        contact_number = 0; % Assign 0 or handle as a single contact electrode
                        fprintf('            WARNING: Could not parse electrode group from label: %s. Using full label as group and contact_num as 0.\n', cleaned_label);
                    end
                end

                % Add to a temporary structure for further processing
                individual_channels_info(end+1) = struct('label', current_label, ...
                                                         'original_idx', channel_num, ...
                                                         'status', current_status, ...
                                                         'group', electrode_group, ...
                                                         'contact_num', contact_number, ... 
                                                         'x', x_coord, 'y', y_coord, 'z', z_coord, ...
                                                         'is_gray_matter', false); 
            else
                fprintf('    EXCLUDING Initial Channel: %s (Status: %s, EKG/ECG: %d, Loc_NA: %d)\n', ...
                    current_label, current_status_good_bad, is_ecg_channel, is_na_electrode_loc);
            end
        end

        if isempty(individual_channels_info)
            fprintf('    WARNING: No valid individual channels found after initial exclusion for Patient %d, Trial %d. Skipping this trial.\n', Patient, trial);
            clearvars -except Patient Patient_initials patient_data_address inter_ictal_trials ictal_or_inter modality ...
                               TARGET_FS LP_FILTER_FREQ HP_FILTER_FREQ HP_FILTER_ORDER LINE_FREQ EPOCH_DURATION_SEC ...
                               MAX_BIPOLAR_DISTANCE_MM output_folder;
            continue;
        end

        % --- Perform Gray Matter Exclusion for individual channels ---
        fprintf('    Step 1/5: Determining gray matter contacts for individual channels...\n');
        
        num_initial_valid_channels = length(individual_channels_info);
        for k = 1:num_initial_valid_channels
            x_coord = individual_channels_info(k).x;
            y_coord = individual_channels_info(k).y;
            z_coord = individual_channels_info(k).z;
            currnt_elec_coord = [x_coord, y_coord, z_coord];
            
            % Check for NaN coordinates (from 'n/a' or other issues during str2double)
            if any(isnan(currnt_elec_coord))
                fprintf('        WARNING: MNI coordinates for channel %s are NaN. Skipping for atlas check and treating as non-gray matter.\n', ...
                    individual_channels_info(k).label);
                individual_channels_info(k).is_gray_matter = false; % Explicitly mark as false
                continue; 
            end
            try
                [atlas_info] = mni2atlas(currnt_elec_coord, [2 3], '1mm'); % Using [2 3] for 'tissue' and 'lobe'
                
                if length(atlas_info) >= 2 && ~isempty(atlas_info(2).label)
                    tissue_label = atlas_info(2).label{1,1};
                    if ~contains(tissue_label, 'White', 'IgnoreCase', true)
                        individual_channels_info(k).is_gray_matter = true; % Mark as gray matter
                    else
                        fprintf('        INDIVIDUAL EXCLUDING Channel %s (MNI: [%.1f %.1f %.1f]): Classified as White Matter.\n', ...
                            individual_channels_info(k).label, x_coord, y_coord, z_coord);
                    end
                else
                    fprintf('        WARNING: Atlas tissue classification not available for channel %s (MNI: [%.1f %.1f %.1f]). Treating as non-gray matter.\n', ...
                        individual_channels_info(k).label, x_coord, y_coord, z_coord);
                end
            catch ME
                fprintf('        ERROR in mni2atlas for channel %s (MNI: [%.1f %.1f %.1f]): %s. Treating as non-gray matter.\n', ...
                    individual_channels_info(k).label, x_coord, y_coord, z_coord, ME.message);
                individual_channels_info(k).is_gray_matter = false;
            end
        end
        % Filter individual_channels_info to keep only gray matter contacts
        individual_channels_info = individual_channels_info([individual_channels_info.is_gray_matter]);

        if isempty(individual_channels_info)
            fprintf('    WARNING: No gray matter contacts found for Patient %d, Trial %d after MNI to atlas conversion. Skipping this trial.\n', Patient, trial);
            clearvars -except Patient Patient_initials patient_data_address inter_ictal_trials ictal_or_inter modality ...
                               TARGET_FS LP_FILTER_FREQ HP_FILTER_FREQ HP_FILTER_ORDER LINE_FREQ EPOCH_DURATION_SEC ...
                               MAX_BIPOLAR_DISTANCE_MM output_folder;
            continue; 
        end
        fprintf('    Kept %d gray matter contacts out of %d initially valid channels.\n', length(individual_channels_info), num_initial_valid_channels);

        % --- Load raw data for all initially selected good channels (which are now GM contacts) ---
        cfg_load = [];
        cfg_load.dataset = signal_address;
        cfg_load.continuous = 'yes';
        
        channels_to_load_labels = {individual_channels_info.label}; 
        
        if isempty(channels_to_load_labels)
             fprintf('    WARNING: No valid individual channels with good location info to load data for Patient %d, Trial %d. Skipping this trial.\n', Patient, trial);
             clearvars -except Patient Patient_initials patient_data_address inter_ictal_trials ictal_or_inter modality ...
                               TARGET_FS LP_FILTER_FREQ HP_FILTER_FREQ HP_FILTER_ORDER LINE_FREQ EPOCH_DURATION_SEC ...
                               MAX_BIPOLAR_DISTANCE_MM output_folder;
            continue;
        end

        cfg_load.channel = unique(channels_to_load_labels); 
        fprintf('    Loading raw data for %d gray matter contacts for bipolarization...\n', length(cfg_load.channel));
        data_raw = ft_preprocessing(cfg_load);

        % Create a mapping from individual_channels_info labels to data_raw.label indices
        [~, raw_to_loaded_idx] = ismember({individual_channels_info.label}, data_raw.label);
        
        if any(raw_to_loaded_idx == 0)
            fprintf('    CRITICAL WARNING: Mismatch between selected GM channels and loaded data labels. Some channels might be missing from ft_preprocessing output. Skipping this trial.\n');
             clearvars -except Patient Patient_initials patient_data_address inter_ictal_trials ictal_or_inter modality ...
                               TARGET_FS LP_FILTER_FREQ HP_FILTER_FREQ HP_FILTER_ORDER LINE_FREQ EPOCH_DURATION_SEC ...
                               MAX_BIPOLAR_DISTANCE_MM output_folder;
            continue;
        end
        % --- Create Bipolar Channels based on Euclidean Distance and 10mm constraint ---
        fprintf('    Step 2/5: Creating bipolar montage based on closest MNI coordinates (max %gmm distance)...\n', MAX_BIPOLAR_DISTANCE_MM);

        num_gm_contacts = length(individual_channels_info);
        
        bipolar_signal_data = [];
        channels_id = {};
        channel_info_epoch = struct('label', {}, 'status', {}, 'group', {});
        electrode_loc_info_epoch = struct('name', {}, 'x', {}, 'y', {}, 'z', {});

        if num_gm_contacts < 2
            fprintf('        WARNING: Less than 2 gray matter contacts available for bipolarization. Skipping bipolarization.\n');
        else
            % Generate all unique pairs and their distances, filtered by MAX_BIPOLAR_DISTANCE_MM
            eligible_pairs = struct('idx1', {}, 'idx2', {}, 'distance', {});
            for i = 1:num_gm_contacts-1
                for j = i+1:num_gm_contacts
                    contact1_info = individual_channels_info(i);
                    contact2_info = individual_channels_info(j);
                    
                    dist = sqrt((contact2_info.x - contact1_info.x)^2 + ...
                                (contact2_info.y - contact1_info.y)^2 + ...
                                (contact2_info.z - contact1_info.z)^2);
                    
                    if dist <= MAX_BIPOLAR_DISTANCE_MM % Apply the 10mm distance constraint
                        eligible_pairs(end+1) = struct('idx1', i, 'idx2', j, 'distance', dist);
                    end
                end
            end

            if isempty(eligible_pairs)
                fprintf('        WARNING: No bipolar pairs found within %gmm distance for Patient %d, Trial %d. Skipping bipolarization.\n', MAX_BIPOLAR_DISTANCE_MM, Patient, trial);
                clearvars -except Patient Patient_initials patient_data_address inter_ictal_trials ictal_or_inter modality ...
                               TARGET_FS LP_FILTER_FREQ HP_FILTER_FREQ HP_FILTER_ORDER LINE_FREQ EPOCH_DURATION_SEC ...
                               MAX_BIPOLAR_DISTANCE_MM output_folder;
                continue; % Skip the rest of the bipolar processing for this trial
            end

            % Sort eligible pairs by distance
            [~, sorted_indices] = sort([eligible_pairs.distance]);
            sorted_pairs = eligible_pairs(sorted_indices);

            % Select non-overlapping bipolar channels (greedy approach)
            used_contact_indices = false(1, num_gm_contacts); % Track which contacts have been used

            for p = 1:length(sorted_pairs)
                current_pair = sorted_pairs(p);
                contact_idx1 = current_pair.idx1;
                contact_idx2 = current_pair.idx2;

                if ~used_contact_indices(contact_idx1) && ~used_contact_indices(contact_idx2)
                    contact1 = individual_channels_info(contact_idx1);
                    contact2 = individual_channels_info(contact_idx2);
                    
                    % Get data for the two contacts from the loaded data_raw
                    idx1_in_loaded_data = raw_to_loaded_idx(contact_idx1); % Use the pre-computed map
                    idx2_in_loaded_data = raw_to_loaded_idx(contact_idx2);

                    if isempty(idx1_in_loaded_data) || isempty(idx2_in_loaded_data) || idx1_in_loaded_data == 0 || idx2_in_loaded_data == 0
                        fprintf('        WARNING: Could not find loaded data for contact pair %s-%s. Skipping bipolar pair.\n', contact1.label, contact2.label);
                        continue;
                    end

                    signal1 = data_raw.trial{1}(idx1_in_loaded_data, :);
                    signal2 = data_raw.trial{1}(idx2_in_loaded_data, :);
                    
                    bipolar_signal = signal2 - signal1; 
                    bipolar_signal_data = [bipolar_signal_data; bipolar_signal];

                    % Bipolar Channel Naming
                    name1 = strrep(strrep(contact1.label, 'EEG ', ''), '-Ref', '');
                    name2 = strrep(strrep(contact2.label, 'EEG ', ''), '-Ref', '');
                    bipolar_name = [name1, '-', name2];
                    channels_id{end+1} = bipolar_name;

                    % Bipolar Coordinate Averaging
                    bipolar_x = (contact1.x + contact2.x) / 2;
                    bipolar_y = (contact1.y + contact2.y) / 2;
                    bipolar_z = (contact1.z + contact2.z) / 2;
                    
                    electrode_loc_info_epoch(end+1) = struct('name', bipolar_name, ...
                                                                      'x', bipolar_x, ...
                                                                      'y', bipolar_y, ...
                                                                      'z', bipolar_z);
                    
                    % Bipolar Status Logic 
                    status1 = contact1.status; % This now comes from channel_info.status_description
                    status2 = contact2.status; % This now comes from channel_info.status_description
                    bipolar_status = 'n/a'; 
                    if (strcmp(status1, 'resect,soz') || strcmp(status2, 'resect,soz')) || ...
                       (strcmp(status1, 'soz') && strcmp(status2, 'resect')) || ...
                       (strcmp(status1, 'resect') && strcmp(status2, 'soz'))
                        bipolar_status = 'resect,soz';
                    elseif (strcmp(status1, 'resect,soz') && strcmp(status2, 'n/a')) || ...
                           (strcmp(status2, 'n/a') && strcmp(status1, 'resect,soz'))
                        bipolar_status = 'resect,soz';
                    elseif (strcmp(status1, 'soz') && strcmp(status2, 'n/a')) || ...
                           (strcmp(status2, 'soz') && strcmp(status1, 'n/a'))
                        bipolar_status = 'soz';
                    elseif (strcmp(status1, 'resect') && strcmp(status2, 'n/a')) || ...
                           (strcmp(status2, 'resect') && strcmp(status1, 'n/a'))
                        bipolar_status = 'resect';                        
                    elseif strcmp(status1, status2) 
                        bipolar_status = status1;
                    else
                        fprintf('        WARNING: Unhandled status combination for %s-%s: %s vs %s. Setting status to n/a.\n', ...
                            contact1.label, contact2.label, status1, status2);
                        bipolar_status = 'n/a';
                    end

                    channel_info_epoch(end+1) = struct('label', bipolar_name, ...
                                                               'status', bipolar_status, ...
                                                               'group', contact1.group); % Assign group from first contact
                    % Mark contacts as used
                    used_contact_indices(contact_idx1) = true;
                    used_contact_indices(contact_idx2) = true;
                end
            end % End of sorted_pairs loop
        end % End of num_gm_contacts check

        if isempty(bipolar_signal_data)
            fprintf('    WARNING: No bipolar channels formed for Patient %d, Trial %d. Skipping save for this trial.\n', Patient, trial);
            clearvars -except Patient Patient_initials patient_data_address inter_ictal_trials ictal_or_inter modality ...
                               TARGET_FS LP_FILTER_FREQ HP_FILTER_FREQ HP_FILTER_ORDER LINE_FREQ EPOCH_DURATION_SEC ...
                               MAX_BIPOLAR_DISTANCE_MM output_folder;
            continue;
        end

        % --- Update data_raw for pre-processing with bipolar data ---
        % Replace data_raw with the new bipolar data and labels
        data_raw.trial{1} = bipolar_signal_data;
        data_raw.label = channels_id;
        data_raw.fsample = header.Fs; % Original sampling rate before resampling
        data_raw.time{1} = (0:size(bipolar_signal_data,2)-1) / data_raw.fsample; % Update time axis

        % --- Pre-processing steps for bipolar data ---
        % 1. Resample signals to TARGET_FS (200 Hz) with anti-aliasing low-pass filter at LP_FILTER_FREQ (80 Hz)
        fprintf('    Step 3/5: Resampling bipolar data from %.2f Hz to %d Hz with anti-aliasing LPF at %d Hz...\n', ...
            data_raw.fsample, TARGET_FS, LP_FILTER_FREQ);
        cfg_resample = [];
        cfg_resample.resamplefs = TARGET_FS;
        cfg_resample.detrend = 'no'; % Demeaning will be done later as a final step
        cfg_resample.lpfilter = 'yes';
        cfg_resample.lpfreq = LP_FILTER_FREQ;
        data_resampled = ft_resampledata(cfg_resample, data_raw);

        % 2. Power-line interference reduction
        fprintf('    Step 4/5: Reducing power-line interference (60Hz harmonics via DFT filter) on bipolar data...\n');
        
        % Remove 3rd harmonic (180 Hz)
        data_resampled.trial{1} = ft_preproc_dftfilter(data_resampled.trial{1}, data_resampled.fsample, [LINE_FREQ * 3]);
        % Remove 2nd harmonic (120 Hz)
        data_resampled.trial{1} = ft_preproc_dftfilter(data_resampled.trial{1}, data_resampled.fsample, [LINE_FREQ * 2]);
        % Remove fundamental (60 Hz)
        data_resampled.trial{1} = ft_preproc_dftfilter(data_resampled.trial{1}, data_resampled.fsample, [LINE_FREQ]);

        % 3. High-pass filter at HP_FILTER_FREQ (0.5 Hz) (FIR filter, order 1000)
        fprintf('    Applying high-pass filter at %g Hz (FIR order %d) on bipolar data...\n', HP_FILTER_FREQ, HP_FILTER_ORDER);
        data_resampled.trial{1} = ft_preproc_highpassfilter(data_resampled.trial{1}, data_resampled.fsample, ...
                                                            HP_FILTER_FREQ, HP_FILTER_ORDER, 'fir');

        % 4. Subtract the mean value of each signal channel (already bipolar)
        fprintf('    Demeaning each bipolar signal channel...\n');
        data_resampled.trial{1} = data_resampled.trial{1} - mean(data_resampled.trial{1}, 2);

        % Now, the processed data is in data_resampled.trial{1}.
        data_processed = data_resampled;

         % --- MODIFIED SECTION: Segmenting the processed data into multiple 1-minute epochs ---
        fprintf('    Step 5/5: Segmenting data into %d-second epochs and saving...\n', EPOCH_DURATION_SEC);

        samples_per_epoch = EPOCH_DURATION_SEC * data_processed.fsample;
        current_full_signal = data_processed.trial{1};
        signal_length_samples = size(current_full_signal, 2);
        
        % Calculate the total number of full epochs we can extract
        num_epochs = floor(signal_length_samples / samples_per_epoch);
        
        if num_epochs == 0
            fprintf('        WARNING: Processed bipolar recording length (%.2f sec) is less than the required %d seconds for even one epoch. Skipping save for this trial.\n', ...
                signal_length_samples / data_processed.fsample, EPOCH_DURATION_SEC);
        else
            fprintf('        Found %.2f minutes of data, extracting %d epoch(s).\n', ...
                signal_length_samples / data_processed.fsample / 60, num_epochs);
                
            % Prepare channel_status_description once, as it's the same for all epochs in this trial
            channel_status_description = {channel_info_epoch.status};
    
            % Loop through and save each epoch
            for epoch_idx = 1:num_epochs
                start_sample = (epoch_idx - 1) * samples_per_epoch + 1;
                end_sample = epoch_idx * samples_per_epoch;
                
                % Extract the data for the current epoch
                signal_epoch = current_full_signal(:, start_sample:end_sample);
                
                % --- Saving the processed signals and channel information for this epoch ---
                % Create a unique output filename for each epoch
                output_filename = sprintf('HUP%s_SEEG_InterIctal_Processed_Trial%02d_Epoch%02d_Bipolar_GM_only.mat', Patient_initials, trial, epoch_idx);
                fprintf('        Saving Epoch %d/%d to: %s\n', epoch_idx, num_epochs, output_filename);
                save([output_folder, output_filename], 'signal_epoch', 'channels_id', ...
                     'channel_info_epoch', 'electrode_loc_info_epoch', ...
                     'channel_status_description', 'TARGET_FS', '-v7.3');
            end
        end
        % --- END OF MODIFIED SECTION ---
        fprintf('    Finished processing Patient %d (HUP%s), Trial %d.\n\n', Patient, Patient_initials, trial);
        
        % Clear variables specific to this trial to free up memory before the next iteration
        clearvars -except output_folder Patient Patient_initials patient_data_address inter_ictal_trials ictal_or_inter modality ...
                               TARGET_FS LP_FILTER_FREQ HP_FILTER_FREQ HP_FILTER_ORDER LINE_FREQ EPOCH_DURATION_SEC ...
                               MAX_BIPOLAR_DISTANCE_MM; 
    end % End of trial loop
end % End of patient loop

fprintf('Script completed.\n');