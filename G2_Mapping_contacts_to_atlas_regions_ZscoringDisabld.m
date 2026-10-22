clc;
close all;
clear all;

addpath(genpath('/scratch/user/uqhkarim/Normative_HCTSA/Coordi_to_atlases'))
addpath(genpath('/scratch/user/uqhkarim/Normative_HCTSA/NIfTI_20140122'))

% Define output folder for saving results
output_folder = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';
% Ensure the output folder exists
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Loading Frauscher data
Frauscher_data = load('/scratch/user/uqhkarim/Normative_HCTSA/Scripts/Frauscher backup/HCTSA_W.mat');
Frauscher_annotations = load('/scratch/user/uqhkarim/Normative_HCTSA/Dataset/MNI_Frauscher/MatlabFile/MatlabFile.mat','ChannelPosition','ChannelRegion','ChannelName','RegionName');


% Loading DS4100 data
% Load the DS4100 annotations, including clinical info for all contacts
% DS4100_data = load('/home/uqhkarim/Normative_HCTSA/Scripts/LittsScripts/old/HCTSA.mat');
% DS4100_annotations_path = '/home/uqhkarim/Normative_HCTSA/Scripts/LittsScripts/old/ConcatenatedChannelMetadata_pruned.mat';
DS4100_data = load('/scratch/user/uqhkarim/Normative_HCTSA/Scripts/LittsScripts/HCTSA.mat');
DS4100_annotations_path = '/scratch/user/uqhkarim/Normative_HCTSA/Scripts/LittsScripts/ConcatenatedChannelMetadata_more_epochs.mat';

if exist(DS4100_annotations_path, 'file')
    DS4100_annotations = load(DS4100_annotations_path, ...
                                  'all_electrode_loc_info_epoch', ...
                                  'all_channel_status_description', ...
                                  'all_clinical_info_epoch', ...
                                  'all_channels_id', ...
                                  'all_channel_info_epoch');
else
    error('DS4100 ConcatenatedChannelMetadata_pruned.mat not found at %s. Please ensure it exists.', DS4100_annotations_path);
end

pattern = '\d+\.?\d*|%'; % Pattern to remove numbers and percentage signs

% --- Process Frauscher_atlases ---
fprintf('Processing Frauscher atlases...\n');
tmp_locations_frauscher = Frauscher_annotations.ChannelPosition;
Frauscher_atlases_processed = cell(1, size(tmp_locations_frauscher, 1)); % Pre-allocate
for i = 1:size(tmp_locations_frauscher, 1)
    % Using [2] for 'tissue' - based on prior context/discussion on mni2atlas usage
    [atlas_info] = mni2atlas([tmp_locations_frauscher(i,1) tmp_locations_frauscher(i,2) tmp_locations_frauscher(i,3)], [2], '2mm');
    [a] = mni2atlas([tmp_locations_frauscher(i,1) tmp_locations_frauscher(i,2) tmp_locations_frauscher(i,3)], [1:6], '2mm');
    ccc
    if ~isempty(atlas_info.label)
        cleaned_text = regexprep(atlas_info.label{1,1}, pattern, '');
        Frauscher_atlases_processed{i} = strtrim(cleaned_text);
    else
        Frauscher_atlases_processed{i} = 'No ROI'; % Assign a placeholder if no ROI found
    end
end
fprintf('Frauscher atlas processing complete.\n');


% --- Process DS4100_atlases for ALL contacts ---
fprintf('Processing DS4100 atlases for all contacts...\n');
tmp_locations_ds4100 = DS4100_annotations.all_electrode_loc_info_epoch;
DS4100_atlases_processed = cell(1, length(tmp_locations_ds4100)); % Pre-allocate
for i = 1:length(tmp_locations_ds4100)
    [atlas_info] = mni2atlas([tmp_locations_ds4100(i).x tmp_locations_ds4100(i).y tmp_locations_ds4100(i).z], [2], '2mm');
    if ~isempty(atlas_info.label)
        cleaned_text = regexprep(atlas_info.label{1,1}, pattern, '');
        DS4100_atlases_processed{i} = strtrim(cleaned_text);
    else
        DS4100_atlases_processed{i} = 'No ROI'; % Assign a placeholder if no ROI found
    end
end
fprintf('DS4100 atlas processing complete for all contacts.\n');

% --- NO FILTERING OF DS4100 CONTACTS BASED ON OUTCOME ---
% As per your request, we will now process all DS4100 contacts
% and save them to allow for more flexible future analyses (e.g., multivariate regression).
DS4100_TS_DataMat_all_contacts = DS4100_data.TS_DataMat;

% Initialize output storage for normalized data for both datasets
% Renamed variables to reflect 0-1 range normalization
DS4100 = NaN(size(DS4100_TS_DataMat_all_contacts));
% DS4100_abs_normalized_per_roi = NaN(size(DS4100_TS_DataMat_all_contacts));

Frauscher = NaN(size(Frauscher_data.TS_DataMat));
% Frauscher_abs_normalized_per_roi = NaN(size(Frauscher_data.TS_DataMat));

fprintf('Proceeding with all %d DS4100 contacts and %d Frauscher contacts for analysis.\n', ...
    size(DS4100_TS_DataMat_all_contacts, 1), size(Frauscher_data.TS_DataMat, 1));


% --- Identify unique ROIs across both datasets ---
% Use DS4100_atlases_processed (containing all contacts) for ROI identification
all_rois = unique([Frauscher_atlases_processed, DS4100_atlases_processed]);
% Remove 'No ROI' or empty strings if they exist and are not desired for comparison
all_rois = all_rois(~ismember(all_rois, {'No ROI', ''}));

fprintf('Identified %d unique ROIs for analysis.\n', length(all_rois));

% --- Loop through each unique ROI for parcel-wise 0-1 range normalization ---
for k = 1:length(all_rois)
    current_roi = all_rois{k};
    fprintf('Processing ROI: %s...\n', current_roi);

    % --- Process Frauscher data for the current ROI ---
    frauscher_roi_indices = find(strcmp(Frauscher_atlases_processed, current_roi));
    if ~isempty(frauscher_roi_indices)
        frauscher_roi_data = Frauscher_data.TS_DataMat(frauscher_roi_indices, :);
        
        % Normalize Frauscher data per ROI across features (rows between 0 and 1)
%         frauscher_roi_data=normalize(frauscher_roi_data);
%         frauscher_roi_data(frauscher_roi_data==inf)=nan;
%         frauscher_roi_data(frauscher_roi_data==-inf)=nan; 
%         frauscher_roi_data_normalized = normalize(frauscher_roi_data, 2, 'range');

        % Not normalize the data
%          frauscher_roi_data_normalized = frauscher_roi_data;
        
        % Store normalized Frauscher data
        Frauscher(frauscher_roi_indices, :) = frauscher_roi_data;
        fprintf('  Processed %d Frauscher contacts in ROI: %s.\n', size(frauscher_roi_data, 1), current_roi);
    else
        fprintf('  No Frauscher data found for ROI: %s. Skipping Frauscher processing for this ROI.\n', current_roi);
    end

    % --- Process DS4100 data for the current ROI ---
    ds4100_roi_indices_in_all_data = find(strcmp(DS4100_atlases_processed, current_roi));
    if ~isempty(ds4100_roi_indices_in_all_data)

        ds4100_roi_data = DS4100_TS_DataMat_all_contacts(ds4100_roi_indices_in_all_data, :);
        
        % Normalize DS4100 data per ROI across features (rows between 0 and 1)
%         ds4100_roi_data_to_analyze=normalize(ds4100_roi_data_to_analyze);
%         ds4100_roi_data_to_analyze(ds4100_roi_data_to_analyze==inf)=nan;
%         ds4100_roi_data_to_analyze(ds4100_roi_data_to_analyze==-inf)=nan;        
%         ds4100_roi_data_normalized = normalize(ds4100_roi_data_to_analyze, 2, 'range');

        % Not normalize the data
%         ds4100_roi_data_normalized = ds4100_roi_data_to_analyze;
        
        % Store normalized DS4100 data
        DS4100(ds4100_roi_indices_in_all_data, :) = ds4100_roi_data;
        fprintf('  Processed %d DS4100 contacts in ROI: %s.\n', size(ds4100_roi_data, 1), current_roi);
    else
        fprintf('  No DS4100 data found for ROI: %s. Skipping DS4100 processing for this ROI.\n', current_roi);
    end

end

fprintf('\nNormalization per ROI across features complete for both DS4100 and Frauscher datasets.\n');

% --- Save Results ---
% Renamed filename to reflect that it contains all contacts data for both datasets
% output_analysis_filename_all_contacts = fullfile(output_folder, 'Data_AllContacts_BothDatasets.mat');
output_analysis_filename_all_contacts = fullfile(output_folder, 'Data_AllContacts_BothDatasets_epochs.mat');

fprintf('Saving normalized data (all contacts) for both datasets to %s...\n', output_analysis_filename_all_contacts);

% Save the normalized matrices for both datasets, along with relevant annotation variables
save(output_analysis_filename_all_contacts, ...
     'DS4100', ...
     'Frauscher', ...
     'DS4100_TS_DataMat_all_contacts', ... % Original data matrix for DS4100
     'Frauscher_data', ... % Original Frauscher data structure
     'DS4100_atlases_processed', ... % Atlases for all DS4100 contacts
     'Frauscher_atlases_processed', ... % Atlases for all Frauscher contacts
     'all_rois', ...
     'DS4100_annotations', ... % Save the entire DS4100 annotations structure
     'Frauscher_annotations', ... % Save the entire Frauscher annotations structure
     '-v7.3');

fprintf('Results for both datasets saved successfully.\n');