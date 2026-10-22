clc; clear; close all;

% --- Setup Folders & Paths ---
addpath(genpath('/home/uqhkarim/Normative_HCTSA/Coordi_to_atlases'))
addpath(genpath('/home/uqhkarim/Normative_HCTSA/NIfTI_20140122'))

output_folder = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/FeatureExtraction/';
% feature_file = fullfile(output_folder, 'DS4100_Features_New.mat');
feature_file = fullfile(output_folder, 'DS4100_Features_New_per_part.mat');


% frauscher_file = fullfile(output_folder, 'Frauscher_Features_New.mat');
frauscher_file = fullfile(output_folder, 'Frauscher_Features_New_per_part.mat');

% --- Load the extracted features ---
fprintf('Loading extracted Power and Coherence features...\n');
load(feature_file, 'DS4100_Power', 'DS4100_Coherence_Mean', 'DS4100_ann', 'band_names');
load(frauscher_file, 'Frauscher_Power', 'Frauscher_Coherence_Mean', 'Frauscher_annotations');

% Combine Power and Coherence into a single feature matrix for each dataset
% Resulting Matrix: [Power_B1...B5, Coherence_B1...B5]
X_ds4100 = [DS4100_Power, DS4100_Coherence_Mean];
X_frauscher = [Frauscher_Power, Frauscher_Coherence_Mean];

num_features = size(X_ds4100, 2);
pattern = '\d+\.?\d*|%'; % Pattern to clean ROI labels

%% 1. ROI Mapping (mni2atlas)
% --- Process Frauscher ROIs ---
fprintf('Mapping Frauscher contacts to MNI Atlas...\n');
tmp_loc_f = Frauscher_annotations.ChannelPosition;
Frauscher_ROIs = cell(size(tmp_loc_f, 1), 1);
for i = 1:size(tmp_loc_f, 1)
    [atlas_info] = mni2atlas([tmp_loc_f(i,1) tmp_loc_f(i,2) tmp_loc_f(i,3)], [2], '2mm');
    if ~isempty(atlas_info.label)
        Frauscher_ROIs{i} = strtrim(regexprep(atlas_info.label{1,1}, pattern, ''));
    else
        Frauscher_ROIs{i} = 'No ROI';
    end
end

% --- Process DS4100 ROIs ---
fprintf('Mapping DS4100 contacts to MNI Atlas...\n');
tmp_loc_ds = DS4100_ann.all_electrode_loc_info_epoch;
DS4100_ROIs = cell(length(tmp_loc_ds), 1);
for i = 1:length(tmp_loc_ds)
    [atlas_info] = mni2atlas([tmp_loc_ds(i).x tmp_loc_ds(i).y tmp_loc_ds(i).z], [2], '2mm');
    if ~isempty(atlas_info.label)
        DS4100_ROIs{i} = strtrim(regexprep(atlas_info.label{1,1}, pattern, ''));
    else
        DS4100_ROIs{i} = 'No ROI';
    end
end

%% 2. Atlas-Based Z-Scoring
DS4100_Power=X_ds4100(:,1:5);
DS4100_Coherence_Mean=X_ds4100(:,6:10);

%% 3. Save Final Mapped Results
% Added Frauscher raw data and concatenated matrix to the save list
% output_filename = fullfile(output_folder, 'DS4100_ZScores_Power_Coherence.mat');
% output_filename = fullfile(output_folder, 'DS4100_ZScores_Power_Coherence_per_part.mat');
output_filename = fullfile(output_folder, 'Data_AllContacts_BothDatasets_epochs_PowCoh_per_part.mat');
save(output_filename, ...
     'DS4100_Power',...
     'DS4100_Coherence_Mean',...
     'DS4100_ROIs', ...
     'DS4100_ann', ...
     'X_frauscher', ...
     'Frauscher_Power', ...
     'Frauscher_Coherence_Mean', ...
     'Frauscher_ROIs', ...
     'band_names', ...
     '-v7.3');

fprintf('Processing complete. Final results (including Atlas data) saved to %s\n', output_filename);