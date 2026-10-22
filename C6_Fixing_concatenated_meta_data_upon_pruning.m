clc;
close all
clear all;
load('ConcatenatedChannelMetadata_more_epochs.mat');
load('HCTSA_AllPreprocessedSignals_epochs_pruned.mat','removed_indices');
all_channel_info_epoch(removed_indices)=[];
all_channel_status_description(removed_indices)=[];
all_channels_id(removed_indices)=[];
all_clinical_info_epoch(removed_indices)=[];
all_electrode_loc_info_epoch(removed_indices)=[];
clearvars removed_indices
save('ConcatenatedChannelMetadata_epochs_pruned.mat')