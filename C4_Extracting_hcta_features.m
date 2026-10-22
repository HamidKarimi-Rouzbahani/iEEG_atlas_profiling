clc;
clear all;
close all;
addpath(genpath('/home/uqhkarim/Normative_HCTSA/Scripts/Toolboxes/hctsa-main'))

% Initialize a new hctsa analysis using these data and the default feature library:
TS_Init('HCTSA_AllPreprocessedSignals_epochs.mat','hctsa');

% Compute all missing values in HCTSA.mat:
TS_Compute('True');
