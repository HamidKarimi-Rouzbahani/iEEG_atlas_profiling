clc;
close all;
clear all;

% --- Parameters ---
N_top = 40; % Number of top features to label

% --- Setup Paths ---
base_folder = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';
within_file = fullfile(base_folder, 'SubjectSpecific/WithinSubject_Final_Results.mat');
across_file = fullfile(base_folder, 'Classifier_EpileptogenicPrediction_Results_RandomForest_PatientSplit_epochs.mat');

% Path to the external Feature Names file
external_names_file = '/home/uqhkarim/Normative_HCTSA/Scripts/LittsScripts/Feature_names.mat';

files_to_plot = {within_file, across_file};
titles = {'Within-Subject Importance', 'Across-Subject Importance'};
var_names = {'all_imp_matrix', 'all_importances_matrix'}; 

for f = 1:length(files_to_plot)
    if ~exist(files_to_plot{f}, 'file')
        fprintf('Skipping %s (file not found)\n', files_to_plot{f});
        continue; 
    end
    
    % Load Data Matrix
    data = load(files_to_plot{f});
    mat = data.(var_names{f});
    [n_subs, n_feats] = size(mat);
    
    % --- FIXED: Load External Feature Names from Table ---
    if exist(external_names_file, 'file')
        ext_data = load(external_names_file);
        feat_names = cell(n_feats, 1);
        
        if isfield(ext_data, 'Feature_names')
            FT = ext_data.Feature_names;
            % Use height() for tables instead of length()
            num_names = height(FT); 
            
            for i = 1:min(n_feats, num_names)
                % Extract from CodeString column
                val = FT.Name{i}; 
                
                % Handle the "cell within a cell" if necessary
                if iscell(val)
                    feat_names{i} = char(val{1});
                else
                    feat_names{i} = char(val);
                end
            end
        end
        
        % Fill any gaps if n_feats > num_names
        for i = 1:n_feats
            if isempty(feat_names{i})
                feat_names{i} = sprintf('F%d', i);
            end
        end
    else
        fprintf('Warning: External names file not found. Using generic names.\n');
        feat_names = arrayfun(@(x) sprintf('Feat_%d', x), 1:n_feats, 'UniformOutput', false);
    end

    % --- 1. Clean Matrix & Synchronize Names ---
    nan_cols = all(isnan(mat), 1);
    mat_clean = mat(:, ~nan_cols);
    feat_names_clean = feat_names(~nan_cols);
    
    % --- 2. Normalization for Visibility ---
    v_min = nanpercentile(mat_clean(:), 1);
    v_max = nanpercentile(mat_clean(:), 98);
    mat_norm = (mat_clean - v_min) ./ (v_max - v_min);
    mat_norm(mat_norm > 1) = 1; mat_norm(mat_norm < 0) = 0;
    
    % --- 3. Calculate Marginal Traces & Find Top N ---
    avg_feat = nanmean(mat_clean, 1);
    avg_subj = nanmean(mat_clean, 2);
    
    [~, sorted_idx] = sort(avg_feat, 'descend');
    top_indices_clean = sorted_idx(1:min(N_top, length(avg_feat)));
    
    % --- 4. PLOTTING ---
    figure('Color', 'w', 'Name', titles{f}, 'Position', [50 50 1600 950]);
    
    % Main Matrix
    axMain = subplot('Position', [0.10 0.22 0.75 0.63]);
    imagesc(mat_norm);
    colormap(jet); 
    ylabel('Subjects (Index)');
    set(gca, 'XTick', []); 
    title(titles{f}, 'FontSize', 15, 'Units', 'normalized', 'Position', [0.5, 1.18, 0]);
    
    % Top Labels
    hold on;
    for i = 1:length(top_indices_clean)
        idx = top_indices_clean(i);
        name_str = strrep(feat_names_clean{idx}, '_', ' ');
        text(idx, 0.4, name_str, 'Rotation', 45, 'FontSize', 9, ...
            'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'Color', 'k', 'Interpreter', 'none');
        line([idx idx], [0.5 size(mat_norm,1)+0.5], 'Color', [1 1 1 0.25], 'LineWidth', 0.5);
    end
    
    % Vertical Trace (Left)
    axLeft = subplot('Position', [0.03 0.22 0.06 0.63]);
    plot(avg_subj, 1:length(avg_subj), 'k', 'LineWidth', 1.5);
    set(gca, 'YDir', 'reverse', 'YTickLabel', [], 'XAxisLocation', 'top');
    ylabel('Avg Imp (a.u.)'); grid on;
    ylim([0.5, size(mat_norm,1)+0.5]);
    
    % Horizontal Trace (Bottom)
    axBottom = subplot('Position', [0.10 0.08 0.75 0.10]);
    plot(avg_feat, 'Color', [0.6 0 0], 'LineWidth', 1.2);
    xlabel(sprintf('Filtered Feature Index (Total: %d)', size(mat_norm,2))); 
    ylabel('Avg (a.u.)');
    grid on; xlim([0.5, size(mat_norm,2)+0.5]);
    
    cb = colorbar(axMain, 'Position', [0.87 0.22 0.02 0.63]);
    ylabel(cb, 'Normalized Importance (a.u.)');
    
    linkaxes([axMain, axBottom], 'x');
    linkaxes([axMain, axLeft], 'y');
end

% --- Robust Percentile Helper ---
function p = nanpercentile(x, prc)
    x = x(~isnan(x)); if isempty(x), p = 0; return; end
    x = sort(x(:)); n = length(x);
    idx = (prc/100) * (n-1) + 1;
    if idx <= 1, p = x(1); elseif idx >= n, p = x(n);
    else, i = floor(idx); f = idx - i; p = x(i) + f * (x(i+1) - x(i)); end
end