clc; close all; clearvars;

% --- 1. Setup Folders & Load Data ---
base_folder = '/home/uqhkarim/Normative_HCTSA/Scripts/LittDataset/';
input_folder = fullfile(base_folder, 'FeatureExtraction/');
input_data_filename = fullfile(input_folder, 'DS4100_Power_Coh_Normalized_Final_per_part.mat');

fprintf('Loading data...\n');
load(input_data_filename, 'DS4100_ann');

% --- 2. Extract Patient-Level Clinical Table ---
if isfield(DS4100_ann, 'all_clinical_info_epoch')
    clin_info_all = DS4100_ann.all_clinical_info_epoch;
else
    clin_info_all = DS4100_ann.all_clinical_info; 
end

patient_ids_all = erase(string({clin_info_all.participant_id}), 'sub-')';
[~, first_idx] = unique(patient_ids_all, 'stable');
clin_subset = clin_info_all(first_idx);

vars_to_keep = {'outcome', 'lesion_status', 'sex', 'target', 'therapy', ...
                'hand', 'implant', 'age', 'age_onset', 'engel'};
T = struct2table(clin_subset);
T = T(:, vars_to_keep);

% --- 3. Feature Engineering & Cleaning ---
T.duration = double(T.age) - double(T.age_onset);
T.age = []; T.age_onset = []; 

for i = 1:width(T)
    var_name = T.Properties.VariableNames{i};
    if strcmpi(var_name, 'duration')
        T.(var_name) = double(T.(var_name));
        T.(var_name)(isnan(T.(var_name))) = nanmedian(T.(var_name));
    else
        col_str = string(T{:,i});
        col_str(col_str == "" | ismissing(col_str) | strcmpi(col_str, "n/a") | ...
                strcmpi(col_str, "unknown") | col_str == "NaN") = "Unknown";
        T.(var_name) = categorical(col_str);
    end
end

% --- 4. Filtering: Drop conditions with < 2 samples ---
[counts, categories] = groupcounts(T.outcome);
valid_cats = categories(counts >= 2);
T = T(ismember(T.outcome, valid_cats), :);
T.outcome = removecats(T.outcome);

% --- 5. Correlation Matrix with Values & Asterisks ---
T_numeric = zeros(height(T), width(T));
for i = 1:width(T)
    if iscategorical(T{:,i})
        T_numeric(:,i) = double(grp2idx(T{:,i}));
    else
        T_numeric(:,i) = T{:,i}; 
    end
end

[rho, pval] = corr(T_numeric, 'Type', 'Spearman', 'Rows', 'pairwise');
abs_rho = abs(rho);
mask = tril(ones(size(rho)));
plot_rho = abs_rho; plot_rho(mask == 0) = NaN;

figure('Color', 'w', 'Name', 'Correlation Matrix', 'Position', [150 150 850 750]);
imagesc(plot_rho, 'AlphaData', ~isnan(plot_rho)); 
cb = colorbar; ylabel(cb, 'Absolute Spearman Rho | \rho |');
colormap(flipud(hot(256))); caxis([0 1]); axis square;

var_labels = strrep(T.Properties.VariableNames, '_', ' ');
set(gca, 'XTick', 1:width(T), 'XTickLabel', var_labels, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:width(T), 'YTickLabel', var_labels);
title(sprintf('%d \\times %d Clinical Correlations (Lower Triangular, * = p < 0.05)', width(T), width(T)));

for i = 1:size(rho,1)
    for j = 1:i
        sig_text = '';
        if pval(i,j) < 0.05 && i ~= j, sig_text = '*'; end
        text(j, i, sprintf('%.2f%s', abs_rho(i,j), sig_text), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
    end
end

% --- 6. Significant Pair Deep-Dive ---
fprintf('\n--- Significant Clinical Pairings (p < 0.05) & Value Overlaps ---\n');
[row_sig, col_sig] = find(pval < 0.05 & tril(ones(size(pval)), -1));
for k = 1:length(row_sig)
    v1 = T.Properties.VariableNames{row_sig(k)};
    v2 = T.Properties.VariableNames{col_sig(k)};
    fprintf('\nSignificant: %s vs %s (|r|=%.2f)\n', v1, v2, abs_rho(row_sig(k), col_sig(k)));
    if iscategorical(T.(v1)) && iscategorical(T.(v2))
        tbl = table(T.(v1), T.(v2), 'VariableNames', {v1, v2});
        summary_pair = groupsummary(tbl, {v1, v2});
        [~, s_idx] = sort(summary_pair.GroupCount, 'descend');
        for n = 1:min(3, height(summary_pair))
             fprintf('   - %d patients: [%s: %s] & [%s: %s]\n', ...
                summary_pair.GroupCount(s_idx(n)), v1, string(summary_pair{s_idx(n),1}), ...
                v2, string(summary_pair{s_idx(n),2}));
        end
    end
end