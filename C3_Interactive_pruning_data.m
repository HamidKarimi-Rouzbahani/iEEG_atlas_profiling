% MATLAB code to prune HCTSA data interactively

clc;
clear;
close all;

% --- Configuration ---
% input_filename = 'HCTSA_AllPreprocessedSignals_1epoch.mat';
% output_filename = ['HCTSA_AllPreprocessedSignals_1epoch_pruned.mat'];
input_filename = 'HCTSA_AllPreprocessedSignals_epochs.mat';
output_filename = ['HCTSA_AllPreprocessedSignals_epochs_pruned.mat'];
signals_per_display = 5; % Display 5 signals at a time (as requested)

% --- Load the data ---
if ~exist(input_filename, 'file')
    fprintf(2, 'Error: Input file ''%s'' not found in the current directory.\n', input_filename);
    fprintf(2, 'Please ensure the file is in the same directory as this script or provide the full path.\n');
    return;
end

fprintf('Loading data from %s...\n', input_filename);
data = load(input_filename);

% The user has consistently confirmed that signals can be plotted using plot(timeSeriesData{1,4}).
% This unequivocally indicates that 'timeSeriesData' is a cell array,
% and the numeric signal data is directly contained within each cell.
% 'keywords' and 'labels' are described as '1*N structure', so we'll assume
% they are struct arrays (e.g., labels(i).label).

keywords = data.keywords;
labels = data.labels;
timeSeriesData = data.timeSeriesData; % This should be a 1xN cell array based on user feedback

N = length(timeSeriesData); % Get the number of signals
fprintf('Loaded %d time series signals.\n', N);

% --- DIAGNOSTIC: Inspecting 'labels' variable after loading ---
% This section provides critical information about the 'labels' variable's structure.
fprintf('\n--- DEBUG: Inspecting ''labels'' variable ---\n');
fprintf('    Class of labels: %s\n', class(labels));
fprintf('    Size of labels: [%s]\n', mat2str(size(labels)));
if isstruct(labels)
    fprintf('    ''labels'' is a struct array. Fieldnames: {%s}\n', strjoin(fieldnames(labels), ', '));
    if N > 0 && ~isempty(labels)
        fprintf('    Content of first labels element (labels(1)):\n');
        % Using evalc to capture the output of disp for a cleaner print to console
        fprintf('%s', evalc('disp(labels(1))')); 
    else
        fprintf('    labels array is empty or has no elements.\n');
    end
elseif iscell(labels)
    fprintf('    ''labels'' is a cell array.\n');
    if N > 0 && ~isempty(labels)
        fprintf('    Content of first cell (labels{1}): ');
        if ischar(labels{1}) || (isstring(labels{1}) && isscalar(labels{1}))
            fprintf('''%s''\n', labels{1});
        else
            fprintf('%s (Class: %s)\n', mat2str(labels{1}), class(labels{1}));
        end
    else
        fprintf('    labels cell array is empty or has no elements.\n');
    end
else
    fprintf('    ''labels'' is neither a struct array nor a cell array.\n');
end
fprintf('-------------------------------------------\n');


removed_indices = []; % Initialize an empty array to store indices of removed signals

% --- Interactive Pruning Loop ---
current_idx = 1;

while current_idx <= N
    % Determine the range of signals to display in the current batch
    end_idx_display = min(current_idx + signals_per_display - 1, N);
    
    fprintf('\n--- Displaying signals from %d to %d (out of %d) ---\n', current_idx, end_idx_display, N);
    
    % Create a new figure for plotting and set it to be docked
    figure('Name', sprintf('Signals %d-%d', current_idx, end_idx_display), 'NumberTitle', 'off');
    set(gcf, 'WindowStyle', 'docked'); % Make the figure dock into MATLAB
    
    sgtitle(sprintf('Signals %d-%d (total %d)', current_idx, end_idx_display, N)); % Super title for the figure
    
    % Extract the current batch of signals. This slice preserves the type (cell array).
    signals_in_batch = timeSeriesData(current_idx:end_idx_display);
    
    for i = 1:length(signals_in_batch)
        subplot(signals_per_display, 1, i); % Arrange plots in a single column (as requested)
        
        current_signal_data = [];
        current_signal_label = 'No label found'; % Default label if none is successfully retrieved

        % --- Data Access for Signal Plotting ---
        % Based on user's direct plotting example: plot(timeSeriesData{1,4}),
        % assuming timeSeriesData is a cell array where each cell directly contains the numeric data.
        if iscell(signals_in_batch) % Confirm it's a cell array batch
            current_signal_data = signals_in_batch{i}; % Access the numeric data directly from the cell
        else 
            % Fallback for unexpected timeSeriesData type (e.g., struct array with a .data field)
            current_struct = signals_in_batch(i);
            if isfield(current_struct, 'data') && ~isempty(current_struct.data)
                current_signal_data = current_struct.data;
            end
        end
        
        % --- Label Access ---
        label_global_idx = current_idx + i - 1;
        if label_global_idx <= length(labels)
            if isstruct(labels) && isfield(labels(label_global_idx), 'label') && ~isempty(labels(label_global_idx).label)
                current_signal_label = labels(label_global_idx).label;
            elseif iscell(labels) && ~isempty(labels{label_global_idx}) % If labels is a cell array
                % Assuming if it's a cell array, the cell directly contains the string label
                if ischar(labels{label_global_idx}) || (isstring(labels{label_global_idx}) && isscalar(labels{label_global_idx}))
                    current_signal_label = labels{label_global_idx};
                else % If cell content is not a char array or scalar string
                    current_signal_label = sprintf('Non-string label (Class: %s)', class(labels{label_global_idx}));
                end
            else
                % Debug print if label access fails for a specific element
                fprintf(2, '    DEBUG: Label for global index %d could not be retrieved.\n', label_global_idx);
                if isstruct(labels)
                    fprintf(2, '        Labels is struct, but field "label" missing or empty for this element.\n');
                    if label_global_idx <= length(labels) && isfield(labels(label_global_idx), 'label')
                         fprintf(2, '        Content of labels(%d).label (if exists): %s\n', label_global_idx, mat2str(labels(label_global_idx).label));
                    end
                    if label_global_idx <= length(labels)
                        fprintf(2, '        Available fields for labels(%d): {%s}\n', label_global_idx, strjoin(fieldnames(labels(label_global_idx)), ', '));
                    end
                elseif iscell(labels)
                    fprintf(2, '        Labels is cell, but cell %d is empty or not string.\n', label_global_idx);
                    if label_global_idx <= length(labels) && ~isempty(labels{label_global_idx})
                        fprintf(2, '        Content of labels{%d}: %s (Class: %s)\n', label_global_idx, mat2str(labels{label_global_idx}), class(labels{label_global_idx}));
                    end
                else
                    fprintf(2, '        Labels variable is neither a struct nor a cell array, or is empty.\n');
                end
            end
        else
            fprintf(2, '    DEBUG: Label index %d out of bounds for labels (length %d).\n', label_global_idx, length(labels));
        end

        % Squeeze the data to remove singleton dimensions (e.g., 1x1x12000 becomes 1x12000 or 12000x1)
        current_signal_data_to_plot = squeeze(current_signal_data);

        % --- Plotting Logic ---
        % Check all conditions for valid plotting data *after* squeezing.
        % 'isvector' check is removed here as plot() handles 2D numeric arrays (e.g., NxM) by plotting columns.
        plot_condition_met = isnumeric(current_signal_data_to_plot) && ...
                             ~isempty(current_signal_data_to_plot) && ...
                             all(isfinite(current_signal_data_to_plot(:))); % Use (:) to flatten for all/any checks

        if plot_condition_met
            plot(current_signal_data_to_plot); % Plot the data
            title_str = sprintf('Row %d: %s', i, current_signal_label);
            title(title_str, 'Interpreter', 'none', 'FontSize', 8);
        else
            plot(NaN); % Plot an empty graph if data is unsuitable
            title_str = sprintf('Row %d: Data unsuitable or missing', i);
            
            % --- CRITICAL DIAGNOSTIC OUTPUT for debugging empty plots (printed to command window) ---
            fprintf(2, '    DEBUG: Signal %d (global index %d) data issue for plotting.\n', i, current_idx + i - 1);
            fprintf(2, '        Data before squeeze: Size: [%s], Class: %s\n', mat2str(size(current_signal_data)), class(current_signal_data));
            fprintf(2, '        Data AFTER squeeze: Size: [%s], Class: %s\n', mat2str(size(current_signal_data_to_plot)), class(current_signal_data_to_plot));
            fprintf(2, '        Conditions (after squeeze): isnumeric=%d, ~isempty=%d, all(isfinite)=%d\n', ...
                isnumeric(current_signal_data_to_plot), ~isempty(current_signal_data_to_plot), all(isfinite(current_signal_data_to_plot(:))));
            if ~isempty(current_signal_data_to_plot)
                fprintf(2, '        NaNs: %d, Infs: %d\n', sum(isnan(current_signal_data_to_plot(:))), sum(isinf(current_signal_data_to_plot(:))));
                if isnumeric(current_signal_data_to_plot) && ~isvector(current_signal_data_to_plot) && ~isscalar(current_signal_data_to_plot)
                    fprintf(2, '        NOTE: Data is numeric but not a 1D vector after squeeze (it''s a 2D matrix). Plotting columns.\n');
                end
            end
            % --- END CRITICAL DIAGNOSTIC OUTPUT ---

            title(title_str, 'Interpreter', 'none', 'FontSize', 8, 'Color', 'red');
        end
        grid on;
    end
    
    drawnow; % Update the figure immediately
    
    % Get user input
    user_input = input('Enter row numbers to remove (e.g., 1 3 5) from CURRENT batch, or press Enter to keep all and continue: ', 's');
    
    % Parse user input
    if isempty(user_input)
        % User pressed Enter, keep all signals in this batch
        fprintf('Keeping all signals in this batch.\n');
        current_idx = end_idx_display + 1; % Move to the next batch
    else
        % Convert input string to an array of numbers
        rows_to_remove_in_batch = str2num(user_input); %#ok<ST2NM>
        
        % Validate input
        if isempty(rows_to_remove_in_batch) || any(rows_to_remove_in_batch <= 0) || any(rows_to_remove_in_batch > length(signals_in_batch))
            fprintf(2, 'Invalid input. Please enter valid row numbers within the current batch (1-%d).\n', length(signals_in_batch));
            % Do not close figure here, let user see it again before next input.
            continue; % Re-display the same batch for valid input
        end
        
        % Calculate global indices for removal
        global_indices_to_remove = (current_idx - 1) + rows_to_remove_in_batch;
        removed_indices = [removed_indices, global_indices_to_remove]; %#ok<AGROW>
        
        fprintf('Marked signals for removal (global indices): %s\n', mat2str(global_indices_to_remove));
        
        % Close the current figure (it's docked so it won't pop out anyway)
        close(gcf);
        
        % Increment current_idx to skip the processed batch and move to the next
        current_idx = end_idx_display + 1;
    end
    
    % Close the figure after processing each batch, unless input was invalid and loop continued
    if isvalid(gcf) % Check if figure handle is still valid before closing
        close(gcf);
    end
end

% --- Prune the data ---
% Ensure removed_indices are unique and sorted
removed_indices = unique(removed_indices);

% Store original total count for reference
N_original = N; 

if isempty(removed_indices)
    fprintf('\nNo signals were marked for removal. The output file will be a copy of the input.\n');
    keywords = keywords;
    labels = labels;
    timeSeriesData = timeSeriesData;
else
    fprintf('\nRemoving %d signals from a total of %d...\n', length(removed_indices), N_original);
    
    % Create logical index for signals to keep
    keep_indices = true(1, N_original);
    keep_indices(removed_indices) = false;
    
    % Perform pruning. This works for both struct arrays and cell arrays when indexing with logical arrays.
    keywords = keywords(keep_indices);
    labels = labels(keep_indices);
    timeSeriesData = timeSeriesData(keep_indices);
    
    fprintf('Remaining signals after pruning: %d\n', length(timeSeriesData));
end

% --- Save the pruned data and removed indices ---
fprintf('Saving pruned data to %s...\n', output_filename);
% Use -v7.3 for potentially large files, especially with 12000-sample signals.
save(output_filename, 'keywords', 'labels', 'timeSeriesData', 'removed_indices', '-v7.3'); 
fprintf('Pruning complete. Data saved to %s.\n', output_filename);

% Display a summary of removed indices
if ~isempty(removed_indices)
    fprintf('\nRemoved signal global indices: %s\n', mat2str(removed_indices));
else
    fprintf('\nNo signals were removed.\n');
end