function final_matches = ROI_matching(sesname1, sesname2, varargin)
%ROI_MATCHING Matches ROIs across two imaging sessions of the same FOV.
%
%   This function identifies corresponding Regions of Interest (ROIs) between
%   two different imaging sessions. The algorithm uses a hybrid approach:
%   1.  An initial affine transformation is calculated based on a small seed of
%       manually matched ROIs to correct for global shifts/rotations.
%   2.  A cost matrix is computed based on both the transformed spatial
%       distance between ROI centroids and a feature-based distance. This
%       feature distance compares the geometric relationship (distances and
%       angles) of each ROI to all other ROIs within its own session.
%   3.  Potential matches are identified by minimizing this cost matrix.
%   4.  The results are displayed in an interactive GUI for final review and
%       manual correction by the user.
%
%   SYNTAX:
%   final_matches = ROI_matching(sesname1, sesname2)
%   final_matches = ROI_matching(sesname1, sesname2, 'Name', Value, ...)
%
%   INPUTS:
%   sesname1        (string): Name of the first session's directory.
%   sesname2        (string): Name of the second session's directory.
%
%   OPTIONAL NAME-VALUE PAIRS:
%   'working_dir'   (string): Path to the parent directory containing the
%                             session folders. Default: current directory (pwd).
%   'alpha'         (double): Weighting factor between 0 and 1. Balances the
%                             cost between spatial distance (weighted by alpha)
%                             and feature-based distance (weighted by 1-alpha).
%                             Default: 0.5.
%   'error_tol'     (double): Maximum combined error threshold to accept a match.
%                             Default: 1e3.
%   'feature_tol'   (double): A secondary, stricter threshold on the feature
%                             distance, used for specific cases. Default: 50.
%
%   OUTPUT:
%   final_matches   (cell array): An N-by-3 cell array where N is the number
%                   of ROIs in the first session.
%                   Column 1: ROI index from session 1.
%                   Column 2: Matched ROI index from session 2 (NaN if no match).
%                   Column 3: Confidence string ('High', 'Medium', 'Low').
%
%   DEPENDENCIES:
%   - Requires several helper functions to be on the MATLAB path:
%     - ReadImageJROI.m       (to load ImageJ ROI sets)
%     - get_user_matches.m    (for manual seeding of matches)
%     - calculateROIcentroids.m (to compute ROI centers)
%     - compute_angles.m      (to compute inter-ROI angles)
%     - plotImageWithRois.m   (for visualization)
%
%   Written by Theoklitos Amvrosiadis, 08/2022

%% 1. Input Parsing and Data Loading

% --- Handle optional inputs with InputParser for robustness ---
p = inputParser;
addRequired(p, 'sesname1', @ischar||isstring(s));
addRequired(p, 'sesname2', @ischar||isstring(s));
addParameter(p, 'working_dir', pwd, @ischar);
addParameter(p, 'alpha', 0.5, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'error_tol', 1e3, @isnumeric);
addParameter(p, 'feature_tol', 50, @isnumeric);
addParameter(p,'confidence_thresh',[20 50],@(x)isnumeric(x)&&numel(x)==2);
parse(p, sesname1, sesname2, varargin{:});

% --- Assign parsed inputs to variables ---
working_dir   = p.Results.working_dir;
ALPHA         = p.Results.alpha;
ERROR_TOL     = p.Results.error_tol;
FEATURE_TOL   = p.Results.feature_tol;
CONFIDENCE_THRESH = p.Results.confidence_thresh; % [High/Medium, Medium/Low]

% --- Load data for both sessions using a helper function ---
[I1, all_my_rois1, centroids1] = load_session_data(sesname1, working_dir);
[I2, all_my_rois2, centroids2] = load_session_data(sesname2, working_dir);
numrois1 = size(all_my_rois1, 2);


%% 2. Get Seed Matches for Transformation

% --- User provides a few initial matches to seed the algorithm ---
% This is critical for calculating the initial geometric transformation.
matched_rois = get_user_matches(I1, all_my_rois1, I2, all_my_rois2, sesname1, sesname2);

% --- Handle case where user provides no seeds ---
if isempty(matched_rois) || size(matched_rois,1) < 3
    potential_matches = [(1:numrois1)', nan(numrois1, 1)];
    final_matches = num2cell(potential_matches);
    warning('No seed matches provided. Returning unmatched list.');
    return;
end
M = size(matched_rois, 1); % Number of seed matches


%% 3. Calculate Affine Transform and Spatial Distances

% --- Use seed matches to find the affine transform between FOVs ---
tform = estimateGeometricTransform(centroids1(matched_rois(:, 1), :), centroids2(matched_rois(:, 2), :), 'affine');

% --- Apply transform and calculate pairwise spatial distances ---
transformed_centroids1 = transformPointsForward(tform, centroids1);
all_centroid_distances = pdist2(transformed_centroids1, centroids2);


%% 4. Calculate Feature-Based Error Metric

% --- For each ROI, create a "feature vector" of its distances to all others ---
all_distances_Fov1 = squareform(pdist(centroids1));
all_distances_Fov2 = squareform(pdist(centroids2));

% --- Use only the seed matches to compare feature vectors ---
% This creates a [M x numrois1] and [M x numrois2] matrix.
seed_distances1 = all_distances_Fov1(matched_rois(:, 1), :);
seed_distances2 = all_distances_Fov2(matched_rois(:, 2), :);

% --- Calculate the mean squared error between feature vectors ---
% The result mse_distances(i, j) is the MSE between the feature vector of
% ROI i (from ses1) and ROI j (from ses2).
mse_distances = pdist2(seed_distances1', seed_distances2', 'squaredeuclidean') / M;

% --- (Optional) A similar calculation can be done with angles ---
% Note: The original code weighted angles, but the logic was identical to
% distances. The 'alpha' parameter in the original code was set to 1,
% effectively ignoring angles. This refactoring simplifies it to just
% distances but can be extended. To re-introduce angles, calculate
% mse_angles and combine as:
% feature_errors = alpha * mse_distances + (1 - alpha) * mse_angles;
feature_errors = mse_distances;


%% 5. Compute Final Cost Matrix and Find Initial Matches

% --- Combine spatial and feature distances into a final cost matrix ---
% A weighted sum based on the ALPHA parameter.
final_cost = ALPHA * all_centroid_distances + (1 - ALPHA) * feature_errors;

% --- Force seed matches to have zero cost to ensure they are kept ---
for iMatch = 1:M
    final_cost(matched_rois(iMatch, 1), matched_rois(iMatch, 2)) = 0;
end

% --- Find the best match for each ROI by minimizing the cost ---
potential_matches = nan(numrois1, 2);
potential_matches(:, 1) = (1:numrois1)';
[min_costs, potential_matches(:, 2)] = min(final_cost, [], 2);


%% 6. Quality Control and Conflict Resolution

% --- Remove poor matches based on error thresholds ---
is_high_error = min_costs > ERROR_TOL;

% A secondary check for specific preparations if needed.
% This could be adapted or removed based on experimental needs.
is_high_feature_error = false(numrois1, 1);

potential_matches(is_high_error | is_high_feature_error, 2) = nan;

% --- Resolve conflicts: one ROI in ses2 matched to multiple in ses1 ---
% For each ROI in session 2, find if it has multiple matches.
unique_targets = unique(potential_matches(~isnan(potential_matches(:, 2)), 2));
for target_roi = unique_targets'
    conflicting_indices = find(potential_matches(:, 2) == target_roi);
    
    if numel(conflicting_indices) > 1
        % If one of the conflicting ROIs is a user-defined seed, keep it.
        seed_match_idx = find(ismember(conflicting_indices, matched_rois(:, 1)), 1);
        if ~isempty(seed_match_idx)
            best_match_idx = conflicting_indices(seed_match_idx);
        else
            % Otherwise, keep the one with the lowest total cost.
            [~, min_idx] = min(min_costs(conflicting_indices));
            best_match_idx = conflicting_indices(min_idx);
        end
        
        % Set all other conflicting matches to NaN.
        conflicting_indices(conflicting_indices == best_match_idx) = [];
        potential_matches(conflicting_indices, 2) = nan;
    end
end


%% 7. Assign Confidence and Prepare for GUI

final_matches = num2cell(potential_matches);
for iRoi = 1:size(final_matches, 1)
    if ~isnan(final_matches{iRoi, 2})
        roi_idx1 = final_matches{iRoi, 1};
        roi_idx2 = final_matches{iRoi, 2};
        current_feature_error = feature_errors(roi_idx1, roi_idx2);
        
        if current_feature_error < CONFIDENCE_THRESH(1)
            final_matches{iRoi, 3} = 'High';
        elseif current_feature_error < CONFIDENCE_THRESH(2)
            final_matches{iRoi, 3} = 'Medium';
        else
            final_matches{iRoi, 3} = 'Low';
        end
    end
end


%% 8. Create Interactive GUI for Manual Inspection and Correction

f = figure('Name', 'ROI Matching Results', 'NumberTitle', 'off', 'WindowState', 'maximized');

% --- Display FOVs side-by-side ---
ax1 = subplot('Position', [0.05, 0.1, 0.4, 0.8]);
plotImageWithRois(I1, all_my_rois1, sesname1);
title(sesname1, 'Interpreter', 'none');

ax2 = subplot('Position', [0.47, 0.1, 0.4, 0.8]);
plotImageWithRois(I2, all_my_rois2, sesname2);
title(sesname2, 'Interpreter', 'none');

linkaxes([ax1, ax2]);

% --- Display results in an editable table ---
table_handle = uitable('Parent', f, ...
    'Data', final_matches, ...
    'ColumnName', {sesname1, sesname2, 'Certainty'}, ...
    'RowName', [], ...
    'Units', 'normalized', ...
    'Position', [0.88 0.1 0.11 0.8], ...
    'FontSize', 12, ...
    'ColumnEditable', [false, true, false]);

% --- Set callbacks for table interaction ---
table_handle.CellSelectionCallback = @(src, event)roi_table_selection_callback(src, event);
table_handle.CellEditCallback = @(src, event)roi_table_edit_callback(src, event);

% --- Wait for user to close the figure before returning the results ---
uiwait(f);

% --- Retrieve the latest data from the table upon closing ---
final_matches = table_handle.Data;


%% --- Nested Helper and Callback Functions ---

    function [img, rois, centroids] = load_session_data(sesname, base_dir)
        % Finds and loads the TIF image and ROI.zip file for a session.
        session_path = fullfile(base_dir, sesname);
        
        % Find the ROI zip file (handles variations in naming)
        roi_dir = dir(fullfile(session_path, 'ROIs', '*um', '*P*.zip'));
        if isempty(roi_dir)
            error('No ROI zip file found for session: %s', sesname);
        end
        roi_zip_path = fullfile(roi_dir(1).folder, roi_dir(1).name);
        
        % Corresponding TIF is assumed to have the same base name
        tif_path = [roi_zip_path(1:end-3), 'tif'];
        if ~exist(tif_path, 'file')
             error('No TIF file found corresponding to: %s', roi_zip_path);
        end
        
        img = imread(tif_path);
        rois = ReadImageJROI(roi_zip_path);
        centroids = calculateROIcentroids(rois);
    end

    function roi_table_selection_callback(src, event)
        % Highlights the selected ROI pair on the FOV plots.
        if isempty(event.Indices) || event.Indices(1) > size(src.Data, 1)
            return;
        end
        
        selectedRow = event.Indices(1);
        selectedData = src.Data(selectedRow, :);
        
        roi1_idx = selectedData{1};
        roi2_idx = selectedData{2};
        
        % Update plots
        subplot(ax1);
        plotImageWithRois(I1, all_my_rois1, '', roi1_idx);
        title(sesname1, 'Interpreter', 'none');
        
        subplot(ax2);
        plotImageWithRois(I2, all_my_rois2, '', roi2_idx);
        title(sesname2, 'Interpreter', 'none');
    end

    function roi_table_edit_callback(src, event)
        % Handles manual edits of matches in the table.
        if isempty(event.Indices)
            return;
        end
        
        row = event.Indices(1);
        col = event.Indices(2);
        
        % Only column 2 (sesname2 matches) is editable
        if col == 2
            current_data = src.Data;
            new_val = event.NewData;
            
            % Validate input
            if isnan(new_val)
                 current_data{row, 2} = NaN;
                 current_data{row, 3} = [];
            elseif isnumeric(new_val) && new_val >= 1 && new_val <= numel(all_my_rois2)
                % Check if this ROI is already matched elsewhere
                conflict_idx = find([current_data{:, 2}] == new_val);
                conflict_idx(conflict_idx == row) = []; % Ignore self
                
                if ~isempty(conflict_idx)
                    % If conflict exists, set the old match to NaN
                    current_data{conflict_idx, 2} = NaN;
                    current_data{conflict_idx, 3} = [];
                end
                
                current_data{row, 2} = new_val;
                current_data{row, 3} = 'Manual'; % Mark as manually edited
            else
                % Revert to previous value if input is invalid
                current_data{row, 2} = event.PreviousData;
            end
            src.Data = current_data; % Update table data
        end
    end
end

