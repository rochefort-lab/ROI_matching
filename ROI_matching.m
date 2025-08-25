% Written by Theoklitos Amvrosiadis, 08/2022

% Matches ROIs across two imaging sessions with the same FOV
% Requires an initial seed of 3 pairs of matched neurons


function final_matches = ROI_matching(sesname1, sesname2, working_dir)

if nargin < 3
    working_dir = 'C:\Users\theox\Desktop\Experiments\Working';
end

% Choose FOV1 ROIset
DIR1 = dir(fullfile(working_dir, [sesname1, '\ROIs\*um\*P03.zip']));
if size(DIR1, 1) > 1
    DIR1 = DIR1(1);
end
Fov1_path = fullfile(DIR1.folder, DIR1.name);
image_path1 = [Fov1_path(1:end-3), 'tif'];
I1 = imread(image_path1);
all_my_rois1 = ReadImageJROI(Fov1_path);


% Choose FOV2 ROIset
DIR2 = dir(fullfile(working_dir, [sesname2, '\ROIs\*um\*P03.zip']));
if size(DIR2, 1) > 1
    DIR2 = DIR2(1);
end
if strcmp(sesname2, '20220719_Thy1214')
    DIR2 = dir(['C:\Users\theox\Desktop\Experiments\Working\', sesname2, '\ROIs\*um\*P10.zip']);
end

Fov2_path = fullfile(DIR2.folder, DIR2.name);
image_path2 = [Fov2_path(1:end-3), 'tif'];
I2 = imread(image_path2);
all_my_rois2 = ReadImageJROI(Fov2_path);

numrois1 = size(all_my_rois1, 2);

% Match ROIs
% Provide a seed of matches to start the process
matched_rois = get_user_matches(I1, all_my_rois1, I2, all_my_rois2, sesname1, sesname2);
if isempty(matched_rois)
    potential_matches(:, 1) = (1:numrois1)';
    potential_matches(:, 2) = nan;

    final_matches = num2cell(potential_matches);
    
    return
end

centroids1 = calculateROIcentroids(all_my_rois1);
centroids2 = calculateROIcentroids(all_my_rois2);

% Calculate all distances between ROIs
all_distances_Fov1 = squareform(pdist(centroids1));
all_distances_Fov2 = squareform(pdist(centroids2));

all_distances1 = all_distances_Fov1(matched_rois(:, 1), :);
all_distances2 = all_distances_Fov2(matched_rois(:, 2), :);

tform = estimateGeometricTransform(centroids1(matched_rois(:, 1), :), centroids2(matched_rois(:, 2), :), 'affine');
transformed_centroids1 = transformPointsForward(tform, centroids1);

transformed_centroids2 = centroids2;

all_centroid_distances = pdist2(transformed_centroids1, transformed_centroids2);

% Calculate angles between ROIs
angles1 = compute_angles(centroids1);
angles2 = compute_angles(centroids2);

angles_diff1 = angles1(matched_rois(:, 1), :);
angles_diff2 = angles2(matched_rois(:, 2), :);

alpha = 1; % Adjust this value between 0 and 1 to balance the weight of distance and angle information

M = size(matched_rois, 1); % Number of seed matches

% Calculate squared Euclidean distances between the 'feature vectors' (columns)
% The result sq_dist(i, j) = sum((all_distances1(:,i) - all_distances2(:,j)).^2)
sq_dist_distances = pdist2(all_distances1', all_distances2', 'squaredeuclidean'); % size: numrois1 x numrois2
sq_dist_angles = pdist2(angles_diff1', angles_diff2', 'squaredeuclidean');       % size: numrois1 x numrois2

% Calculate the mean squared error by dividing by the number of seed matches (M)
mse_distances = sq_dist_distances / M;
mse_angles = sq_dist_angles / M;

distance_errors = alpha * mse_distances + (1 - alpha) * mse_angles;

final_distances = all_centroid_distances + distance_errors;

for imatch = 1:size(matched_rois, 1)
    final_distances(matched_rois(imatch, 1), matched_rois(imatch, 2)) = 0;
end

% assignment = munkres(final_distances);

potential_matches = nan(numrois1, 2); % Pre-allocate with NaNs
% Assign matches based on minimum errors
potential_matches(:, 1) = (1:numrois1)';
% [~, potential_matches(:, 2)] = min(distance_errors, [], 2);
[~, potential_matches(:, 2)] = min(final_distances, [], 2);
% potential_matches(:, 2) = potential_matches(:, 2);
% potential_matches(matched_rois(:, 1), 2) = matched_rois(:, 2);


% Apply quality check and remove spurious matches
% Set a tolerance for the error above which matches should be considered
% erroneous
error_tolerance = 1e3;

for iRoi = 1:size(potential_matches, 1)
    if final_distances(potential_matches(iRoi, 1), potential_matches(iRoi, 2)) > error_tolerance ||...
            ((contains(sesname1, 'Thy1214') || contains(sesname1, 'Thy1217') || contains(sesname1, 'Bl424')) && distance_errors(potential_matches(iRoi, 1), potential_matches(iRoi, 2)) >= 50)
        potential_matches(iRoi, 2) = nan;
    end
end

locked_pairs = matched_rois; 

% Check for ROIs that have been matched more than once and keep only the
% best match
for iRoi = 1:size(potential_matches, 1)
    if ~isnan(potential_matches(iRoi, 2)) && sum(potential_matches(:, 2) == potential_matches(iRoi, 2)) > 1
        conflicting_idy = find(potential_matches(:, 2) == potential_matches(iRoi, 2));
        
        if ~logical(sum(ismember(conflicting_idy, locked_pairs(:, 1))))
            dist_err_conf = nan(1, numel(conflicting_idy));
            for iconflict = 1:numel(conflicting_idy)
                dist_err_conf(iconflict) = final_distances(conflicting_idy(iconflict), potential_matches(conflicting_idy(iconflict), 2));
            end
            roi_losers = conflicting_idy(dist_err_conf ~= min(dist_err_conf));
            potential_matches(roi_losers, 2) = nan;
        else
            roi_losers = conflicting_idy(conflicting_idy ~= conflicting_idy(ismember(conflicting_idy, locked_pairs(:, 1))));
            potential_matches(roi_losers, 2) = nan;
        end        
    end
end

% Attach a confidence value
final_matches = num2cell(potential_matches);
for iRoi = 1:size(final_matches, 1)
    if ~isnan(final_matches{iRoi, 2}) && distance_errors(final_matches{iRoi, 1}, final_matches{iRoi, 2}) < 20
        final_matches{iRoi, 3} = 'High';
    elseif ~isnan(final_matches{iRoi, 2}) && distance_errors(final_matches{iRoi, 1}, final_matches{iRoi, 2}) < 50
        final_matches{iRoi, 3} = 'Medium';
    elseif ~isnan(final_matches{iRoi, 2})
        final_matches{iRoi, 3} = 'Low';
        % final_matches{iRoi, 2} = nan;
        % final_matches{iRoi, 3} = nan;
    end
end

% Create a figure
figure
ax1 = subplot('Position', [0.001, 0.25, 0.32, 0.64]);
plotImageWithRois(I1, all_my_rois1, sesname1);
ax2 = subplot('Position', [0.321, 0.25, 0.32, 0.64]);
plotImageWithRois(I2, all_my_rois2, sesname2);
linkaxes

% Create a table
table_handle = uitable('Data', final_matches, 'ColumnName', {sesname1, sesname2, 'Certainty'}, 'RowName', [], 'Units', 'normalized', 'Position', [0.8 0.1 0.2 0.8]);
table_handle.FontSize = 12;
table_handle.ColumnEditable = [false, true, false];  % Only column 2 is editable

% Set table CellSelectionCallback

table_handle.CellEditCallback = @(src, event)roi_table_edit_callback(src, event, I1, all_my_rois1, I2, all_my_rois2, sesname1, sesname2);
table_handle.CellSelectionCallback = @(src, event)roi_table_selection_callback(src, event, I1, all_my_rois1, I2, all_my_rois2, sesname1, sesname2);

% Set the figure title
sgtitle('ROI Matching Results');

uiwait(gcf);  % Wait for user to close the figure before returning the results


function roi_table_edit_callback(src, event, I1, all_my_rois1, I2, all_my_rois2, sesname1, sesname2)
    
    selectedRowIndex = event.Indices(1);
    selectedColumnIndex = event.Indices(2);

    % Only allow editing of column 2
    if selectedColumnIndex == 2
        newValue = event.NewData;

        if ~isnan(newValue) && newValue >= 1 && newValue <= numel(all_my_rois2)
            % Find if newValue already exists in column 2
            conflictIdx = find(cell2mat(final_matches(:, 2)) == newValue);

            % Set the conflicting cell to NaN
            if ~isempty(conflictIdx)
                final_matches{conflictIdx, 2} = NaN;
                final_matches{conflictIdx, 3} = [];
            end

            % Update the edited cell
            final_matches{selectedRowIndex, 2} = newValue;
            src.Data = final_matches;
        else
            final_matches{selectedRowIndex, 2} = NaN;
            final_matches{selectedRowIndex, 3} = [];
            src.Data = final_matches;
        end
    end
    
    src.Data = final_matches;

end

function roi_table_selection_callback(src, event, I1, all_my_rois1, I2, all_my_rois2, sesname1, sesname2)
    
    if isempty(event.Indices)
        return;
    end
    selectedRowIndex = event.Indices(1);

    selectedROI1 = final_matches{selectedRowIndex, 1};
    selectedROI2 = final_matches{selectedRowIndex, 2};

    subplot(ax1);
    plotImageWithRois(I1, all_my_rois1, sesname1, selectedROI1)
    subplot(ax2);
    plotImageWithRois(I2, all_my_rois2, sesname2, selectedROI2)
    linkaxes
end

end