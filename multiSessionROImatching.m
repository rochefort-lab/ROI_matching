function commonNeuronIDs = multiSessionROImatching(sesList, varargin)
%MULTISESSIONROIMATCHING Finds common ROIs across multiple sessions using caching.
%
%   This function identifies ROIs that are present across a given list of
%   imaging sessions. It uses the first session in the list as the reference
%   and performs pairwise matching against it. A caching system is used to
%   store previously computed pairwise matches, significantly speeding up
%   repeated analyses with overlapping session lists.
%
%   SYNTAX:
%   commonNeuronIDs = multiSessionROImatching(sesList)
%   commonNeuronIDs = multiSessionROImatching(sesList, 'Name', Value, ...)
%
%   INPUTS:
%   sesList         (cell array): A cell array of session name strings. The
%                   first session (sesList{1}) is treated as the reference
%                   for all pairwise comparisons.
%
%   OPTIONAL NAME-VALUE PAIRS:
%   'working_dir'   (string): Path to the parent directory containing the
%                             session folders. Defaults to the current
%                             directory (`pwd`). The cache file is also
%                             stored here.
%
%   OUTPUT:
%   commonNeuronIDs (cell array): A cell array with the same number of
%                   elements as `sesList`.
%                   - commonNeuronIDs{1}: A vector of ROI indices from the
%                     reference session that were found in ALL other sessions
%                     in the input `sesList`.
%                   - commonNeuronIDs{i} (for i > 1): A vector of the
%                     corresponding matched ROI indices in session `sesList{i}`.
%
%   CACHE BEHAVIOR:
%   A cache file named '<ref_session_name>_roi_matching_cache.mat' is
%   created in `working_dir`. It stores the pairwise mappings between the
%   reference session and every other session it has been compared with. On
%   subsequent runs, the function loads this cache and only computes
%   mappings for new, unprocessed sessions, saving significant time.
%
%   DEPENDENCIES:
%   - ROI_matching.m: The function for pairwise ROI matching.
%
%   Written by Theoklitos Amvrosiadis

%% 1. Setup and Input Parsing

p = inputParser;
addRequired(p, 'sesList', @(x) iscell(x) && ~isempty(x));
addParameter(p, 'working_dir', pwd, @ischar);
parse(p, sesList, varargin{:});

working_dir = p.Results.working_dir;
nSessions = numel(sesList);
refSessionName = sesList{1};

%% 2. Cache Handling: Load and Validate

cacheFilename = sprintf('%s_roi_matching_cache.mat', refSessionName);
cacheFilepath = fullfile(working_dir, cacheFilename);

[mappingsMap, processedSessions] = load_and_validate_cache(cacheFilepath, refSessionName);

%% 3. Identify and Process New Sessions

% Determine which sessions from the input list are not in the cache
sessionsToProcess = setdiff(sesList, processedSessions, 'stable');
newMappingsMade = false;

if ~isempty(sessionsToProcess)
    fprintf('Found %d new session(s) to process.\n', numel(sessionsToProcess));
    newMappingsMade = true;
    
    for i = 1:numel(sessionsToProcess)
        currentSessionName = sessionsToProcess{i};
        
        % Skip if it's the reference session itself
        if strcmp(currentSessionName, refSessionName), continue; end
        
        fprintf('--> Matching: %s (ref) <--> %s\n', refSessionName, currentSessionName);
        
        try
            matches = ROI_matching(refSessionName, currentSessionName, 'working_dir', working_dir);
            mapping_i = parse_roi_matches(matches);
            
            if isempty(mapping_i)
                fprintf('    No valid matches found.\n');
            else
                fprintf('    Found %d valid matches.\n', size(mapping_i, 1));
            end
            
            % Store the new mapping and update processed list
            mappingsMap(currentSessionName) = mapping_i;
            processedSessions{end+1} = currentSessionName;
            
        catch ME
            warning('ROI_matching failed between "%s" and "%s": %s. Skipping this session.', ...
                    refSessionName, currentSessionName, ME.message);
            % Store an empty mapping to prevent reprocessing on next run
            mappingsMap(currentSessionName) = zeros(0, 2);
            processedSessions{end+1} = currentSessionName;
        end
    end
else
    fprintf('All %d requested sessions were found in the cache.\n', nSessions);
end

%% 4. Find Common ROIs for the Requested Session List

fprintf('Calculating intersection of common ROIs for the %d provided sessions...\n', nSessions);
commonNeuronIDs = cell(nSessions, 1);

% If only one session is provided, all its ROIs are "common"
if nSessions < 2
    if isKey(mappingsMap, refSessionName) % Should ideally get number of ROIs from data
        numROIsRef = size(mappingsMap(refSessionName),1); % Placeholder
        commonNeuronIDs{1} = (1:numROIsRef)';
    else
        warning('Cannot determine number of ROIs for single reference session.');
        commonNeuronIDs{1} = [];
    end
    return;
end

% Initialize the set of common ROIs with all matched ROIs from the first comparison session
firstComparisonSession = sesList{2};
if isKey(mappingsMap, firstComparisonSession)
    initialMap = mappingsMap(firstComparisonSession);
    commonRefIndices = initialMap(:, 1);
else
    warning('Mapping for session "%s" not found. Cannot find common neurons.', firstComparisonSession);
    return; % Exit, leaving commonNeuronIDs as empty cells
end

% Iteratively take the intersection with subsequent sessions
for i = 3:nSessions
    currentSessionName = sesList{i};
    if isKey(mappingsMap, currentSessionName)
        currentMap = mappingsMap(currentSessionName);
        commonRefIndices = intersect(commonRefIndices, currentMap(:, 1));
    else
        warning('Mapping for session "%s" not found. Cannot find common neurons.', currentSessionName);
        commonRefIndices = []; % This session is missing, so intersection is empty
    end
    
    % Early exit if no common neurons remain
    if isempty(commonRefIndices)
        break;
    end
end

fprintf('Found %d common ROIs across all provided sessions.\n', numel(commonRefIndices));

%% 5. Save Updated Cache

if newMappingsMade
    try
        refSessionNameCache = refSessionName; % Use a different variable name for saving
        save(cacheFilepath, 'refSessionNameCache', 'processedSessions', 'mappingsMap');
        fprintf('Cache file updated successfully: %s\n', cacheFilename);
    catch ME
        warning('Failed to save updated cache file "%s": %s', cacheFilepath, ME.message);
    end
end

%% 6. Assemble Final Output

if isempty(commonRefIndices)
    return; % Exit, leaving commonNeuronIDs as empty cells
end

commonNeuronIDs{1} = sort(commonRefIndices); % Store sorted reference indices

% For each session, find the corresponding target indices
for i = 2:nSessions
    currentSessionName = sesList{i};
    currentMap = mappingsMap(currentSessionName);
    
    % Use ismember for efficient lookup
    [isFound, locationInMap] = ismember(commonRefIndices, currentMap(:, 1));
    
    if ~all(isFound)
        warning('Internal inconsistency: a common ROI was not found in the map for "%s".', currentSessionName);
    end
    
    % Pre-allocate and populate
    targetIndices = nan(numel(commonRefIndices), 1);
    targetIndices(isFound) = currentMap(locationInMap(isFound), 2);
    commonNeuronIDs{i} = targetIndices;
end

end

%% --- Local Helper Functions ---

function [mappingsMap, processedSessions] = load_and_validate_cache(cacheFilepath, refSessionName)
    % Initializes or loads and validates the cache file.
    
    % Default state: empty cache
    mappingsMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    processedSessions = {refSessionName}; % Reference is always "processed"

    if exist(cacheFilepath, 'file')
        fprintf('Loading existing cache: %s\n', cacheFilepath);
        try
            cacheData = load(cacheFilepath);
            
            % --- Validation checks ---
            if ~isfield(cacheData, 'refSessionNameCache') || ~strcmp(cacheData.refSessionNameCache, refSessionName)
                warning('Cache reference session mismatch. Ignoring cache.');
                return;
            end
            if ~isfield(cacheData, 'mappingsMap') || ~isa(cacheData.mappingsMap, 'containers.Map')
                warning('Invalid or missing mappingsMap in cache. Ignoring cache.');
                return;
            end
            if ~isfield(cacheData, 'processedSessions') || ~iscell(cacheData.processedSessions)
                warning('Invalid or missing processedSessions in cache. Ignoring cache.');
                return;
            end
            
            % If valid, load the data
            mappingsMap = cacheData.mappingsMap;
            processedSessions = cacheData.processedSessions;
            fprintf('Cache loaded for %d previously processed sessions.\n', numel(processedSessions));
            
        catch ME
            warning('Failed to load cache file "%s": %s. Starting fresh.', cacheFilepath, ME.message);
            % Reset to default state
            mappingsMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
            processedSessions = {refSessionName};
        end
    else
        fprintf('No cache file found. A new cache will be created.\n');
    end
end

function mapping_matrix = parse_roi_matches(matches_cell)
    % Converts the cell array from ROI_matching into a clean [N x 2] numeric matrix.
    if ~iscell(matches_cell) || size(matches_cell, 2) < 2
        mapping_matrix = zeros(0, 2);
        return;
    end

    % Find rows with valid, non-NaN numeric pairs in the first two columns
    isValidRow = cellfun(@(c1, c2) isnumeric(c1) && isscalar(c1) && ...
                                  isnumeric(c2) && isscalar(c2) && ~isnan(c2), ...
                                  matches_cell(:,1), matches_cell(:,2));
    
    if ~any(isValidRow)
        mapping_matrix = zeros(0, 2);
        return;
    end
    
    % Convert valid rows to a matrix
    mapping_matrix = cell2mat(matches_cell(isValidRow, 1:2));
end
