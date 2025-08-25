function matches = get_user_matches(I1, all_my_rois1, I2, all_my_rois2, sesname1, sesname2)

% Initialize variables
matches = [];
list1_selection = [];
list2_selection = [];

% Create a figure for the GUI
h = figure('Position', [50 50 1600 900], 'Name', 'Match ROIs', 'NumberTitle', 'off');

% Plot the images with ROIs
ax1 = subplot('Position', [0.001, 0.25, 0.32, 0.64]);
plotImageWithRois(I1, all_my_rois1, sesname1);
ax2 = subplot('Position', [0.321, 0.25, 0.32, 0.64]);
plotImageWithRois(I2, all_my_rois2, sesname2);


% Create the listboxes for the two ROI sets
list1 = uicontrol('Style', 'listbox', 'String', 1:size(all_my_rois1, 2), 'Position', [1125 200 250 400], 'Max', 2, 'Min', 0);
list2 = uicontrol('Style', 'listbox', 'String', 1:size(all_my_rois2, 2), 'Position', [1425 200 250 400], 'Max', 2, 'Min', 0);

% Set the listbox callback functions
set(list1, 'Callback', @list1_callback);
set(list2, 'Callback', @list2_callback);

% Create a button to add a match
add_button = uicontrol('Style', 'pushbutton', 'String', 'Add match', 'Position', [1125 100 200 60], 'Callback', @add_match);

% Create a button to confirm the matches
confirm_button = uicontrol('Style', 'pushbutton', 'String', 'Confirm', 'Position', [1425 100 200 60], 'Callback', @confirm_callback);

% Create a button to confirm the matches
no_matches_button = uicontrol('Style', 'pushbutton', 'String', 'No matches', 'Position', [1225 40 200 60], 'Callback', @nomatch_callback);

% Create a text box to display matched pairs
match_text = uicontrol('Style', 'text', 'Position', [1225 600 400 200], 'HorizontalAlignment', 'left', 'FontSize', 14);


% Wait for user input
uiwait(h);

    function list1_callback(src, ~)
        list1_selection = src.Value;
        highlight_roi()
    end

    function list2_callback(src, ~)
        list2_selection = src.Value;
        highlight_roi()
    end

    function add_match(~, ~)
        if ~isempty(list1_selection) && ~isempty(list2_selection)
            matches = [matches; list1_selection, list2_selection];
            list1_selection = [];
            list2_selection = [];
            update_match_text();
        end
    end

    function update_match_text()
        match_str = 'Matched pairs: ';
        for i = 1:size(matches, 1)
            match_str = [match_str, sprintf('ROI %d -> ROI %d\n', matches(i, 1), matches(i, 2))];
        end
        set(match_text, 'String', match_str);
    end

    function confirm_callback(~, ~)
        if size(matches, 1) >= 3
            uiresume(h);
            close(h);
        else
            msgbox('Please select at least 3 pairs of matched ROIs.', 'Not enough matches', 'error');
        end
    end

    function nomatch_callback(~, ~)
        matches = [];
        uiresume(h);
        close(h);
    end

    function highlight_roi()
        if ~isempty(list1_selection)
            axes(ax1)
            plotImageWithRois(I1, all_my_rois1, sesname1, list1_selection);
        else
            axes(ax1)
            plotImageWithRois(I1, all_my_rois1, sesname1, []);
        end
        if ~isempty(list2_selection)
            axes(ax2)
            plotImageWithRois(I2, all_my_rois2, sesname2, list2_selection);
        else
            axes(ax2)
            plotImageWithRois(I2, all_my_rois2, sesname2, []);
        end
    end


end
