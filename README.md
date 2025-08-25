# ROI Matching for Longitudinal 2-Photon Imaging

This MATLAB repository provides a function to match Regions of Interest (ROIs) from neuronal imaging (e.g., 2-photon calcium imaging) across two different sessions recorded from the same Field of View (FOV).

## Algorithm Overview

The matching process is semi-automated and uses a hybrid strategy to achieve robust results even with minor FOV shifts or rotations between sessions.

1.  **Manual Seeding**: The user first manually identifies a small number (typically 3-5) of obvious "anchor" neuron pairs.
2.  **Affine Transformation**: These seed pairs are used to calculate an affine geometric transformation (`estimateGeometricTransform`) that globally aligns the ROI centroids from session 1 onto the coordinate space of session 2.
3.  **Cost Matrix Calculation**: A cost matrix is then generated to score every possible pairing of an ROI from session 1 with an ROI from session 2. The cost for a pair $(i, j)$ is a weighted sum of two metrics:
    * **Spatial Distance**: The Euclidean distance between the transformed centroid of ROI $i$ and the centroid of ROI $j$.
    * **Feature Distance**: A measure of how geometrically similar the local environment of ROI $i$ is to that of ROI $j$. This is calculated by creating a "feature vector" for each ROI, consisting of its distances to all other ROIs in its own session. The mean squared error between the feature vectors of ROI $i$ and ROI $j$ gives the feature distance.
4.  **Initial Matching & Filtering**: Potential matches are found by selecting the pairs with the minimum cost. Matches are then filtered to remove pairs with excessively high costs and to resolve conflicts where multiple ROIs from session 1 are matched to a single ROI in session 2 (the lowest-cost pairing is kept).
5.  **Interactive GUI**: The results are presented in a graphical user interface showing the two FOVs side-by-side with an editable table of matches. The user can click through pairs, visually inspect them, and manually correct any errors before closing the figure to output the final list of matches.

---

## Dependencies

* **MATLAB** (R2020b or newer recommended)
* **Image Processing Toolbox**
* **Statistics and Machine Learning Toolbox**

### Helper Functions

The main function `ROI_matching.m` relies on several helper functions which must be located in the `helpers/` directory and included in the MATLAB path.

* `ReadImageJROI.m`: Reads `.zip` files containing ROI sets exported from ImageJ. (Note: This is a common function, often from the SiMView or NoRMCorre packages, or custom lab code).
* `get_user_matches.m`: A GUI function for the user to select the initial seed matches.
* `calculateROIcentroids.m`: Computes the geometric center of each ROI.
* `plotImageWithRois.m`: Displays the mean fluorescence image with ROI outlines.
* `compute_angles.m`: (Optional) Computes inter-ROI angles if used in the cost metric.

---

## Installation and Usage

1.  Clone the repository:
    ```bash
    git clone [https://github.com/your-username/ROI-Matching.git](https://github.com/your-username/ROI-Matching.git)
    ```
2.  Open MATLAB.
3.  Add the repository and its `helpers` subfolder to your MATLAB path:
    ```matlab
    addpath(genpath('/path/to/ROI-Matching'));
    ```
4.  Run the main function with your session data.

### Example

```matlab
% Define session names and the working directory
ses1 = '20220815_Mouse1_Ses1';
ses2 = '20220816_Mouse1_Ses2';
work_dir = '/path/to/your/experiment/folder';

% Run the matching function
final_matches = ROI_matching(ses1, ses2, 'working_dir', work_dir);

% Run with custom parameters
final_matches_custom = ROI_matching(ses1, ses2, ...
    'working_dir', work_dir, ...
    'alpha', 0.6, ...       % Give slightly more weight to spatial distance
    'error_tol', 1500);     % Be more lenient with the error threshold
```
The function will launch an interactive window. After you review and close the window, the `final_matches` variable will be populated in your workspace.

---

## Multi-Session Matching (`multiSessionROImatching.m`)

This function extends the pairwise matching to find neurons that are consistently tracked across multiple sessions. It uses the `ROI_matching.m` function as its core engine and implements a caching system to avoid re-computing existing matches, making it highly efficient for longitudinal analysis.

### Workflow and Caching

1.  **Reference Session**: The first session in the input list is designated as the "reference". All other sessions are matched against this reference.
2.  **Cache Check**: The function looks for a `.mat` cache file (e.g., `ref_session_name_roi_matching_cache.mat`). If a valid cache exists, it loads the stored pairwise mappings.
3.  **Process New Sessions**: It compares the list of requested sessions against the list of processed sessions from the cache. Only new, unprocessed sessions are run through the `ROI_matching.m` function.
4.  **Update Cache**: The results of any new pairings are added to the mapping data structure, and the cache file is overwritten with the updated information.
5.  **Calculate Intersection**: The function then uses the complete set of required mappings (from cache or newly computed) to find the set of reference ROIs that have a valid match in *every single one* of the requested sessions.
6.  **Format Output**: The final output is a cell array providing the list of common ROI indices for the reference session, and the corresponding indices for each of the other sessions.

### Usage Example

```matlab
% Define a list of sessions for analysis
all_sessions = {'20220815_Mouse1_Ses1', ...  % This will be the reference
                '20220816_Mouse1_Ses2', ...
                '20220817_Mouse1_Ses3', ...
                '20220818_Mouse1_Ses4'};

% Define the working directory
work_dir = '/path/to/your/experiment/folder';

% --- First run ---
% This will compute matches for Ses2, Ses3, and Ses4 against Ses1
% and create a new cache file.
common_neurons = multiSessionROImatching(all_sessions, 'working_dir', work_dir);

% --- Subsequent run with a subset of sessions ---
% This run will be much faster as it will load all mappings from the cache.
% It will only compute the new intersection for the requested subset.
subset_sessions = {'20220815_Mouse1_Ses1', ...
                   '20220817_Mouse1_Ses3'};
common_subset = multiSessionROImatching(subset_sessions, 'working_dir', work_dir);

% --- Subsequent run with a new session added ---
% This will load the existing cache and only compute the new mapping for Ses5.
sessions_with_new = {'20220815_Mouse1_Ses1', ...
                     '20220816_Mouse1_Ses2', ...
                     '20220819_Mouse1_Ses5'}; % New session
common_new = multiSessionROImatching(sessions_with_new, 'working_dir', work_dir);
```



## File Structure

The repository is organized as follows:

```
ROI-Matching/
├── ROI_matching.m             # Main function
├── helpers/                   # Directory for all dependency functions
└── README.md                  # This documentation
```

It is recommended to create a separate directory for your experimental data.
