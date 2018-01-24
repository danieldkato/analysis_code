function condition_histograms(params_file)





%% Load parameters and metadata:

% Load params from params file:
params = loadjson(params_file);
output_directory = params.output_directory;

% Load some useful imaging session meadata:
grab_metadata = loadjson(params.grab_metadata);
mouse = grab_metadata.mouse;
date = grab_metadata.date;

% Load condition information like color codes, etc:
C_struct = loadjson(params.conditions_path);
Conditions = C_struct.conditions;
Conditions = cell2mat(Conditions); % format as an array of structs, rather than cell array

% Compute some quantities that will be needed later on:
frame_rate = grab_metadata.frame_rate;
pre_frames = ceil(frame_rate * params.pre_sec);
post_frames = ceil(frame_rate * params.post_sec);


%% For each condition, get the peak responses of all neurons:

% Split all the data up by trial:
T = trialize_data(params.rawF_path, ... 
    params.galvo_path, ...
    params.timer_path, ...
    params.ardu_path, ...
    params.conditions_path, ...
    params.grab_metadata, ...
    params.pre_sec, ...
    params.post_sec, ...
    params.show_inflection_points);

% Validate that all trials include measurements from same number of ROIs;
% if so, define the number of ROIs:
n_ROIs_each_trial = arrayfun(@(x) size(x.dFF, 1), T.Trials);
check_n_ROIs = circshift(n_ROIs_each_trial, 1);
if isequal(n_ROIs_each_trial, check_n_ROIs)
    num_ROIs = n_ROIs_each_trial(1);
else
    error('Not all trials include observations from the same number of ROIs. Please check integrity of input data');
end

Conditions = split_trials_by_condition(T.Trials, Conditions); % Split trials by condition
Conditions = get_condition_means(Conditions); % Get means for each condition

% Get the max projection of each mean dF/F traces after stimulus onset:
for c = 1:length(Conditions)
    Conditions(c).max_project = max(Conditions(c).Mean, 2);
end


%% Run an ANOVA to test if the distributions are significantly different:

% Reshape data in a format that can be handed to ANOVA function:
total_observations = length(Conditions) * num_ROIs;
y = nan(total_observations,1); % initialize data vector
groups = cell(total_observations,1); % initialize groups vector
for cc = 1:length(Conditions)
    
    % Define start and stop indices into data vector:
    start_idx = (cc-1) * num_ROIs + 1;
    end_idx = cc * num_ROIs;
    
    % Populate data vector:
    y(start_idx:end_idx) = Conditions(c).max_project(:);
    
    % Populate groups vector:
    groups(start_idx:end_idx) = repmat({Conditions(c).name}, num_ROIs, 1);
end

% Perform the ANOVA:
lm = fitlm(groups, y);
