function compare_conditions_within_session(params_file)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-26


%% OVERVIEW:
% This function compares distributions of peak responses to each
% stimulus/trial condition within a single imaging session. This includes:
%
% 1) regressing peak dF/F responses against condition (equivalent to an
%    ANOVA; see OUTPUT for more detail) to test whether distributions of
%    peak responses to different conditions are significantly different.
%    This function saves the resulting LinearModel object and associated
%    F-statistic and p-value.
%
% 2) plotting histograms of peak responses to each condition side-by-side. 
%
% 3) plotting scatter plots of peak responses to selected pairs of
%    conditions


%% REQUIREMENTS:
% 1) trialize_data.m
% 2) read_ardulines.m
% 3) get_start_frames.m
% 4) readContinuousDAT.m, available at https://github.com/gpierce5/BehaviorAnalysis/blob/master/readContinuousDAT.m (commit 71b3a3c)
% 5) LocalMinima.m, available at //10.112.43.46/mnt/homes/dan/code_libraries/clay/LocalMinima.m
% 6) get_neuron_from_trials.m
% 7) split_trials_by_condition.m
% 8) get_condition_means.m
% 9) match_trials_2_conditions.m
% 10) plot_mean_peristim_trace.m
% 11) plot_individual_peristim_traces.m


%% INPUTS:
% 1) params_file - path to a JSON-formatted parameters file. The encoded
%    JSON object should include the following fields:
%
%       raw_path - path to a file containing an n x f raw activity matrix for the
%       imaging session being analyzed, where n is the number of neurons
%       identified and f is the number of frames recorded in the session. Can
%       be either a CSV or an HDF5. If the latter, the activity matrix must be
%       saved in the first dataset.
%
%       galvo_path - path to a LabView .dat file containing a trace of the
%       analog voltage signal used to drive the slow scan-mirror galvanometer
%       during the grab. Contains information about frame start times.
%
%       timer_path - path to a LabView .dat file containing a trace of the
%       analog trial timer signal recorded during the grab. In the current
%       protocol, trial onset is immediately preceded by a 50 ms, 5 V TTL
%       pulse sent from the Arduino responsible for controlling stimulus
%       hardware. Contains information about trial start times.

%       ardu_path - a path to a .txt file containing serial output from an Arduino
%       running an ArduFSM protocol for a single behavior session. In accordance
%       with the general ArduFSM framework, each line of output either
%       acknowledges the receipt of instructions from the host PC, asserts
%       upcoming trial parameters, reports recorded behavior parameters, or
%       signals the start of a trial. More information about the ArduFSM
%       framework can be found at https://github.com/cxrodgers/ArduFSM. 

%       conditions_path - path to a JSON file containing information about the
%       trial conditions to be analyzed. This should contain one top-level
%       list called "conditions", each element of which should itself be a
%       JSON object. Each of these objects should minimimally include a the
%       following fields:
%
%           name - char vector specifying a human-readable condition name
%
%           abbreviation - char vector specifying an abbreviated condition
%           name (useful for plotting)
%
%           color - RGB triple specifying the color code for the
%           corresponding condition
%
%           params - a struct with one sub-field for each parameter that
%           defines the condition; each of these sub-fields should be named
%           after the corresponding trial parameter reported in the arduino
%           serial output file.
%
%       An example element of conditions might be as follows:
%
%           {"name":"stepper only",
%           "abbreviation": "W",
%           "color": [1.00, 0.00, 0.00],
%           "params": {
%               "STPRIDX": 1,
%               "SPKRIDX": 0
%               }
%           }

%   grab_metadata_path - path to the a JSON file containing metadata file
%   about the movie being analyzed. This file should include a field called
%   'frame_rate' specifying the frame rate of the movie in frames per second.

%   pre_sec - number of seconds before stimulus onset to include in
%   peri-stimulus window.

%   post_frames - number of seconds after stimulus onset to include in
%   peri-stimulus window.

%   show_inflection_points - boolean flag specifying whether or not to
%   plot galvo trace, timer trace, and frame and trial start times. 

%   outputDirectory - directory where all created figures should be saved.


%% OUTPUTS:
% This function has no formal return, but saves to secondary storage the
% following: 

% 1) compare_conditions_within_session.mat - .mat file containing a
%    struct called compare_conditions_within_session, which includes the
%    following fields:
%
%       linear_model - a MATLAB LinearModel object containing the results
%       of fitting the following model:
%
%           response ~ -1 + C1 + C2 + C3 ...
%
%       where each element of response is the peak of one neuron's mean
%       peri-stimulus dF/F trace for one condition, and C1, C2, C3, etc.,
%       are categorical binary variables specifying the trial/stimulus
%       condition. This model does not assess interactive effects, and is
%       only to be used for assessing whether there are any significant
%       differences between the distributions of responses to different
%       trial/stimulus conditions.
%
%       F - F-statistic associated with the LinearModel object described
%       above (i.e., that associated with an ANOVA on all conditions).
%
%       p - p-value associated with the LinearModel object described
%       above (i.e., that associated with an ANOVA on all conditions).
%
% 2) compare_conditions_within_sessions.fig - figure including subplots
%    with histograms of peak responses to each stimulus/trial condition.


%% Load parameters and metadata:

% Load params from params file:
params = loadjson(params_file);
output_directory = params.output_directory;

if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

% Load some useful imaging session meadata:
grab_metadata = loadjson(params.grab_metadata);
mouse = grab_metadata.mouse;
date = grab_metadata.date;

% Load condition information like name, hardware parameters, color codes, etc:
C_struct = loadjson(params.conditions_path);
Conditions = cell2mat(C_struct.conditions);

% Compute some quantities that will be needed later on:
frame_rate = grab_metadata.frame_rate;
pre_stim_frames = ceil(frame_rate * params.pre_sec);


%% Do initial processing:

% Split data by trial:
T = trialize_data(params.rawF_path, ... 
    params.galvo_path, ...
    params.timer_path, ...
    params.ardu_path, ...
    params.conditions_path, ...
    params.grab_metadata, ...
    params.pre_sec, ...
    params.post_sec, ...
    params.show_inflection_points);

Conditions = split_trials_by_condition(T.Trials, Conditions); % Split trials by condition
Conditions = get_condition_means(Conditions); % Get the mean peristimulus dF/F trace of every neuron for every condition
Conditions = get_condition_amplitudes(Conditions, pre_stim_frames); % Get the amplitude of the mean peristimulus dF/F trace of every neuron for every condition


%% Run linear regression to test if the distributions corresponding to each condition are significantly different:

% Get the number of ROIs in the dataset:
n_ROIs_each_trial = arrayfun(@(x) size(x.dFF, 1), T.Trials);
check_n_ROIs = circshift(n_ROIs_each_trial, 1);
if isequal(n_ROIs_each_trial, check_n_ROIs)
    num_ROIs = n_ROIs_each_trial(1);
else
    error('Not all trials include observations from the same number of ROIs. Please check integrity of input data');
end

% Write response data and condition info into matrices:
total_observations = length(Conditions) * num_ROIs;
response = nan(total_observations,1); % initialize response vector
design_matrix = zeros(total_observations, length(Conditions)); % initialize design matrix
for c = 1:length(Conditions)
    
    % Define start and stop indices into data vector:
    start_idx = (c-1) * num_ROIs + 1;
    end_idx = c * num_ROIs;
    
    % Populate response vector:
    response(start_idx:end_idx) = Conditions(c).amplitudes;
    
    % Populate design matrix:
    design_matrix(start_idx:end_idx, c) = 1;
end

% Convert data matrices to table; this will make the output easier to interpret:
cond_names = arrayfun(@(x) x.name, Conditions, 'UniformOutput', false); % get condition names
cond_names = cellfun(@(x) strrep(x, ' ', '_'), cond_names, 'UniformOutput', false); % edit condition names as acceptable variable names
Tbl = array2table(design_matrix, 'VariableNames', cond_names);
Tbl.response = response; % fitlm by default assumes that last variable is response variable

% Define modelspec input argument to fitlm; this must include -1 to
% specify that the model should not include an intercept term (otherwise fitlm
% will raise a rank deficiency error):
model_spec = ['response ~ -1 + ' strjoin(cond_names, ' + ')];

% Perform the regression:
lm = fitlm(Tbl, model_spec);
[p, F, r] = coefTest(lm);

% Put the model and stats in a struct that can be saved to secondary storage:
compare_conditions_within_session.linear_model = lm;
compare_conditions_within_session.p = p;
compare_conditions_within_session.F = F;
compare_conditions_within_session.r = r;
stats_path = [output_directory filesep 'compare_conditions_within_session_stats.mat'];
save(stats_path, 'compare_conditions_within_session');


%% Plot the distributions of each condition:
fig_title = ['Mouse ' mouse ', session ' date ', dF/F responses by condition'];
f = plot_condition_histograms(Conditions, 'amplitudes', fig_title);
fig_path = [output_directory filesep 'compare_conditions_within_session_hist.fig'];
savefig(f, fig_path);


%% Create scatterplots for selected pairs of conditions:



%% Save metadata:
Metadata.inputs(1).path = params.rawF_path;
Metadata.inputs(2).path = params.galvo_path;
Metadata.inputs(3).path = params.timer_path;
Metadata.inputs(4).path = params.ardu_path;
Metadata.inputs(5).path = params.conditions_path;
Metadata.inputs(6).path = params.grab_metadata;
Metadata.params.pre_onset_period = params.pre_sec;
Metadata.params.post_onset_period = params.post_sec;
Metadata.outputs(1).path = stats_path;
Metadata.outputs(2).path = fig_path;
Metadata.outputs(3).path = h1_path;
Metadata.outputs(4).path = h2_path;
Metadata.outputs(5).path = h3_path;

write_metadata(Metadata, [output_directory filesep 'compare_conditions_within_session.json']);
