function mean_neuron_fig = plot_mean_neuron(params_file)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-05-25


%% OVERVIEW:
% This function creates a figure of the mean peri-stimulus dF/F response
% averaged across all neurons in the imaged population.

% If the data include responses to multiple stimulus or trial conditions,
% then individual activity traces will be color-coded according to the
% mapping defined either in the data file itself or in a condition settings
% file specified in the function call (see below for detail).


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
% 12) write_metadata.m
% 13) getLastCommit.m
% 14) get_sha1.m


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

%   outputDirectory - directory where output figure should be saved.


%% OUTPUTS:
% This function has no formal return, but saves to secondary storage a
% figure the includes the mean peri-stimulus dF/F response for each
% condition, averaged across all neurons in the imaged population. 


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


%% Trialize data and split by condition:

% Split all the data up by trial:
[Trials, Meta] = trialize_data(params.rawF_path, ... 
    params.galvo_path, ...
    params.timer_path, ...
    params.ardu_path, ...
    params.conditions_path, ...
    params.grab_metadata, ...
    params.pre_sec, ...
    params.post_sec, ...
    params.show_inflection_points);

% Validate that the stimulus duration is the same for every trial, and if
% not, throw a warning and skip drawing stimulus window:
all_trial_durations = [Trials.STIMDUR];
all_trial_durations_check = circshift(all_trial_durations, 1);
if isequal(all_trial_durations, all_trial_durations_check)
    equal_stimdurs = true;
end    
if equal_stimdurs
    stim_dur = all_trial_durations(1)/1000; % convert from milliseconds to seconds
else    
    stim_dur = [];
    warning('Not all stimuli have the same duration; stimulus period windows will be ommitted from plots.');
end

% Validate that all trials include measurements from same number of ROIs;
% if so, define the number of ROIs:
n_ROIs_each_trial = arrayfun(@(x) size(x.dFF, 1), Trials);
check_n_ROIs = circshift(n_ROIs_each_trial, 1);
if isequal(n_ROIs_each_trial, check_n_ROIs)
    num_ROIs = n_ROIs_each_trial(1);
else
    error('Not all trials include observations from the same number of ROIs. Please check integrity of input data');
end

Conditions = split_trials_by_condition(Trials, Conditions); % Split trials by condition


%% Compute the "mean neuron" for each condition:

Conditions = get_condition_means(Conditions); % Get mean dFF trace for each neuron for each condition:

% Average across neurons for each condition:
for c = 1:length(Conditions)
    Conditions(c).SEM = std(Conditions(c).Mean, 0, 1) / sqrt(size(Conditions(c).Mean, 1));
    Conditions(c).Mean = mean(Conditions(c).Mean, 1);
end


%% Plot:
fig_full_path = [output_directory filesep mouse '_' date '_mean_neuron_dFF_response.fig'];
fig_full_title = {'Mean neuron dF/F response by condition'; ['\fontsize{10}Mouse ' mouse]; ['\fontsize{10}Session ' date]; ['\fontsize{10}Num ROIS = ' num2str(num_ROIs)]};
mean_neuron_fig= plot_mean_peristim_trace(Conditions, pre_frames, post_frames, stim_dur, frame_rate, fig_full_title);
savefig(mean_neuron_fig, fig_full_path);


%% Write metadata:
Metadata.inputs(1).path = params.rawF_path;
Metadata.inputs(2).path = params.galvo_path;
Metadata.inputs(3).path = params.timer_path;
Metadata.inputs(4).path = params.ardu_path;
Metadata.inputs(5).path = params.conditions_path;
Metadata.inputs(6).path = params.grab_metadata;
Metadata.params.pre_onset_period = params.pre_sec;
Metadata.params.post_onset_period = params.post_sec;
Metadata.outputs(1).path = fig_full_path;

write_metadata(Metadata, [output_directory filesep 'plot_mean_neuron_metadata.json']);