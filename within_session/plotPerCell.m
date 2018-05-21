function [meanPaths, rawPaths] = plotPerCell(params_file)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-25


%% OVERVIEW:
% This function creates peri-stimulus activity plots for every unit recorded
% from a given acquisition session. For every unit, this function creates 2
% plots:

% 1) A plot with mean activity traces plus SEM bars, and
% 2) A plot with activity traces for every individual trial

% All plots are saved to storage in an output directory specified in the
% function call.

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
% This function saves to disk 2 plots for each ROI from the analyzed
% dataset:

% 1) A plot with mean peri-stimulus activity traces plus SEM bars for each condition, and
% 2) A plot with peri-stimulus traces for every individual trial, color-coded by condition

% In addition, this function formally returns:
% 1) meanPaths - an N x 1 cell array of strings containing full paths to
% the saved mean activity figures, and 
% 2) rawPaths - an N x 1 cell array of strings containing full paths to the
% saved individual trial activity figures. 


% TODO:
% 1) Stimulus duration is currently hard-coded into the script. This should
% also be read in dynamically from somewhere else, like trials (although
% trials currently only includes trial start frame and condition, not
% duration, so this would entail changes to how the trial matrix is
% created). In getting the trial duration from trials, one could also
% verify that all trial durations are the same (which is the only situation
% in which the shaded stimulus period rectangle makes sense).

% 2) This function requires that the condition names in trials exactly
% match the condition names in conditions. Should probably throw up an
% error if this doesn't hold. 

% 3) Maybe think about making the conditions input optional, as it's
% probably really not necessary; it could just read out the trial
% conditions from trials, then automatically assign colors and
% abbreviations to each if none are specified in conditions.

% 4) Perhaps make the outputDirectory argument optional, and have it
% default to the current working directory.


%% Load data and metadata

% Load params from params file:
params = loadjson(params_file);
output_directory = params.output_directory;

% Load grab metadata:
grab_metadata = loadjson(params.grab_metadata);
mouse = grab_metadata.mouse;
date = grab_metadata.date;

% Load condition information like color codes, etc:
C_struct = loadjson(params.conditions_path);
Conditions = C_struct.conditions;
Conditions = cell2mat(Conditions);


%% Trialize data for each neuron:

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

% Validate that all trials include measurements from same number of ROIs; if so, define the number of ROIs:
n_ROIs_each_trial = arrayfun(@(x) size(x.dFF, 1), Trials);
check_n_ROIs = circshift(n_ROIs_each_trial, 1);
if isequal(n_ROIs_each_trial, check_n_ROIs)
    num_ROIs = n_ROIs_each_trial(1);
else
    error('Not all trials include observations from the same number of ROIs. Please check integrity of input data');
end

%{
% Reorganize the data in a way that will be easy to iterate through and plot on a neuron-by-neuron basis: 
neurons_trials = neuron_trial(T); % Create a struct with the form neuron > trial
neurons_conditions_trials = neuron_condtn_trial(neurons_trials, Conditions); % Create a struct with the form neuron > condition > individual trial
neurons_conditions_means = neuron_condtn_mean(neurons_conditions_trials); % Create a struct with the form neuron > condition > mean & SEM
%}

%% Get some parameters from trialized neural data that will be useful for plotting :

frame_rate = grab_metadata.frame_rate;
pre_stim_frames = ceil(frame_rate * params.pre_sec);
post_stim_frames = ceil(frame_rate * params.post_sec);
peri_stim_frames = pre_stim_frames + post_stim_frames;

% Confirm that the stimulus duration is the same for every trial, and if
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


%% Prepare output directory and output structs to catch paths to figures: 

% Create cell arrays that will contain full paths to created figures; these will be returned to calling function:
meanPaths = cell(num_ROIs, 1);
rawPaths = cell(num_ROIs, 1);

% Make sure the output directory exists, create it if it doesn't then cd into it:
status = exist(output_directory, 'dir');
if status == 0
    mkdir(output_directory);
end
old = cd(output_directory);
mkdir('plot_each_ROI');


%% For each ROI, create 2 figures: 
% 1) a figure plotting mean traces (with SEM bars) for each condition
% 2) a figure plotting traces for all individual trials (color coded by condition)
n_figures = 0;

for n = 1:num_ROIs
    
    disp(['Plotting ROI ' num2str(n) ' out of ' num2str(num_ROIs)]);
    
    curr_n_trials = get_neuron_from_trials(Trials, n); % Get trialized data just for current neuron:
    curr_n_conditions = split_trials_by_condition(curr_n_trials, Conditions); % split trialized data for current neuron by condition
    curr_n_conditions_means = get_condition_means(curr_n_conditions);
    
    num_str = pad(num2str(n), 3, 'left', '0');
    
    % Plot mean dFF trace and SEM for each condition for the current neuron:
    mean_fig_name =  [mouse '_' date '_ROI_', num_str, '_mean_traces.fig'];
    mean_fig_full_path =  [output_directory filesep 'plot_each_ROI' filesep mean_fig_name];
    mean_fig_title = {'Mean dF/F traces by condition'; ['\fontsize{10}Mouse ' mouse]; ['\fontsize{10}Session ' date]; ['ROI #', num2str(n)]};
    mean_fig = plot_mean_peristim_trace(curr_n_conditions_means, pre_stim_frames, post_stim_frames, stim_dur, frame_rate, mean_fig_title);
    meanPaths{n} = mean_fig_full_path;
    savefig(mean_fig, mean_fig_full_path);
    n_figures = n_figures + 1;
    Metadata.outputs(n_figures).path = mean_fig_full_path;
    close(mean_fig);
    
    % Plot individual traces for each condition for the current neuron:
    individual_fig_name =  [mouse '_' date '_ROI_', num_str, '_individual_traces.fig'];
    individual_fig_full_path =  [output_directory filesep 'plot_each_ROI' filesep individual_fig_name];
    individual_fig_title = {'Individual dF/F traces by condition'; ['\fontsize{10}Mouse ' mouse]; ['\fontsize{10}Session ' date]; ['ROI #', num2str(n)]};
    single_trials_fig = plot_individual_peristim_traces(curr_n_conditions, pre_stim_frames, post_stim_frames, stim_dur, frame_rate, individual_fig_title);
    rawPaths{n} = individual_fig_full_path;
    savefig(single_trials_fig, individual_fig_full_path);
    n_figures = n_figures + 1;
    Metadata.outputs(n_figures).path = individual_fig_full_path;
    close(single_trials_fig);
    
    
end 


%% Create scatterplots:



%% Write metadata:
Metadata.inputs(1).path = params.rawF_path;
Metadata.inputs(2).path = params.galvo_path;
Metadata.inputs(3).path = params.timer_path;
Metadata.inputs(4).path = params.ardu_path;
Metadata.inputs(5).path = params.conditions_path;
Metadata.inputs(6).path = params.grab_metadata;
Metadata.params.pre_onset_period = params.pre_sec;
Metadata.params.post_onset_period = params.post_sec;

write_metadata(Metadata, [output_directory filesep 'plot_each_ROI' filesep 'plot_every_cell_metadata.json']);