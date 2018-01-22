function [meanPaths, rawPaths] = plotPerCell(params_file)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-21


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

% Note that this function itself is ONLY responsible for taking the mean
% and SEM of the data passed as input and plotting it. It dos NOT perform
% any other analysis, e.g. compute dF/F from raw fluorescence values.
% Thus, if the data passed as input are raw fluorescence values, then
% the values plotted will be raw fluorescence values. 


%% REQUIREMENTS:
% 1) trialize_data.m
% 2) read_ardulines.m
% 3) get_start_frames.m
% 4) readContinuousDAT.m, available at https://github.com/gpierce5/BehaviorAnalysis/blob/master/readContinuousDAT.m (commit 71b3a3c)
% 5) LocalMinima.m, available at //10.112.43.46/mnt/homes/dan/code_libraries/clay/LocalMinima.m
% 6) neuronize_trials.m
% 7) neuron_trials_by_condition.m
% 8) neuron_condtn_means.m
% 9) match_trials_2_conditions.m
% 10) plot_mean_peristim_trace.m
% 11) plot_individual_peristim_traces.m


%% INPUTS:
% 1) activity - path to an HDF5 file containing the data to be plotted.
% Data must be parsed as follows: the HDF5 must include one dataset for
% each trial condition or stimulus condition to be analyzed. Each dataset
% should be an N x T x P activity matrix, where N is the number of ROIs to
% be analyzed, T is the number of samples in the peri-stimulus period to be
% plotted, and P is the number of presentations of the corresponding trial
% or stimulus condition. N and and T must be the same for all datasets, but
% P may be different for each.

% In order to properly lay out figure objects like shaded rectangles for
% stimulus epochs, the HDF5 root must have the following attributes:
    
%   a) num_samples_pre_stim: the number of frames before stimulus onset
%   included in the trace for each trial

%   b) num_samples_post_stim: the number of frames after stimulus onset
%   included in the trace for each trial

% In addition to the required attributes described above, the HDF5 may also
% have the following optional attributes for specifying the appearance of
% figures:


% 2) outputDirectory - directory where all created figures should be saved.

% 3) grabMetadata - path to a MATLAB-evaluable .txt file containing
% information about the acquisition session. This must include the
% data acquisition rate in samples per second. 

% 4) conditionSettings - path to a MATLAB-evaluable .txt file defining a 
% C x 1 cell array of structs, where C is the number of distinct trial
% conditions presented during throughout the course of the trials to be
% analyzed. Each struct must have at least the following three fields:

%   a) Name - the name of the trial condition. This must exactly match the
%   trial condition descriptions in the second column of trials, described
%   above.

%   b) Color - color code for the given condition, in any valid MATLAB
%   format for encoding color. Will be used in plotting. 

%   c) abbreviation - abbreviation for the given trial condition. Will be
%   used in creating legends for each figure.


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


%% Get params from params file:
params = loadjson(params_file);
output_directory = params.output_directory;


%% Get some basic information about imaging session that will be included in figures:
session_metadata = loadjson(params.grab_metadata);
mouse = session_metadata.mouse;
date = session_metadata.date;


%% Load grab metadata:
grab_metadata = loadjson(params.grab_metadata);


%% Load condition information like color codes, etc:
C_struct = loadjson(params.conditions_path);
Conditions = C_struct.conditions;
Conditions = cell2mat(Conditions);


%% Trialize data for each neuron:

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

% Reorganize the data in a way that will be easy to iterate through and plot on a neuron-by-neuron basis: 
neurons_trials = neuron_trial(T); % Create a struct with the form neuron > trial
neurons_conditions_trials = neuron_condtn_trial(neurons_trials, Conditions); % Create a struct with the form neuron > condition > individual trial
neurons_conditions_means = neuron_condtn_mean(neurons_conditions_trials); % Create a struct with the form neuron > condition > mean & SEM


%% Get some parameters from trialized neural data that will be useful for plotting :

frame_rate = neurons_trials.frame_rate;
pre_stim_frames = neurons_trials.pre_frames;
post_stim_frames = neurons_trials.post_frames;
peri_stim_frames = pre_stim_frames + post_stim_frames;

%{
disp(['pre_stim_frames = ' num2str(pre_stim_frames)]);
disp(['post_stim_frames = ' num2str(post_stim_frames)]);
%}

% Confirm that the stimulus duration is the same for every trial, and if
% not, throw a warning and skip drawing stimulus window:
all_trial_durations = [T.Trials.STIMDUR];
all_trial_durations_check = circshift(all_trial_durations, 1);
if isequal(all_trial_durations, all_trial_durations_check)
    equal_stimdurs = true;
end    
if equal_stimdurs
    stim_dur = all_trial_durations(1)/1000; % convert from milliseconds to seconds
else    
    stim_dir = [];
    warning('Not all stimuli have the same duration; stimulus period windows will be ommitted from plots.');
end


%% Prepare output directory and output structs to catch paths to figures: 

% Create cell arrays that will contain full paths to created figures; these will be returned to calling function:
num_ROIs = length(neurons_trials.Neurons);
meanPaths = cell(num_ROIs, 1);
rawPaths = cell(num_ROIs, 1);

% Make sure the output directory exists, create it if it doesn't then cd into it:
status = exist(output_directory, 'dir');
if status == 0
    mkdir(output_directory);
end
old = cd(output_directory);


%% For each ROI, create 2 figures: 
% 1) a figure plotting mean traces (with SEM bars) for each condition
% 2) a figure plotting traces for all individual trials (color coded by condition)
n_figures = 1;

for n = 1:num_ROIs
    
    num_str = pad(num2str(n), 3, 'left', '0');
    
    % Plot mean dFF trace and SEM for each condition for the current neuron:
    mean_data = neurons_conditions_means.Neurons(n);
    mean_fig_name =  [mouse '_' date '_ROI_', num_str, '_mean_traces.fig'];
    mean_fig_full_path =  fullfile(output_directory, mean_fig_name);
    mean_fig_title = {'Mean dF/F traces by condition'; ['\fontsize{10}Mouse ' mouse]; ['\fontsize{10}Session ' date]; ['ROI #', num2str(n)]};
    plot_mean_peristim_trace(mean_data, frame_rate, pre_stim_frames, post_stim_frames, Conditions, mean_fig_full_path, stim_dur, mean_fig_title);
    meanPaths{n} = mean_fig_full_path;
    n_figures = n_figures + 1;
    Metadata.outputs(n_figures).path = mean_fig_full_path;
    
    % Plot individual traces for each condition for the current neuron:
    individual_data = neurons_conditions_trials.Neurons(n);
    individual_fig_name =  [mouse '_' date '_ROI_', num_str, '_individual_traces.fig'];
    individual_fig_full_path =  fullfile(output_directory, individual_fig_name);
    individual_fig_title = {'Individual dF/F traces by condition'; ['\fontsize{10}Mouse ' mouse]; ['\fontsize{10}Session ' date]; ['ROI #', num2str(n)]};
    plot_individual_peristim_traces(individual_data, frame_rate, pre_stim_frames, post_stim_frames, Conditions, individual_fig_full_path, stim_dur, individual_fig_title);
    rawPaths{n} = individual_fig_full_path;
    n_figures = n_figures + 1;
    Metadata.outputs(n_figures).path = individual_fig_full_path;
end 


%% Write metadata:
Metadata.inputs(1).path = params.rawF_path;
Metadata.inputs(2).path = params.galvo_path;
Metadata.inputs(3).path = params.timer_path;
Metadata.inputs(4).path = params.ardu_path;
Metadata.inputs(5).path = params.conditions_path;
Metadata.inputs(6).path = params.grab_metadata;
Metadata.params.pre_onset_period = params.pre_sec;
Metadata.params.post_onset_period = params.post_sec;

write_metadata(Metadata, 'plot_every_cell_metadata.json');