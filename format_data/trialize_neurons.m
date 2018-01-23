function N = trialize_neurons(rawF_path, galvo_path, timer_path, ardu_path, conditions_path, grab_metadata, pre_sec, post_sec, show_inflection_points)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-15 


%% I. OVERVIEW: 
% This function returns a struct N containing trialized data for each
% neuron from a given imaging session. See OUTPUT below for more detail.


%% II. REQUIREMENTS:
% 1) get_trial_info.m
% 2) read_ardulines.m
% 3) get_start_frames.m
% 4) readContinuousDAT.m, available at https://github.com/gpierce5/BehaviorAnalysis/blob/master/readContinuousDAT.m (commit 71b3a3c)
% 5) LocalMinima.m, available at //10.112.43.46/mnt/homes/dan/code_libraries/clay/LocalMinima.m


%% III. INPUTS: 
% 1) raw_path - path to a file containing an n x f raw activity matrix for the
%    imaging session being analyzed, where n is the number of neurons
%    identified and f is the number of frames recorded in the session. Can
%    be either a CSV or an HDF5. If the latter, the activity matrix must be
%    saved in the first dataset.

% 2) galvo_path - path to a LabView .dat file containing a trace of the
%    analog voltage signal used to drive the slow scan-mirror galvanometer
%    during the grab. Contains information about frame start times.

% 3) timer_path - path to a LabView .dat file containing a trace of the
%    analog trial timer signal recorded during the grab. In the current
%    protocol, trial onset is immediately preceded by a 50 ms, 5 V TTL
%    pulse sent from the Arduino responsible for controlling stimulus
%    hardware. Contains information about trial start times.

% 4) ardu_path - a path to a .txt file containing serial output from an Arduino
%    running an ArduFSM protocol for a single behavior session. In accordance
%    with the general ArduFSM framework, each line of output either
%    acknowledges the receipt of instructions from the host PC, asserts
%    upcoming trial parameters, reports recorded behavior parameters, or
%    signals the start of a trial. More information about the ArduFSM
%    framework can be found at https://github.com/cxrodgers/ArduFSM. 

% 5) conditions_path - path to a JSON file containing information about the
%    trial conditions to be analyzed. This should contain one top-level
%    list called "conditions", each element of which should itself be a
%    JSON object. Each of these objects should minimimally include a "name"
%    field that states the name of the condition, as well as a "params"
%    field. The "params" field should itself have one sub-field for each
%    parameter that defines the condition; each of these sub-fields should
%    be named after the corresponding trial parameter reported in the
%    arduino serial output file. For an example, see

% 6) grab_metadata_path - path to the a JSON file containing metadata file
%    about the movie being analyzed. This file should include a field called
%    'frame_rate' specifying the frame rate of the movie in frames per second.

% 7) pre_frames - number of seconds before stimulus onset to include in
%    peri-stimulus window.

% 8) pre_frames - number of seconds after stimulus onset to include in
%    peri-stimulus window.

% 9) show_inflection_points - boolean flag specifying whether or not to
%    plot galvo trace, timer trace, and frame and trial start times. 


%% IV. OUTPUTS:
% 1) N - struct containing trialized data for each neuron. N includes the
% following fields:
%
%   frame_rate - frame rate of the movie in frames per second
%
%   pre_frames - number of frames before trial start included in trial window
%
%   post_frames - number of frames after trial start included in trial window
%
%   Neurons - n x 1 array of structs, where n is the number of neurons
%   identified in the movie. Each element of Neurons has the following
%   fields:
%
%       Conditions - c x 1 array of structs, where c is the number of
%       conditions being analyzed in the movie. Each element of Conditions
%       has the following fields:
%
%           Name - char array specifying human-readable condition name 
%
%           Abbreviation - char array specifying concise, human-readable
%           condition name (useful for plotting)
%       
%           Mean - q x 1 vector of the mean dF/F trace for the current
%           neuron for the current condition, where q is the number of
%           frames in the peri-stimulus window
%
%           SEM - q x 1 vector of the SEM trace for the current neuron for
%           the current condition, where q is the number of frames in the
%           peri-stimulus window
%
%           Trial - p x 1 array of structs, where p is the number of
%           trials *of the current condition* delivered during the movie.
%           Each element consists of the following fields:
%
%               dFF - q x 1 vector of the current neuron's dFF activity
%               over the current trial, where q is the length of the
%               peri-stimulus window in frames
%
%               start_frame - frame during which the stimulus was delivered
%
%           as well one field corresponding to each trial parameter
%           reported in the Arduino serial output file.


%% Load neural data:

disp('Loading neural activity data...');
[rawF_dir, rawF_name, rawF_ext] = fileparts(rawF_path);

if strcmp(rawF_ext, '.csv')
    RawF = csvread(rawF_path);
elseif strcmp(rawF_ext, '.h5')
    info = h5info(rawF_path);
    
    % Throw a warning if there's more than one dataset in the HDF5:
    if length(info.Datasets) > 1
        warning('More than one dataset detected in input HDF5; will use first one.');
    end
    
    % Throw an error if the dataset has the wrong number of dimensions:
    if length(info.Datasets(1).Dataspace.Size) ~= 2
        error('Input HDF5 dataset does not have 2 dimensions; please check double-check input file path and formatting.');
    end
    
    dataset_name = info.Datasets(1).Name;
    RawF = h5read(rawF_path, ['/' dataset_name], [1 1], [Inf Inf]);
end

n_cells = size(RawF, 1);
n_frames = size(RawF, 2);
disp([num2str(n_cells) ' ROIs detected.']);
disp([num2str(n_frames) ' detected from imaging data.']);   

%% Load grab metadata and compute trial window:
disp('Loading grab metadata...');
Grab = loadjson(grab_metadata);
frame_rate = Grab.frame_rate;
pre_frames = ceil(pre_sec * frame_rate);
post_frames = ceil(post_sec * frame_rate);

N.frame_rate = frame_rate;
N.pre_frames = pre_frames;
N.post_frames = post_frames;

disp(['Frame rate = ' num2str(frame_rate)]);
disp(['Frames before stim onset = ' num2str(pre_frames)]);
disp(['Frames after stim onset = ' num2str(post_frames)]);


%% Load condition info:
disp('Loading condition definitions...');
C = loadjson(conditions_path);

% Create a c x 1 cell array, where c is the number of conditions. Each
% element is a struct describing one of the conditions to be analyzed. Each
% struct includes the fields 'name', 'abbreviation', and 'color'. Most
% importantly, each struct also includes a field called 'params', itself a
% struct that in turn has one sub-field corresponding to each trial
% parameter that defines the condition:
Conditions = C.conditions;  


%% Get trial info:

% Create a t-element array of structs, where t is the number of trials that
% occur during the movie. Each element has a field 'start_frame' that
% specifies the start frame of the corresponding trial, as well as a field
% for each trial parameter reported in the serial output from the Arduino.
disp('Getting trial parameters and start frames...');
Trials = get_trial_info(galvo_path, timer_path, ardu_path, grab_metadata, show_inflection_points); 


%% Omit trials whose pre-stimulus period extends before the first frame of the movie or whose post-stimulus period extends beyond the last frame of the movie:
too_close_to_start = [Trials.start_frame] < (1 + pre_frames);
too_close_to_end = [Trials.start_frame] > (n_frames - post_frames);
Trials = Trials(~too_close_to_start & ~too_close_to_end);

disp(['Omitting ' num2str(length(find(too_close_to_start))) ' trial(s) whose pre-stimulus period extends before the beginning of the movie.']);
disp(['Omitting ' num2str(length(find(too_close_to_end))) ' trial(s) whose post-stimulus period extends beyond the end of the movie.']);


%% Determine which trials belong to each condition:
disp('Matching conditions to trials... ');
for c = 1:length(Conditions)
    
    filter = ones(1, length(Trials));
    
    params = fieldnames(Conditions{c}.params);  
    for p = 1:length(params) 
        param_name = params{p};
        param_value = Conditions{c}.params.(param_name);
        filter_update = [Trials.(param_name)] == param_value;
        %{
        disp(['param_name = ' param_name]);
        disp(['param value = ' num2str(param_value)]);
        disp('Trials.(param_name) :');
        disp(class([Trials.(param_name)]));
        disp([Trials.(param_name)]);
        disp(filter_update);
        %}
        filter = filter & filter_update;
    end
    Conditions{c}.Trials = find(filter);   

    % Warn the user if no trials of the current condition are found:
    if isempty(Conditions{c}.Trials)
        warning(['No trials of condition ' Conditions{c}.name ' found; condition ' Conditions{c}.name ' will be excluded from analysis.']);
    end
    
end

% Filter out any conditions that are not delivered during this movie:
condition_delivered = cellfun(@(x) ~isempty(x.Trials), Conditions);
Conditions = Conditions(condition_delivered);


%% Iterate through each neuron:

disp('Getting trialized data for each neuron...');
for n = 1:n_cells
    
    % We're going to summarize each condition for the current neuron:
    for c = 1:length(Conditions)
        
        Neurons(n).Conditions(c).Name = Conditions{c}.name;
        Neurons(n).Conditions(c).Abbreviation = Conditions{c}.abbreviation;
        
        % For each trial within the current condition, extract dF/F traces for the current neuron:
        for t = 1:length(Conditions{c}.Trials)
            
            % Get the absolute trial number:
            trial_num = Conditions{c}.Trials(t);
            start_frame = Trials(trial_num).start_frame;
            
            % Get the raw data for the trial:
            raw = RawF(n, start_frame-pre_frames:start_frame+post_frames-1);
            baseline = mean(RawF(n, start_frame-pre_frames:start_frame-1),2);
            
            % Compute the dFF trace for the trial:
            Neurons(n).Conditions(c).Trials(t).dFF = (raw-baseline)/baseline;
           
            % Add current trial parameters to N; this is redundant, but it
            % might be convenient to have the data in this format at some
            % point later on:
            curr_trial_params = fieldnames(Trials(trial_num));
            for p = 1:length(curr_trial_params)
                Neurons(n).Conditions(c).Trials(t).(curr_trial_params{p}) = Trials(trial_num).(curr_trial_params{p});
            end
            
        end
        
        % Compute mean dF/F trace for this neuron for this condition:
        all_trials = vertcat(Neurons(n).Conditions(c).Trials.dFF);
        %disp(['Neuron ' num2str(n) ', condition ' num2str(c) ' number of trials: ' num2str(size(all_trials,1))]);
        Neurons(n).Conditions(c).Mean = mean(all_trials,1);
        
        % Compute SEM trace for this neuron for this condition:
        Neurons(n).Conditions(c).SEM = std(all_trials,1) / sqrt(size(all_trials,1));
        
    end
end

N.Neurons = Neurons;