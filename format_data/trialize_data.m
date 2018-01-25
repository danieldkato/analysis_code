function T = trialize_data(rawF_path, galvo_path, timer_path, ardu_path, conditions_path, grab_metadata, pre_sec, post_sec, show_inflection_points)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-25


%% I. OVERVIEW: 
% This function returns a struct T containing trialized data for a given
% imaging session. See OUTPUT below for more detail.


%% II. REQUIREMENTS:
% 1) read_ardulines.m
% 2) get_start_frames.m
% 3) readContinuousDAT.m, available at https://github.com/gpierce5/BehaviorAnalysis/blob/master/readContinuousDAT.m (commit 71b3a3c)
% 4) LocalMinima.m, available at //10.112.43.46/mnt/homes/dan/code_libraries/clay/LocalMinima.m


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
%    JSON object. Each of these objects should minimimally include a the
%    following fields:
%
%       name - char vectory specifying a human-readable condition name
%
%       params - a struct with one sub-field for each parameter that
%       defines the condition; each of these sub-fields should be named
%       after the corresponding trial parameter reported in the arduino
%       serial output file.
%
%    Each element may include other fields that may be useful for plotting,
%    etc., as well. An example element of conditions might be as follows:
%
%       {"name":"stepper only",
%        "abbreviation": "W",
%        "color": [1.00, 0.00, 0.00],
%        "params": {
%           "STPRIDX": 1,
%           "SPKRIDX": 0
%           }
%       }

% 6) grab_metadata_path - path to the a JSON file containing metadata file
%    about the movie being analyzed. This file should include a field called
%    'frame_rate' specifying the frame rate of the movie in frames per second.

% 7) pre_sec - number of seconds before stimulus onset to include in
%    peri-stimulus window.

% 8) post_sec - number of seconds after stimulus onset to include in
%    peri-stimulus window.

% 9) show_inflection_points - boolean flag specifying whether or not to
%    plot galvo trace, timer trace, and frame and trial start times. 


%% IV. OUTPUTS:
% 1) T - struct containing trialized data for a session. T includes the
% following fields:

%   frame_rate - frame rate of the movie in frames per second
%
%   pre_frames - number of frames before trial start included in trial window
%
%   post_frames - number of frames after trial start included in trial window
%
%   Trials - t x 1 array of structs, where t is the number of analyzable
%   trials delivered during the movie (where analyzable means the
%   pre-stimulus period does not extend before the start of the movie and
%   the post-stimulus period does not extend beyond the end of the movie).
%   Each element of Trials has the following fields:

%       dFF - an n x p matrix containing the population neural activity
%       during the corresponding peri-stimulus period, where n is the
%       number of neurons and p is the number of frames in the
%       peri-stimulus period.
%
%       Condition - char vector giving the name of the stimulus condition.

%       In addition, each element of Trials includes one field
%       corresponding to each trial parameter, reported in the Arduino
%       serial output file. E.g., if trial t reported as having a
%       STPKIDX of 1 and a SPKRIDX of 0 in the Arduino serial output file,
%       then Trials will include:

%           Trials(t).STPRIDX = 1;
%           Trials(t).SPKRIDX = 0;


%% Load data and metadata:

% Load neural data:
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

% Get number of neurons and frames:
n_cells = size(RawF, 1);
n_frames = size(RawF, 2);
disp([num2str(n_cells) ' ROIs detected in imaging data.']);   

% Load grab metadata and compute trial window:
disp('Loading grab metadata...');
Grab = loadjson(grab_metadata);
frame_rate = Grab.frame_rate;
pre_frames = ceil(pre_sec * frame_rate);
post_frames = ceil(post_sec * frame_rate);

% Add peri-stim window params to T:
T.frame_rate = frame_rate;
T.pre_frames = pre_frames;
T.post_frames = post_frames;
disp(['Frame rate = ' num2str(frame_rate)]);
disp(['Frames before stim onset = ' num2str(pre_frames)]);
disp(['Frames after stim onset = ' num2str(post_frames)]);

% Load condition info:
disp('Loading condition definitions...');
C = loadjson(conditions_path);

% Create Conditions, a c x 1 array, where c is the number of
% conditions. Each element is a struct describing one of the conditions to
% be analyzed. Each struct includes the fields 'name', 'abbreviation', and
% 'color'. Most importantly, each struct also includes a field called
% 'params', itself a struct that in turn has one sub-field corresponding to
% each trial parameter that defines the condition:
disp('trialize_data.m');
Conditions = C.conditions;  
Conditions = cell2mat(Conditions);


%% Extract and assemble trial info:

% Get trial parameters from ardulines:
Trials = read_ardulines(ardu_path);

% Get trial start frames from galvo and timer signals:
Frames = get_start_frames(galvo_path, timer_path, grab_metadata, show_inflection_points);

% Throw some debug messages confirming whether the number of frames and trials from different sources match:
disp([num2str(n_frames) ' frames detected from imaging data.']);
disp([num2str(length(Trials)) ' trials detected from Arduino.']);
disp([num2str(length(find(~isnan(Frames)))) ' trials detected in timer signal.']);

% Make sure that the number of trials extracted from read_ardulines matches that from get_start_frames; if not, throw an error:
if length(Trials) ~= length(Frames)
    error(['Number of trials in ardulines file (' num2str(length(Trials)) ') does not match number of trials detected in timer signal (' num2str(length(Frames)) '). Please check data integrity and parameters used to detect trial starts in register_trials_2_frames.']);
end

% Add the frame start numbers to T:
for t = 1:length(Trials)
    Trials(t).start_frame = Frames(t);
end


%% Filter out trials or peri-stimulus periods that extend beyond the beginning or end of the movie:

% Omit trials whose pre-stimulus period extends before the first frame of the movie or whose post-stimulus period extends beyond the last frame of the movie:
too_close_to_start = [Trials.start_frame] < (1 + pre_frames);
too_close_to_end = [Trials.start_frame] > (n_frames - post_frames);
disp('isnan([Trials.start_frame])');
disp(class(isnan([Trials.start_frame])));
disp(isnan([Trials.start_frame]));
disp('too_close_to_start');
disp(class(too_close_to_start));
disp(too_close_to_start);
exclude = isnan([Trials.start_frame]) | (too_close_to_start | too_close_to_end);
Trials = Trials(~exclude);

% Throw some debug messages stating if any trials were ommitted:
pre_vid_nans = find(~isnan(Frames), 1, 'first') - 1;
post_vid_nans = find(~isnan(fliplr(Frames)), 1, 'first') - 1;
n_too_close_to_start = length(find(too_close_to_start));
n_too_close_to_end = length(find(too_close_to_end));

disp(['Omitting ' num2str(pre_vid_nans + n_too_close_to_start) ' trial(s) whose pre-stimulus period extends before the beginning of the movie.']);
disp(['Omitting ' num2str(post_vid_nans + n_too_close_to_end) ' trial(s) whose post-stimulus period extends beyond the end of the movie.']);


%% Add neural data to each trial:
for t = 1:length(Trials)
    
    start_frame = Trials(t).start_frame;
    
    raw = RawF(:, start_frame-pre_frames:start_frame+post_frames-1);
    baseline = mean(RawF(:, start_frame-pre_frames:start_frame-1),2);
    dFF = (raw-baseline)./baseline;
    
    Trials(t).trial_num = t;
    Trials(t).dFF = dFF;
    Trials(t).Condition = [];
end

T.Trials = Trials;
