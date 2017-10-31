function [T] = register_trials_2_frames(params_file)

% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II.REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2017-10-30


%% I. OVERVIEW:
% This function returns a T x 1 array of structs containing the starting
% frame number and condition for every trial delivered during a given grab,
% where T is the number of trials delivered during the grab.


%% II. REQUIREMENTS:

% This function requires the following data/metadata:
% 1) A LabView .dat file of the analog voltage signal used to drive the slow scan-mirror galvanomter
% 2) A LabView .dat file of an analog trial timer signal recorded during the grab 
% 3) A text file containing Arduino feedback received over a serial port
% 4) A metadata file named `2p_metadata.json` that includes a `frame_rate`
%    field specifying the frame rate of the corresponding grab in frames per 
%    second. Must be located in the same directory as the raw grab data.

% This function requires the following software:
% 1) The MATLAB function readContinuousDAT, available at https://github.com/gpierce5/BehaviorAnalysis/blob/master/readContinuousDAT.m (commit 71b3a3c)
% 2) The MATLAB function LocalMinima, available at //10.112.43.46/mnt/homes/dan/code_libraries/clay/LocalMinima.m
% 3) The MATLAB function read_ardulines, available at https://github.com/danieldkato/trial_registration/blob/master/read_ardulines.m
% 4) The MATLAB toolbox JSONlab, available at https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files


%% III. INPUTS:
% 1) params_file - path to a JSON file defining trial registration
%    parameters. This file must minimally include the following fields: 

%   a) galvo_path - path to a trace of the analog voltage signal used to drive
%   the slow scan-mirror galvanometer during the grab. Should be saved as a
%   LabView .dat file. Contains information about frame start times.

%   b) timer_path - path to a trace of the analog trial timer signal recorded
%   during the grab. In the current protocol, trial onset is immediately
%   preceded by a 10 ms, 5 V TTL pulse sent from the Arduino responsible for
%   controlling stimulus hardware. Should be saved as a LabView .dat file.
%   Contains information about trial start times.

%   c) ardu_path - path to a .txt file containing serial output received from
%   an Arduino over the course of the grab. Contains information about trial type.

%   d) grab_path - string argument containing path to the raw TIFF of
%   the grab being analyzed. This directory MUST contain a metadata JSON
%   file called `2p_metadata.json` that includes the field `frame_rate`, 
%   which gives the grab frame rate in frames per second.

%   e) output_path - optional path to the directory where the output matrix should be
%   saved.

%   f) show_inflection_points - optional boolean argument controlling whether
%   or not to plot the galvo and timer traces along with identified frame and
%   trial start times. 


%% IV. OUTPUTS: 
% This function returns a T x 1 struct array containing the starting frame
% number and parameters for every trial delivered during the grab, where T
% is the number of complete trials delivered during the grab. This
% parameters information includes trial duration and a concise,
% human-readable condition name, as well as other, protocol-specific
% parameters.


%% TODO:
% 1) Also write Trials as a .mat file in addition to JSON?

% 2) Add back in writing of metdata - I removed this because I wanted to
% switch over to writing metadata as a JSON and had to deprecate old
% writeMetadata fucntion, but now I have to replace it. 


%% Load parameters:
params = loadjson(params_file);
grab_path = params.galvo_path;
timer_path = params.timer_path;
ardu_path = params.ardu_path;
condition_settings = params.condition_settings;
output_path = params.output_path;
show_inflection_points = params.show_inflection_points;

% Read image grab metadata to get framerate:
[grab_directory, nm, ext] = fileparts(grab_path);
list = dir(grab_directory);
grab_metadata_path = list(arrayfun(@(a) strcmp(a.name, '2p_metadata.json'), list)); % look for a file called '2p_metadata.json' in the same directory as the raw grab file

% Raise an error if meta.txt is not found:
if length(grab_metadata_path) == 0
    error('Metadata file for grab not found; make sure that 2p_metadata.json is located in the same directory as raw TIFF.');
else
    disp('metadata path name:');
    disp(grab_metadata_path(1).name);
    grab_metadata_fid = fopen(fullfile(grab_directory, grab_metadata_path(1).name), 'r');
    disp('grab_metadata_fid');
    disp(grab_metadata_fid);
    grab_metadata = loadjson(grab_metadata_path);

    % Raise an error if meta.txt does not contain a variable called frame_rate:
    if (isfield(grab_metadata,'frame_rate'))
        frame_rate = grab_metadata.frame_rate;
        T.frame_rate = frame_rate;
    else
        error('Frame rate not found; make sure that grab metadata file includes field ''frame_rate'' whose value is frame rate in Hz.');
    end

end

% Read galvo header to get sample rate:
galvo_fid = fopen(galvo_path, 'r', 'b');
[header_size, header] = SkipHeader(galvo_fid);
sample_rate = str2double(header{7}(18:end));
    

%% Set output display parameters:
if nargin< 6
    show_inflection_points = 0;
end

if nargin < 5
    output_path = cd;
end
    
    
%% Load galvo, timer and Arduino data:
galvo_signal = readContinuousDAT(galvo_path); % Load the galvanometer data from the raw .dat file into an s x 1 vector, where s is number of samples taken during grab 
trial_timer_signal = readContinuousDAT(timer_path); % Load the trial timer data from the raw .dat file into an s x 1 vector, where s is the number of samples taken during a grab
T.trials = read_ardulines(ardu_path, condition_settings); %% Get an ordered list of trial types from arudlines
    
    
%% Get the galvo signal sample number of every frame start:

% Compute some input parameters for LocalMinima, called below:
frame_period = 1/frame_rate;
min_distance_galvo = frame_period * sample_rate; % The function LocalMinima will include only the lowest of any local minima found within this many samples of each other
min_distance_galvo = min_distance_galvo * .9; % Fudge factor; the true number of samples between certain pairs of frame-start times is slightly less than the theoretical value
galvo_threshold = -1.6; % Whatever units gavloTrace is expressed in (Volts, I think); the function LocalMinima will exclude any local minima higher than this value; for the time being, I just got this from eyeballing a sample galvo trace, but I may ultimately need more sophisticated ways of getting this if there's any variability

% Get a vector of every galvo signal sample at which a frame begins:
frame_start_samples = LocalMinima(galvo_signal, min_distance_galvo, galvo_threshold);
    
    
%% Get the trial timer signal sample number of every trial start:

% Get every sample index in timerTrace corresponding to the onset of a
% new trial; trial onsets are indicated by local maxima, so run
% LocalMinima on -trialTrace:

% Compute some parameters for LocalMinima, called below:
minITI = 3; % Seconds; again, it would be better if there were a way to do this dynamically
min_distance_timer = minITI * sample_rate;
timer_threshold = -4; % timerTrace units (Volts, I think); I just eyeballed this for now, but I should probably find a way to get this dynamically.

% Get a vector of every timer signal sample at which a trial begins: 
trial_start_samples = LocalMinima(-trial_timer_signal, min_distance_timer, timer_threshold);
    
    
%% Show traces if requested by user:

% Plot the local minima on top the galvo trace if desired; this can be
% handy just to double check that reasonable parameters for LocalMinima
% have been chosen, but may be cumbersome if processing large batches
% of data.

%{
if show_inflection_points == 1
    figure;
    hold on;
    t = (1:1:length(galvo_signal));
    plot(galvo_signal);
    plot(t(frame_start_samples), galvo_signal(frame_start_samples), 'r.');
    plot(trial_timer_signal);
    plot(t(trial_start_samples), trial_timer_signal(trial_start_samples), 'r.'); 
end
%}    
    
%% Omit any trials delivered before the first frame or after the last frame: 

% Find indices of first and last trials within movie:
first_trial_in_movie = find(trial_start_samples>=min(frame_start_samples),'first');
last_trial_in_movie = find(trial_start_samples<max(frame_start_samples),'last');

% Omit trials that fall outside of the movie:
trial_start_samples = trial_start_samples( trial_start_samples>=min(frame_start_samples) & trial_start_samples<=max(frame_start_samples) );
T.trials = T.trials(first_trial_in_movie:last_trial_in_movie); % Also omit these trials from struct array Trials:

    
%% Match every trial to the frame within which it started:

trial_start_frames = NaN(length(T),1);

% For each trial start sample, find the maximum frame start sample below it:
for i = 1:length(trial_start_frames)
    [M, I] = max(frame_start_samples( frame_start_samples <= trial_start_samples(i) ));
    trial_start_frames(i) = I + 1; % have to add 1 because there's one frame that completes before the first local minimum
end

    
%% Add start frame to T.trials struct:
T.trials = arrayfun(@(s,f) setfield(s,'start_frame',f),T.trials,trial_start_frames);


%{
trial_matrix = cell(length(trial_start_samples), 3);
trial_matrix(:, 1) = trial_start_frames;
disp(size(trial_start_samples));
disp(size(Trials));
trial_matrix(:, 2:3) = Trials; 
%}    
    
%% Write T to secondary storage: 

savejson('',T,'trial_info.json');

%{
% Check that the output path exists:
status = exist(output_path);
if status == 0
    mkdir(output_path);
end

% Create and open the output file for writing:
filename = fullfile(output_path, 'trialMatrix.csv');
fileID = fopen(filename, 'w');

% Write header:
if exist('grabPath', 'var')
    fprintf(fileID, strcat(['Grab,', strrep(grab_path,'\','\\'), '\n']));
else
    fprintf(fileID, strcat(['Grab, none specified \n']));
end
fprintf(fileID, strcat(['Galvanometer trace, ', strrep(galvo_path,'\','\\'), '\n']));
fprintf(fileID, strcat(['Trial timer signal trace, ', strrep(timer_path,'\','\\'), '\n']));
fprintf(fileID, strcat(['Arduino output, ', strrep(ardu_path,'\','\\'), '\n']));
fprintf(fileID, '\n');
fprintf(fileID, 'Trial start frame number, Trial type, Trial duration (ms) \n');

% Write body:
formatSpec = '%d, %s, %d \n';
for i = 1:size(trial_matrix, 1)
    fprintf(fileID, formatSpec, trial_matrix{i,:});
end

fclose(fileID);
%}    
    
%% Write metadata:

%{
inputs = {{'galvanometer trace', galvo_path};
          {'timer signal trace', timer_path};
          {'arduino feedback', ardu_path}
    };

outputs = {{'trial matrix', strcat([output_path, 'triaMatrix.csv'])}};
parameters = {};

old = cd(output_path);
writeMetadata('trial_registration', 'sampleDomain', inputs, outputs, parameters);
cd(old);
%}
    
end