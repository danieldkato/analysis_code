function F = register_trials_2_frames(galvo_path, timer_path, grab_path, show_inflection_points)

% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II.REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-13


%% I. OVERVIEW:
% This function returns an f-element vector containing the starting frame
% number for every trial delivered during a given grab, where f is the
% number of trials delivered over the coruse of the session. Any trials
% that occur before the first frame or after the last frame of the movie
% are represented by a NaN.


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


%% III. INPUTS:
% 1) galvo_path - path to a trace of the analog voltage signal used to drive
%    the slow scan-mirror galvanometer during the grab. Should be saved as a
%    LabView .dat file. Contains information about frame start times.

% 2) timer_path - path to a trace of the analog trial timer signal recorded
%    during the grab. In the current protocol, trial onset is immediately
%    preceded by a 10 ms, 5 V TTL pulse sent from the Arduino responsible for
%    controlling stimulus hardware. Should be saved as a LabView .dat file.
%    Contains information about trial start times.

% 3) grab_path - path to the raw TIFF of the movie being analyzed. The
%    directory containing the raw TIFF must also contain a JSON file called
%    '2P_metadata.json', which includes a field called 'frame_rate'
%    specifying the frame rate of the movie in frames per second.

% 4) show_inflection_points - boolean flag specifying whether or not to
%    plot galvo trace, timer trace, and frame and trial start times. 


%% IV. OUTPUTS: 
% 1) F - an f-element vector of start frames for each trial delivered
%    throughout the session, where f is the number of trials delivered
%    throughout the session. Any trials delivered before the first frame or
%    after the last frame will be represented by a NaN.


%% TODO:
% 1) Also write Trials as a .mat file in addition to JSON?

% 2) Add back in writing of metdata - I removed this because I wanted to
% switch over to writing metadata as a JSON and had to deprecate old
% writeMetadata fucntion, but now I have to replace it. 


%% Load parameters:
%{
params = loadjson(params_file);
grab_path = params.grab_path;
galvo_path = params.galvo_path;
timer_path = params.timer_path;
ardu_path = params.ardu_path; 
condition_settings = params.condition_settings;
output_path = params.output_path;
show_inflection_points = params.show_inflection_points;
%}

% Read galvo header to get sample rate:
galvo_fid = fopen(galvo_path, 'r', 'b');
[header_size, header] = SkipHeader(galvo_fid);
sample_rate = str2double(header{7}(18:end));
disp(sample_rate);

% Read image grab metadata to get framerate:
[grab_directory, nm, ext] = fileparts(grab_path);
ls = dir(grab_directory);
grab_metadata_path = ls(arrayfun(@(a) strcmp(a.name, '2P_metadata.json'), ls)); % look for a file called '2p_metadata.json' in the same directory as the raw grab file
disp(grab_metadata_path)

% Raise an error if metadata file is not found or frame rate is not defined within metadata file:
if length(grab_metadata_path) == 0
    error('Metadata file for grab not found; make sure that 2P_metadata.json is located in the same directory as raw TIFF.');
else
    disp(['metadata path name: ' grab_metadata_path(1).name]);
    grab_metadata = loadjson([grab_directory filesep grab_metadata_path(1).name]);

    % Raise an error if meta.txt does not contain a variable called frame_rate:
    if (isfield(grab_metadata,'frame_rate'))
        frame_rate = grab_metadata.frame_rate;
    else
        error('Frame rate not found; make sure that grab metadata file includes field ''frame_rate'' whose value is frame rate in Hz.');
    end
end
    
    
%% Load galvo and timer data:
galvo_signal = readContinuousDAT(galvo_path); % Load the galvanometer data from the raw .dat file into an s x 1 vector, where s is number of samples taken during grab 
galvo_signal = galvo_signal(galvo_signal>-5.0); % Need to filter out pre- and post-movie galvo signal here, otherwise LocalMinima will pick up local minima in the noise in the galvo signal outside of the actual movie

trial_timer_signal = readContinuousDAT(timer_path); % Load the trial timing data from the raw .dat file into an s x 1 vector, where s is the number of samples taken during a grab
% F.trials = read_ardulines(ardu_path, condition_settings); %% Get an ordered list of trial types from arudlines
    
    
%% Get the galvo signal sample number of every frame start:

% Compute some input parameters for LocalMinima, called below:
frame_period = 1/frame_rate;
min_distance_galvo = frame_period * sample_rate; % The function LocalMinima will include only the lowest of any local minima found within this many samples of each other
min_distance_galvo = min_distance_galvo * .9; % Fudge factor; the true number of samples between certain pairs of frame-start times is slightly less than the theoretical value
galvo_threshold = -1.6; % Whatever units gavloTrace is expressed in (Volts, I think); the function LocalMinima will exclude any local minima higher than this value; for the time being, I just got this from eyeballing a sample galvo trace, but I may ultimately need more sophisticated ways of getting this if there's any variability

% Get a vector of every galvo signal sample at which a frame begins:
frame_start_samples = LocalMinima(galvo_signal, min_distance_galvo, galvo_threshold);
disp([num2str(length(frame_start_samples)) ' frames detected']);       


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
disp([num2str(length(trial_start_samples)) ' trials detected']);    


%% Show traces if requested by user:

% Plot the local minima on top the galvo trace if desired; this can be
% handy just to double check that reasonable parameters for LocalMinima
% have been chosen, but may be cumbersome if processing large batches
% of data.

if show_inflection_points == 1
    figure;
    hold on;
    t = (1:1:length(galvo_signal));
    plot(galvo_signal);
    plot(t(frame_start_samples), galvo_signal(frame_start_samples), 'r.');
    plot(trial_timer_signal);
    plot(t(trial_start_samples), trial_timer_signal(trial_start_samples), 'r.'); 
end

    
%% Omit any trials delivered before the first frame or after the last frame: 

% Find indices of first and last trials within movie:
is_in_movie = trial_start_samples>=min(frame_start_samples) & trial_start_samples<=max(frame_start_samples); % t-element binary vector specifying whether each trial falls within the movie
movie_trials = find(is_in_movie); % q-element vector of indices into trial_start_samples correspoding to trials that fall within the movie

    
%% Match every trial to the frame within which it started:

F = NaN(size(trial_start_samples)); % initialize with NaNs; elements corresponding to trials within the movie will be populated, trials that occur before the beginning or after the end of the movie will remain NaNs

% For each trial start sample, find the maximum frame start sample below it:
for i = 1:length(movie_trials)
    curr_trial = movie_trials(i);
    [M, I] = max(frame_start_samples( frame_start_samples <= trial_start_samples(curr_trial) ));
    F(curr_trial) = I + 1; % have to add 1 because there's one frame that completes before the first local minimum
end

    
%% Add start frame to T.trials struct:
%T.trials = arrayfun(@(s,f) setfield(s,'start_frame',f),T.trials,trial_start_frames);


%{
trial_matrix = cell(length(trial_start_samples), 3);
trial_matrix(:, 1) = trial_start_frames;
disp(size(trial_start_samples));
disp(size(Trials));
trial_matrix(:, 2:3) = Trials; 
%}    
    
%% Write T to secondary storage: 

%{
json_path = [output_path filesep 'trial_info.json'];
savejson('',T,json_path);
%}

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
% Inputs:
Metadata.inputs(1).path = grab_metadata_path;
Metadata.inputs(2).path = galvo_path;
Metadata.inputs(3).path = timer_path;
Metadata.inputs(4).path = ardu_path;
Metadata.inputs(5).path = condition_settings;

% Outputs:
Metadata.inputs(1).path = json_path;

metadata_path = [output_path filesep 'trial_registration_metadata.json'];
write_metadata(Metadata, metadata_path);
%}

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