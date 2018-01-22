function F = get_start_frames(galvo_path, timer_path, grab_metadata_path, show_inflection_points)
% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II.REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-21


%% I. OVERVIEW:
% This function returns an f-element vector containing the starting frame
% number for every trial delivered during a given grab, where f is the
% number of trials delivered over the coruse of the session. Any trials
% that occur before the first frame or after the last frame of the movie
% are represented by a NaN.


%% II. REQUIREMENTS:
% 1) readContinuousDAT.m, available at https://github.com/gpierce5/BehaviorAnalysis/blob/master/readContinuousDAT.m (commit 71b3a3c)
% 2) LocalMinima.m, available at //10.112.43.46/mnt/homes/dan/code_libraries/clay/LocalMinima.m


%% III. INPUTS:
% 1) galvo_path - path to a trace of the analog voltage signal used to drive
%    the slow scan-mirror galvanometer during the grab. Should be saved as a
%    LabView .dat file. Contains information about frame start times.

% 2) timer_path - path to a trace of the analog trial timer signal recorded
%    during the grab. In the current protocol, trial onset is immediately
%    preceded by a 10 ms, 5 V TTL pulse sent from the Arduino responsible for
%    controlling stimulus hardware. Should be saved as a LabView .dat file.
%    Contains information about trial start times.

% 3) grab_metadata_path - path to a JSON file containing metadata related
%    to the movie being analyzed. This must include a field called
%    'frame_rate', which sepcifies the frame rate of the movie in frames
%    per second.

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

% Read galvo header to get sample rate:
galvo_fid = fopen(galvo_path, 'r', 'b');
[header_size, header] = SkipHeader(galvo_fid);
sample_rate = str2double(header{7}(18:end));
disp(sample_rate);

% Read grab metadata:
grab_metadata = loadjson(grab_metadata_path);

% Raise an error if meta.txt does not contain a variable called frame_rate:
if (isfield(grab_metadata,'frame_rate'))
    frame_rate = grab_metadata.frame_rate;
else
    error('Frame rate not found; make sure that grab metadata file includes field ''frame_rate'' whose value is frame rate in Hz.');
end
    
    
%% Load galvo data:
galvo_signal = readContinuousDAT(galvo_path); % Load the galvanometer data from the raw .dat file into an s x 1 vector, where s is number of samples taken during grab 

% Need to filter out pre- and post-movie galvo signal here, otherwise LocalMinima will pick up local minima in the noise in the galvo signal outside of the actual movie
useable_galvo = galvo_signal>-5.0;
galvo_signal = galvo_signal(useable_galvo);


%% Load timer data:
trial_timer_signal = readContinuousDAT(timer_path); % Load the trial timing data from the raw .dat file into an s x 1 vector, where s is the number of samples taken during a grab
trial_timer_signal = trial_timer_signal(useable_galvo); % need to exclude samples already excluded from galvo_signal so that galvo and timer signals have the same offset

    
%% Get the sample number of every trial start:

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


%% Get the sample number of every frame start:

% Compute some input parameters for LocalMinima, called below:
frame_period = 1/frame_rate;
min_distance_galvo = frame_period * sample_rate; % The function LocalMinima will include only the lowest of any local minima found within this many samples of each other
min_distance_galvo = min_distance_galvo * .9; % Fudge factor; the true number of samples between certain pairs of frame-start times is slightly less than the theoretical value
galvo_threshold = -1.6; % Whatever units gavloTrace is expressed in (Volts, I think); the function LocalMinima will exclude any local minima higher than this value; for the time being, I just got this from eyeballing a sample galvo trace, but I may ultimately need more sophisticated ways of getting this if there's any variability

% Get a vector of every galvo signal sample at which a frame begins:
frame_start_samples = LocalMinima(galvo_signal, min_distance_galvo, galvo_threshold);
disp([num2str(length(frame_start_samples)) ' frames detected from galvo trace.']);       


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
    plot(t(trial_start_samples), trial_timer_signal(trial_start_samples), 'g.'); 
end

    
%% Omit any trials delivered before the first frame or after the last frame: 

% Find indices of first and last trials within movie:
after_first_frame = trial_start_samples>min(frame_start_samples);
before_last_frame = trial_start_samples<max(frame_start_samples); 
is_in_movie = after_first_frame & before_last_frame; % t-element binary vector specifying whether each trial falls within the movie
movie_trials = find(is_in_movie); % q-element vector of indices into trial_start_samples correspoding to trials that fall within the movie

%{
n_trials_too_early = length(find(~after_first_frame));
n_trials_too_late = length(find(~before_last_frame));

if n_trials_too_early > 0
    disp(['Omitting ' num2str(n_trials_too_early) ' trials that occur before first frame.']);    
end

if n_trials_too_late > 0
    disp(['Omitting ' num2str(n_trials_too_late) ' trials that occur after last frame.'])    
end
%}

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