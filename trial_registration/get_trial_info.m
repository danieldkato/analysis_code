function T = get_trial_info(galvo_path, timer_path, ardu_path, grab_metadata, show_inflection_points)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-15 


%% I. OVERVIEW: 
% This function creates a t-element array of structs containing start frame
% number and trial parameter information and  for each trial delivered over
% the course of a movie, where t is the number of trials delivered during
% the movie.


%% II. REQUIREMENTS:
% This function requires the following data/metadata:
% 1) A LabView .dat file of the analog voltage signal used to drive the slow scan-mirror galvanomter
% 2) A LabView .dat file of a TTL signal that goes high for ~25 ms at the beginning of each trial 
% 3) A text file containing serial feedback from an Arduino running an ArduFSM protocol (an 'ardulines' file)
% 4) A metadata file named `2P_metadata.json` that includes a `frame_rate`
%    field specifying the frame rate of the corresponding grab in frames per 
%    second. Must be located in the same directory as the raw grab data.

% This function requires the following software:
% 1) read_ardulines.m
% 2) get_start_frames.m
% 3) readContinuousDAT.m, available at https://github.com/gpierce5/BehaviorAnalysis/blob/master/readContinuousDAT.m (commit 71b3a3c)
% 4) LocalMinima.m, available at //10.112.43.46/mnt/homes/dan/code_libraries/clay/LocalMinima.m


%% III. INPUTS: 
% 1) galvo_path - path to a trace of the analog voltage signal used to drive
%    the slow scan-mirror galvanometer during the grab. Should be saved as a
%    LabView .dat file. Contains information about frame start times.

% 2) timer_path - path to a trace of the analog trial timer signal recorded
%    during the grab. In the current protocol, trial onset is immediately
%    preceded by a 10 ms, 5 V TTL pulse sent from the Arduino responsible for
%    controlling stimulus hardware. Should be saved as a LabView .dat file.
%    Contains information about trial start times.

% 3) ardu_path - a path to a .txt file containing serial output from an Arduino
%    running an ArduFSM protocol for a single behavior session. In accordance
%    with the general ArduFSM framework, each line of output either
%    acknowledges the receipt of instructions from the host PC, asserts
%    upcoming trial parameters, reports recorded behavior parameters, or
%    signals the start of a trial. More information about the ArduFSM
%    framework can be found at https://github.com/cxrodgers/ArduFSM. 

% 4) grab_path - path to the raw TIFF of the movie being analyzed. The
%    directory containing the raw TIFF must also contain a JSON file called
%    '2P_metadata.json', which includes a field called 'frame_rate'
%    specifying the frame rate of the movie in frames per second.

% 5) show_inflection_points - boolean flag specifying whether or not to
%    plot galvo trace, timer trace, and frame and trial start times. 


%% IV. OUTPUTS:
% 1) T - a t-element array of structs describing the trial parameters and
%    start frame of each trial delivered throughout the input movie, where
%    t is the number of trials delivered *during* the movie. Any trials
%    that occur before the first frame or after the last frame are omitted.

%    Each element of T minimally includes a field called 'start_frame',
%    which states the frame during which the trial begins. In addition,
%    each element of T includes one field for every trial parameter stated
%    in the input ardulines file. Thus, a trial represented in the
%    ardulines file by the following:

%    756793 TRL_104_START
%    756793 TRLP STPRIDX 1
%    756793 TRLP SPKRIDX 0
%    756793 TRLP STIMDUR 2000
%    756793 TRLP REW 0

%    would be represented in T like this:

%    T(104).frame_start = 223728
%    T(104).STPRIDX = 1
%    T(104).SPKRIDX = 0
%    T(104).STIMDUR = 2000
%    T(104).REW = 0


%% Get trial info from ardulines:
T = read_ardulines(ardu_path);


%% Get trial start frames:
F = get_start_frames(galvo_path, timer_path, grab_metadata, show_inflection_points);
disp('Trial start frames:');
disp(F);


%% Make sure that the number of trials extracted from read_ardulines matches that from register_trials_2_frames; if not, throw an error:
if length(T) ~= length(F)
    error(['Number of trials in ardulines file (' num2str(length(T)) ') does not match number of trials detected in timer signal (' num2str(length(F)) '). Please check data integrity and parameters used to detect trial starts in register_trials_2_frames.']);
end


%% Filter out any trials delivered before the first frame or after the last frame of the movie (represented by NaNs in F):
in_movie = ~isnan(F);
T = T(in_movie);
F = F(in_movie);


%% Add the frame start numbers to T:
for t = 1:length(T)
    T(t).start_frame = F(t);
end