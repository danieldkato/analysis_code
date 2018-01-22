function N = neuron_trial(T)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-19


%% I. OVERVIEW: 
% This function returns a struct N containing trialized data for a given
% imaging session split up by neuron. It essentially takes the output T
% from the function trialize_data() and converts it from a trial-major
% structure to a neuron-major structure. See OUTPUT below for more detail.


%% II. REQUIREMENTS:
% 1) MATLAB v >= ???


%% III. INPUTS: 
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


%% IV. OUTPUTS:
% 1) N - struct containing trialized data for a session, split up by
% neuron. N includes the following fields:

%   frame_rate - frame rate of the movie in frames per second
%
%   pre_frames - number of frames before trial start included in trial window
%
%   post_frames - number of frames after trial start included in trial window

%   Neurons - n x 1 array of structs, where n is the number of neurons
%   identified in the video. Each element of Neurons has the following fields:

%       Trials - t x 1 array of structs, where t is the number of
%       analyzable trials delivered during the movie (where analyzable
%       means the pre-stimulus period does not extend before the start of
%       the movie and the post-stimulus period does not extend beyond the
%       end of the movie). Each element of Trials has the following fields:

%           dFF - a 1 x p vector containing the neural activity of the
%           corresponding neuron during the corresponding peri-stimulus period,
%           where p is the number of frames in the peri-stimulus period.
%
%           Condition - char vector giving the name of the stimulus condition.

%           In addition, each element of Trials includes one field
%           corresponding to each trial parameter, reported in the Arduino
%           serial output file. E.g., if trial t reported as having a
%           STPKIDX of 1 and a SPKRIDX of 0 in the Arduino serial output file,
%           then Trials will include:

%               Trials(t).STPRIDX = 1;
%               Trials(t).SPKRIDX = 0;


%% Get some useful peri-stimulus window parameters from input struct T:
N.frame_rate = T.frame_rate;
N.pre_frames = T.pre_frames;
N.post_frames = T.post_frames;


%% Get the number of neurons in the dataset:

% Create a t x 1 vector stating the number of neurons in each trial; SHOULD BE THE SAME FOR ALL TRIALS:
n_per_trial = arrayfun(@(x) size(x.dFF, 1), T.Trials); 

% Confirm that all trials have the same number of neurons; if not, throw an error:
n_per_trial_shift = circshift(n_per_trial, 1);
if isequal(n_per_trial, n_per_trial_shift)
    n_cells = n_per_trial(1);
else
    error('Different trials have data for different numbers of neurons; please check integrity of input struct T.');
end


%% Create and populate Neurons struct:

% Extract the actual trial data from T:
Trials = T.Trials;

% Initialize struct array Neurons:
Neuron = struct();
Neurons = repmat(Neuron, n_cells, 1);

% Populate each element of Neurons with its activity data from all trials:
for n = 1:length(Neurons)
    for t = 1:length(Trials)
        % Get all param names for the current trial:
        trial_field_names = fieldnames(Trials(t));
        not_dFF = cellfun(@(c) ~strcmp(c, 'dFF'), trial_field_names);
        trial_params = trial_field_names(not_dFF);
        
        % Copy trial params to Neurons(n).Trials(t):
        for p = 1:length(trial_params)
            param_name = trial_params{p};
            Neurons(n).Trials(t).(param_name) = Trials(t).(param_name);
        end
        
        % Copy dFF data over to Neurons(n).Trials(t):
        Neurons(n).Trials(t).dFF = Trials(t).dFF(n, :);
    end
end

N.Neurons = Neurons;
