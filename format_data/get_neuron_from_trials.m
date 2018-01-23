function T_out = get_neuron_from_trials(T_in, n)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-22


%% I. OVERVIEW: 
% This function returns a struct T_out contianing trialized data for a
% single neuron. See OUTPUT below for more detail.


%% II. REQUIREMENTS:
% 1) MATLAB v >= ???


%% III. INPUTS: 
% 1) T_in - struct containing trialized data for a session. T includes the
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
%   Each element of Trials must minimally include the following fields:

%       dFF - an n x p matrix containing the population neural activity
%       during the corresponding peri-stimulus period, where n is the
%       number of neurons and p is the number of frames in the
%       peri-stimulus period.

%       In addition, each element of Trials should include one field
%       corresponding to each trial parameter reported in the Arduino
%       serial output file. E.g., if trial t reported as having a
%       STPKIDX of 1 and a SPKRIDX of 0 in the Arduino serial output file,
%       then Trials will include:

%           Trials(t).STPRIDX = 1;
%           Trials(t).SPKRIDX = 0;

% 2) n - index of the neuron to extract dFF data for.


%% IV. OUTPUTS:
% 1) T_out - struct containing trialized data for a session. T includes the
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
%   Each element of Trials must minimally include the following fields:

%       dFF - a 1 x p vector containing the dFF activity of the neuron
%       specified by input argument n during the peri-stimulus period of
%       the corresponding trial, where p is the number of frames in the
%       peri-stimulus period.

%       In addition, each element of Trials should include one field
%       corresponding to each trial parameter reported in the Arduino
%       serial output file. E.g., if trial t reported as having a
%       STPKIDX of 1 and a SPKRIDX of 0 in the Arduino serial output file,
%       then Trials will include:

%           Trials(t).STPRIDX = 1;
%           Trials(t).SPKRIDX = 0;


%%
T_out.frame_rate = T_in.frame_rate;
T_out.pre_frames = T_in.pre_frames;
T_out.post_frames = T_in.post_frames;

for t = 1:length(T_in.Trials)
    
    % Copy dFF just for selected neuron to output struct:
    T_out.Trials(t).dFF = T_in.Trials(t).dFF(n, :);
    
    % Find all non-dFF trial params:
    trial_params = fieldnames(T_in.Trials(t));
    not_dFF = cellfun(@(c) ~strcmp(c, 'dFF'), trial_field_names);
    trial_params = trial_params(not_dFF);
    
    % Copy all non-dF trial params to T_out.Trials(t):
    for p = 1:length(trial_params)
        param_name = trial_params{p};
        T_out.Trials(t).(param_name) = T_in.Trials(t).(param_name);
    end
    
end