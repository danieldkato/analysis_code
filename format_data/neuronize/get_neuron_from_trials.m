function Trials_out = get_neuron_from_trials(Trials_in, n, field)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-05-21


%% I. OVERVIEW: 
% This function returns a struct T_out contianing trialized data for a
% single neuron. See OUTPUT below for more detail.


%% II. REQUIREMENTS:
% 1) MATLAB v >= ???


%% III. INPUTS: 
% 1) Trials_in - t x 1 array of structs, where t is the number of analyzable
%    trials delivered during the movie (where analyzable means the
%    pre-stimulus period does not extend before the start of the movie and
%    the post-stimulus period does not extend beyond the end of the movie).
%    Each element of Trials must minimally include the following fields:

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
% 1) Trials_out - struct containing trialized data for a session. T includes the
% following fields:
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
if nargin < 3
    field = 'dFF';
end

for t = 1:length(Trials)
    % Replace data matrix for all cells with data matrix just for selected neuron:
    data_matrix = Trials_in(t).(field);
    Trials_out(t).(field) = data_matrix(n, :);
end