function Conditions = get_condition_amplitudes(Conditions, start, stop)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-29 


%% I. OVERVIEW: 
% This function takes trialized data that has been split by condition and
% gets the amplitude of each neuron's mean trace for each condition.


%% II. REQUIREMENTS:
% 1) MATLAB v >= ???


%% III. INPUTS: 
% 1) Conditions - c x 1 array of structs, where c is the number of
%    conditions being analyzed. Each element of Conditions must minimally
%    include the following fields:
%
%       Mean - n x p matrix giving mean peri-stimulus dF/F activity traces,
%       where n is the number of neurons and p is the length of the
%       peri-stimulus window in frames. Each row represents a single
%       neuron's mean dF/F response to the corresponding condition.
%
% 2) start - start frame of time window within which to measure each trace's amplitude
%
% 3) stop (optional) - end frame of time window within which to take measure
%    trace's amplitude. If not specified by user, this will be set to the
%    end of each trace.
%
%
%% IV. OUTPUTS:
% 1) Conditions - same as the input, but each element has the following
% additional field:
%
%       amplitudes - n x 1 vector giving the amplitude of each neuron's
%       mean response to the corresponding condition, where n is the number
%       of neurons and p is the length of the peri-stimulus window in
%       frames.
%

%%
for c = 1:length(Conditions)
    
    mean_traces = Conditions(c).Mean;
    baselines = mean(mean_traces(:, 1:start-1), 2);
    
    if nargin < 2
        abs_peaks = max(Conditions(c).Means(:, start:stop), [], 2);
    else
        abs_peaks = max(Conditions(c).Means(:, start:end), [], 2);
    end
    
    Conditions(c).amplitudes = abs_peaks - baselines;
end