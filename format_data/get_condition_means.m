function Conditions_out = get_condition_means(Conditions_in)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-22 


%% I. OVERVIEW: 
% This function takes trialized data that has been split by condition and
% computes the mean and SEM for each condition.


%% II. REQUIREMENTS:
% 1) MATLAB v >= ???


%% III. INPUTS: 
% 1) Conditions_in - c x 1 array of structs, where c is the number of
%    conditions in which the current neuron was observed. Each element
%    of Conditions must minimally include the following fields:
%
%       Name - a concise, human-readable condition name
%
%       Trials - u x 1 array of structs, where u is the number of trials
%       of the correspinding condition delivered during the video. Each
%       element must minimally include of the following fields:
%
%           dFF - p x 1 vector of a neuron's dFF activity over the current
%           trial, where p is the length of the peri-stimulus window in
%           frames
%
%
%% IV. OUTPUTS:
% 1) Conditions_in - c x 1 array of structs, where c is the number of
%    conditions in which the current neuron was observed. Each element
%    of Conditions includes the following fields:
%
%       Name - a concise, human-readable condition name
%
%       Mean - p x 1 vector of a neuron's mean dFF activity over the
%       peri-stimulus period for the corresponding trial or stimulus
%       condition, where p is the length of the peri-stimulus window in
%       frames
%
%       SEM - p x 1 vector of a neuron's SEM over the peri-stimulus period
%       for the corresponding trial or stimulus condition, where p is the
%       length of the peri-stimulus window in frames
%
%       Trials - u x 1 array of structs, where u is the number of trials
%       of the correspinding condition delivered during the video. Each
%       element must minimally include of the following fields:
%
%           dFF - p x 1 vector of a neuron's dFF activity over the current
%           trial, where p is the length of the peri-stimulus window in
%           frames


%%
for c = 1:length(Conditions_in)
    
    % Copy over all fields form Conditions_in to Conditions_out
    condition_fields = fieldnames(Conditions_in);
    for f = 1:length(condition_fields)
        Conditions_out(c).(condition_fields{f}) = Conditions_in(c).(condition_fields{f});
    end
    
    % Compute the mean dFF trace and the SEM trace for the current condition for the current neuron:
    data = vertcat(Conditions_in(c).Trials.dFF);
    mean_dFF_trace = mean(data, 1);
    sem_trace = std(data, 1) / sqrt(size(data,1));
    
    % Copy mean dFF trace and SEM trace to output struct:
    Conditions_out(c).Mean = mean_dFF_trace;
    Conditions_out(c).SEM = sem_trace;
end