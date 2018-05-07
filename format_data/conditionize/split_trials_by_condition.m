function Conditions = split_trials_by_condition(Trials, Conditions)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-22


%% I. OVERVIEW: 
% This function takes a struct describing a series of trials and splits it
% up by condition, yielding a struct Conditions_out that is organized by
% condition at the top level and then by trial at a subordinate level. See
% OUTPUT below for more detail.


%% II. REQUIREMENTS:
% 1) match_trials_to_conditions.m


%% III. INPUTS: 
% 1) Trials - t x 1 array of structs, where t is the number of trials to be
%    analyzed . Each element of Trials should minimally include one field
%    corresponding to each trial parameter reported in the Arduino serial
%    output file. E.g., if trial t reported as having a STPKIDX of 1 and a
%    SPKRIDX of 0 in the Arduino serial output file, then Trials should
%    include:
%
%       Trials(t).STPRIDX = 1; 
%       Trials(t).SPKRIDX = 0;

%   In addition, each element of trials may include a field containing
%   neural activity data of some kind, either for a single neuron or for a
%   population.


% 2) Conditions_in - c x 1 array of structs, where c is the number of
%    stimulus or trial conditions being analyzed. Each element should
%    minimimally include the following fields:
%
%       name - char vectory specifying a human-readable condition name
%
%       params - a struct with one sub-field for each parameter that
%       defines the condition; each of these sub-fields should be named
%       after the corresponding trial parameter reported in the arduino
%       serial output file.
%
%       For example, an element of Conditions might be as follows:
%
%       Conditions(1).name = 'stepper only';
%       Conditions(1).abbreviation = 'W';
%       Conditions(1).params.STPRIDX = 1;
%       Conditions(1).params.SPKRDX = 0;


%% IV. OUTPUTS:
% 2) Conditions_out - same as input struct Conditions_in, but each element
% additionally includes:

%   Trials - u x 1 array of structs, where u is the number of trials that
%   fall under the corresponding condition. Each element of
%   Conditions_out.Trials has the same structure as an element of the input
%   struct Trials.


%%
Conditions = match_trials_to_conditions(Trials, Conditions);

for c = 1:length(Conditions)
    
    % Initialize empty vector of matching trial indices and filter
    matching_trials = [];
    filter = ones(1, length(Trials));
    
    % Go through every trial parameter that defines the current condition
    % and find the trials that match ~all~ of them:
    params = fieldnames(Conditions(c).params);  
    for p = 1:length(params) 
        param_name = params{p};
        param_value = Conditions(c).params.(param_name);
        filter_update = [Trials.(param_name)] == param_value;
        filter = filter & filter_update;
    end
    matching_trials = find(filter);
    
    % Warn the user if no trials of the current condition are found:
    if isempty(matching_trials)
        warning(['No trials of condition ' Conditions(c).name ' found.']);
    end
    
    Conditions(c).Trials = Trials(matching_trials);
end