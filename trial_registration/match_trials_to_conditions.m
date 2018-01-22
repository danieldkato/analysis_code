function Conditions = match_trials_to_conditions(Trials, Conditions)

% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-21


%% I. OVERVIEW: 
% For a given set of trials T and a given set of trial or stimulus
% conditions C, this function returns a struct indicating which trials fall
% under which condition.


%% II. REQUIREMENTS:
% 1) MATLAB >= v.???


%% III. INPUTS: 
% 1) Trials - t x 1 array of structs, where t is the number of trials to be
%    assigned condition labels. Each element of Trials should minimally
%    include one field corresponding to each trial parameter, reported in the
%    Arduino serial output file. E.g., if trial t reported as having a STPKIDX
%    of 1 and a SPKRIDX of 0 in the Arduino serial output file, then Trials
%    will include:
%
%       Trials(t).STPRIDX = 1;
%       Trials(t).SPKRIDX = 0;

% 2) Conditions - c x 1 cell array of structs, where c is the number of
%    stimulus or trial conditions being analyzed. Each element should
%    minimimally include a "name" field that states the name of the
%    condition, as well as a "params" field. The "params" field should
%    itself have one sub-field for each parameter that defines the
%    condition; each of these sub-fields should be named after the
%    corresponding trial parameter reported in the arduino serial output
%    file. For example, an element of Conditions might be as follows:

%       Conditions{1}.name = 'stepper only';
%       Conditions{1}.abbreviation = 'W';
%       Conditions(1).params.STPRIDX = 1;
%       Conditions{1}.params.SPKRDX = 0;


%% IV. OUTPUTS:
% 2) Conditions - same as input struct Conditions, but each element
% includes an additional field called matching_trials, a u x 1 vector
% specifying the indices into Trials that specifies which trials belong to
% the corresponding condition.


%%
for c = 1:length(Conditions)
    
    % Initialize empty vector of matching trial indices and filter
    matching_trials = [];
    filter = ones(1, length(Trials));
    
    % Go through every trial parameter that defines the current condition
    % and find the trials that match ~all~ of them:
    params = fieldnames(Conditions{c}.params);  
    for p = 1:length(params) 
        param_name = params{p};
        param_value = Conditions{c}.params.(param_name);
        filter_update = [Trials.(param_name)] == param_value;
        filter = filter & filter_update;
    end
    matching_trials = find(filter);
    Conditions{c}.matching_trials = matching_trials;
    
    % Warn the user if no trials of the current condition are found:
    if isempty(matching_trials)
        warning(['No trials of condition ' Conditions{c}.name ' found.']);
    end
end