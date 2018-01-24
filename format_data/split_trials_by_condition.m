function Conditions_out = split_trials_by_condition(Trials, Conditions_in)
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
Conditions_matched = match_trials_to_conditions(Trials, Conditions_in);

C = struct();
Conditions_out = repmat(C, 1, length(Conditions_in));

for c = 1:length(Conditions_matched)
    
    % Find all Condition metadata fields (name, abbreviation, color code):
    condition_fields = fieldnames(Conditions_in);
    not_Trials = cellfun(@(c) ~strcmp(c, 'Trials'), condition_fields);
    condition_fields = condition_fields(not_Trials);
    
    % Copy all Condition metadata fields to output struct:
    for f = 1:length(condition_fields)
        Conditions_out(c).(condition_fields{f}) = Conditions_in(c).(condition_fields{f});
    end
    
    % Copy over the trials data from the matching condition:
    matching_trials = Conditions_matched(c).matching_trials;
    Conditions_out(c).Trials = Trials(matching_trials);
end