function M = neuron_condtn_means(N, Conditions)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-21 


%% I. OVERVIEW: 
% This function returns a struct N containing trialized data for each
% neuron from a given imaging session. See OUTPUT below for more detail.


%% II. REQUIREMENTS:
% 1) match_trials_to_conditions.m


%% III. INPUTS: 
% 1) N - struct containing trialized data for each neuron. N includes the
% following fields:
%
%   frame_rate - frame rate of the movie in frames per second
%
%   pre_frames - number of frames before trial start included in trial window
%
%   post_frames - number of frames after trial start included in trial window
%
%   Neurons - n x 1 array of structs, where n is the number of neurons
%   identified in the movie. Each element of Neurons has the following
%   fields:
%
%           Trials - t x 1 array of structs, where t is the number of
%           trials delivered during the movie. Each element of Trials
%           consists of the following fields:
%
%               dFF - q x 1 vector of the current neuron's dFF activity
%               over the current trial, where q is the length of the
%               peri-stimulus window in frames
%
%               start_frame - frame during which the stimulus was delivered
%
%           as well one field corresponding to each trial parameter
%           reported in the Arduino serial output file.

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
% 1) M - struct containing trialized data for each neuron. N includes the
% following fields:
%
%   frame_rate - frame rate of the movie in frames per second
%
%   pre_frames - number of frames before trial start included in trial window
%
%   post_frames - number of frames after trial start included in trial window
%
%   Neurons - n x 1 array of structs, where n is the number of neurons
%   identified in the movie. Each element of Neurons has the following
%   fields:
%
%       Conditions - c x 1 array of structs, where c is the number of
%       conditions in which the current neuron was observed. Each element
%       of Conditions has the following fields:
%
%           Trials - p x 1 array of structs, where p is the number of
%           trials *of the current condition* delivered during the movie.
%           Each element consists of the following fields:
%
%               dFF - q x 1 vector of the current neuron's dFF activity
%               over the current trial, where q is the length of the
%               peri-stimulus window in frames
%
%               start_frame - frame during which the stimulus was delivered
%
%           as well one field corresponding to each trial parameter
%           reported in the Arduino serial output file.
%
%       In addition, each element of M.Neurons.Conditions will include whatever
%       fields are included for the corresponding condition in the input
%       Conditions structure. For plotting purposes, these should include:

%           Name - char array specifying human-readable condition name 
%
%           Abbreviation - char array specifying concise, human-readable
%           condition name (useful for plotting)


%% Get some metadata:
M.frame_rate = N.frame_rate;
M.pre_frames = N.pre_frames;
M.post_frames = N.post_frames;


%% Split information for each neuron by condition:

% For each neuron:
for n = 1:length(N.Neurons)
    
    % Get the indices of which trials belong to each condition:
    matching_conditions = match_trials_to_conditions(Conditions, Neurons(n).Trials);
    
    % For each condition:
    for c = 1:length(matching_conditions)
        
        curr_condition = matching_conditions{c};
        
        % Fetch any trials that belong to the current condition:
        curr_cond_trial_indices = curr_condition.matching_trials;
        curr_cond_trials = Neurons(n).Trials(curr_cond_trial_indices); 
        
        % If the current neuron was not observed in ANY trials of the
        % current condition, then make M.Neurons(n).Conditions(c).Trials an
        % empty array and throw a warining:
        if isempty(curr_cond_trials)
            M.Neurons(n).Conditions(c).Trials = [];
            warning(['Neuron ' num2str(n) ' not observed in any trials of condition ' curr_condition.name]);
        end
        
        % For each trial that falls within the current condition:
        for t = 1:length(curr_cond_trials)
            
            curr_trial = curr_cond_trials(t);
            
            % Copy all of the trial data, including parameters and dFF
            % info, to the output struct:
            curr_trial_params = fieldnames(curr_trial);            
            for p = 1:length(curr_trial_params)
                M.Neurons(n).Conditions(c).Trials(t).(curr_trial_params{p}) = curr_trial.(curr_trial_params{p});
            end            
        end
       
        % Add general condition information to output struct:
        curr_condition_params = fieldnames(curr_condition);
        for pp = 1:length(curr_condition_params)
            M.Neurons(n).Conditions(c).(curr_condition_params{pp}) = curr_condition.(curr_condition_params{pp});
        end
                
    end
    
end