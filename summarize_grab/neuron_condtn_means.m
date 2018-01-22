function M = neuron_condtn_means(N)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-21 


%% I. OVERVIEW: 
% This function returns a struct N containing the mean dFF trace and SEM
% over some peri-stimulus period for each neuron for each condition from a
% given imaging session. See OUTPUT below for more detail.


%% II. REQUIREMENTS:
% 1) MATLAB v >= ???


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
%           Mean - 1 x p vector specifying the mean peri-stimulus dFF trace
%           for the corresponding neuron and condition, where p is the
%           duration of the peri-stimulus window in frames.
%
%           SEM - 1 x p vector specifying the SEM of the dFF trace for the
%           corresponding neuron and condition, where p is the duration of
%           the peri-stimulus window in frames.
%
%       In addition, each element of M.Neurons.Conditions will include any
%       other fields that were included in N.Neurons.Conditions. For
%       plotting purposes, these should include:
%
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

    curr_neuron = N.Neurons(n);
    curr_neuron_conditions = curr_neuron.Conditions;
    
    % For each condition in which the current neuron was observed:
    for c = 1:length(curr_neuron_conditions)
    
        curr_condition = curr_neuron_conditions(c);
            
        % Compute the mean dFF trace and the SEM trace for the current condition for the current neuron:
        data = vertcat(curr_condition.Trials.dFF);
        mean_dFF_trace = mean(data, 1);
        sem_trace = std(data, 1) / sqrt(size(data,1));
    
        % Copy mean dFF trace and SEM trace to output struct:
        M.Neurons(n).Conditions(c).Mean = mean_dFF_trace;
        M.Neurons(n).Conditions(c).SEM = sem_trace;
        
        % Copy over any other condition parameters to output struct:
        curr_condition_params = fieldnames(curr_condition);
        for pp = 1:length(curr_condition_params)
            M.Neurons(n).Conditions(c).(curr_condition_params{pp}) = N.Neurons(n).Conditions(c).(curr_condition_params{pp});
        end
        
    end
    
end