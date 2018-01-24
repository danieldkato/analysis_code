function Conditions_out = get_condition_means(Conditions_in)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-23 


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
%           dFF - n x p matrix of a neuron or neural population's dFF
%           activity over the current trial, where n is the number of
%           neurons and p is the length of the peri-stimulus window in
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
%       Mean - n x p matrix of giving a neuron or neural population's mean
%       dFF activity over the peri-stimulus period for the corresponding
%       trial or stimulus condition, where n is the number of neurons and p
%       is the length of the peri-stimulus window in frames
%
%       SEM - 1 x p vector of a neuron's SEM over the peri-stimulus period
%       for the corresponding trial or stimulus condition, where p is the
%       length of the peri-stimulus window in frames
%
%       Trials - u x 1 array of structs, where u is the number of trials
%       of the correspinding condition delivered during the video. Each
%       element must minimally include of the following fields:
%
%           dFF - n x p matrix of a neuron or neural population's dFF
%           activity over the current trial, where n is the number of
%           neurons and p is the length of the peri-stimulus window in
%           frames


%%

% Get the names of all fields for each element of the input conditions
% struct:
condition_fields = fieldnames(Conditions_in);

for c = 1:length(Conditions_in)
    
    curr_cond_trials = Conditions_in(c).Trials;
    n_trials = length(curr_cond_trials);
    
    % Frist do some validation:
    
    % Verify that all observations for the current condition include the
    % same number of neurons; if so, determine whether this is
    % single-neuron or population data; if not, raise an error.
    dFF_dim1 = arrayfun(@(x) size(x.dFF, 1), curr_cond_trials);
    check_dim1 = circshift(dFF_dim1, 1);
    if isequal(dFF_dim1, check_dim1)
        n_cells = size(curr_cond_trials(1).dFF, 1);
    else
        error(['Different numbers of neurons observed in different trials of condition ' Conditions_in(c).Name '; without explicit identification of neurons between trials, cannot average across trials.']);
    end
    
    % Verify that all observations for the current condition have the same
    % number of frames; if not, raise an error:
    dFF_dim2 = arrayfun(@(x) size(x.dFF, 2), curr_cond_trials);    
    check_dim2 = circshift(dFF_dim2, 1);
    if isequal(dFF_dim2, check_dim2)
        n_cells = size(curr_cond_trials(1).dFF, 2);
    else
        error(['Different numbers of frames observed in different trials of condition ' Conditions_in(c).Name '; cannot compute mean peri-stimulus dFF trace.']);
    end
    
    % Copy over all fields form Conditions_in to Conditions_out
    for f = 1:length(condition_fields)
        Conditions_out(c).(condition_fields{f}) = Conditions_in(c).(condition_fields{f});
    end

    % Compute the mean dFF trace and the SEM trace for the current condition for the current neuron:
    if n_cells == 1
        data = vertcat(curr_cond_trials.dFF);
        mean_dFF_trace = mean(data, 1);
        sem_trace = std(data, 1) / sqrt(size(data,1));
    elseif n_cells > 1
        data = NaN(n_cells, n_frames, n_trials);
        for t = 1:n_trials
            data(:,:,t) = conditions_trials(t).dFF(:,:);
        end 
        mean_dFF_trace = mean(data, 3);
        sem_trace = std(data, 3) / sqrt(size(data,3));
    end
    
    % Copy mean dFF trace and SEM trace to output struct:
    Conditions_out(c).Mean = mean_dFF_trace;
    Conditions_out(c).SEM = sem_trace;
end