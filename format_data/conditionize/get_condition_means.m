function Conditions = get_condition_means(Conditions)
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
% 1) Conditions - c x 1 array of structs, where c is the number of
%    conditions being analyzed. Each element of Conditions must minimally
%    include the following fields:
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
%           frames.
%
%
%% IV. OUTPUTS:
% 1) Conditions - same as the input, but each element has the following
% additional field:
%
%       Mean - n x p matrix giving mean peri-stimulus dF/F activity traces,
%       where n is the number of neurons and p is the length of the
%       peri-stimulus window in frames. Each row represents a single
%       neuron's mean dF/F response to the corresponding condition.


%%

% Get the names of all fields for each element of the input conditions
% struct:
condition_fields = fieldnames(Conditions);

for c = 1:length(Conditions)
    
    curr_cond_trials = Conditions(c).Trials;
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
        error(['Different numbers of neurons observed in different trials of condition ' Conditions(c).Name '; without explicit identification of neurons between trials, cannot average across trials.']);
    end
    
    % Verify that all observations for the current condition have the same
    % number of frames; if not, raise an error:
    dFF_dim2 = arrayfun(@(x) size(x.dFF, 2), curr_cond_trials);    
    check_dim2 = circshift(dFF_dim2, 1);
    if isequal(dFF_dim2, check_dim2)
        n_frames = size(curr_cond_trials(1).dFF, 2);
    else
        error(['Different numbers of frames observed in different trials of condition ' Conditions(c).Name '; cannot compute mean peri-stimulus dFF trace.']);
    end

    % Compute the mean dFF trace and the SEM trace for the current condition for the current neuron:
    if n_cells == 1
        data = vertcat(curr_cond_trials.dFF);
        mean_dFF_trace = mean(data, 1);
        sem_trace = std(data, 0, 1) / sqrt(size(data,1));
    elseif n_cells > 1
        data = NaN(n_cells, n_frames, n_trials);
        for t = 1:n_trials
            data(:,:,t) = curr_cond_trials(t).dFF(:,:);
        end 
        mean_dFF_trace = mean(data, 3);
        sem_trace = std(data, 0, 3) / sqrt(size(data,3));
    end
    
    % Copy mean dFF trace and SEM trace to output struct:
    Conditions(c).Mean = mean_dFF_trace;
    Conditions(c).SEM = sem_trace;
end