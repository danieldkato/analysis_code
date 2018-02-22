function single_trials_figure = plot_individual_peristim_traces(Conditions, pre_stim_frames, post_stim_frames, stim_duration, frame_rate, figtitle)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-25 


%% I. OVERVIEW: 
% This function plots individual dFF traces over some peri-stimulus period
% for a given neuron for each stimulus or trial condition in which the
% neuron was observed.


%% II. REQUIREMENTS:
% 1) MATLAB v >= ???


%% III. INPUTS: 
% 1) data - c x 1 array of structs, where c is the number of
%    conditions in which the current neuron was observed. Each element
%    of Conditions must minimally include the following fields:
%
%       Name - char array specifying human-readable condition name 
%
%       Abbreviation - char array specifying concise, human-readable
%       condition name (useful for plotting)
%
%       Color - 1 x 3 RGB vector specifying the color code for the
%       corresponding condition.
%
%       Trials - 1 x u array of structs, where u is the
%       number of trials of the corresponding conditions. Each element of
%       Trials must minimally include the following fields:
%
%           dFF - 1 x p vector specifying the peri-stimulus dFF trace for
%           the corresponding neuron and trial, where p is the duration of
%           the peri-stimulus window in frames.
%
% 2) pre_stim_frames - number of frames within peri-stimulus period before
%    stimulus onset
%
% 3) post_stim_frames - number of frames within peri-stimulus period
%    including and after stimulus onset
%
% 4) frame_rate - rate at which data was acquired, in frames per second
%
% 5) stim_duration (optional)- stimulus duration in seconds
%
% 6) title (optional) - char vector or cell array of char vectors
%    specifying figure title


%% IV. OUTPUTS:
% 1) single_trials_figure - handle to created figure


%% Create figure:
single_trials_figure = figure();
hold on;


%% Compute some quantities that will be useful for plotting later on:
peri_stim_frames = pre_stim_frames + post_stim_frames;

% Compute appropriate x-tick coordinates (in data units) so that there is one tick at stimulus onset and ticks marks every 2 sec:
seconds_per_tick = 2;     
frames_per_tick = frame_rate*seconds_per_tick;

n_ticks_before_onset = floor((pre_stim_frames)/frames_per_tick);     
ticks_before_onset = (pre_stim_frames+1) - fliplr(frames_per_tick*(0:1:n_ticks_before_onset)); % remember, pre_stim_frames + 1 is the frame where the stim actually starts 

n_ticks_after_onset = floor((post_stim_frames)/frames_per_tick);    
ticks_after_onset = (pre_stim_frames+1) + frames_per_tick*(1:1:n_ticks_after_onset);

xticks = [ticks_before_onset ticks_after_onset];

% Write x-tick labels into a cell array:
labels = (xticks - (pre_stim_frames+1))/frame_rate;
labels = arrayfun(@(a) num2str(a), labels, 'UniformOutput', 0);

% Define some vectors that will be used for plotting shaded SEM areas:
domain = (1:1:pre_stim_frames+post_stim_frames);
sem_X = [domain fliplr(domain)]; 

if ~isempty(stim_duration)
    stim_period = [pre_stim_frames+1 pre_stim_frames+1+stim_duration*frame_rate];
end


%% Plot data for each condition:
n_conditions = length(Conditions);
trials_per_condition = nan(1, n_conditions);

for c = 1:n_conditions
    
        curr_condition = Conditions(c);
    
        % Get the color code for the current condition:
        color = curr_condition.color;
        
        % Plot individual traces:
        for t = 1:length(curr_condition.Trials)
            handles(c) = plot(curr_condition.Trials(t).dFF, 'Color', color, 'LineWidth', 1.25);
        end

        % Record how many trials there were of this condition:
        trials_per_condition(c) = length(curr_condition.Trials);    
end


%% Format whole figure:

% If we know the stimulus duration and it's the same for every stimulus delivery, draw a rectangle covering the stimulus period:
if ~isempty(stim_duration)
    yl = ylim;
    recY = [yl fliplr(yl)];
    %h = fill([stimPeriod(1) stimPeriod(1) stimPeriod(2) stimPeriod(2)], recY, [0.01 0.01 0.01], 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1);
    p = patch([stim_period(1) stim_period(1) stim_period(2) stim_period(2)], recY, [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    uistack(p, 'bottom');
    ylim(yl);
end

% Plot y = 0 for reference:
plot(zeros(1,peri_stim_frames), 'Color', [0.75 0.75 0.75], 'LineWidth', 0.1);

% Position and label x-ticks appropriately:
set(gca, 'XTick', xticks, 'XTickLabel', labels);            

% Label x- and y- axes:
ylabel('dF/F (a.u.)');
xlabel('Time rel. to stim onset (s)');

% Print title:
if ~isempty(figtitle)
    title(figtitle);
end

xl = xlim;
xlim([1 peri_stim_frames]);

% Add legend:     
leg_text = arrayfun(@(y, z) strcat([y.abbreviation, ', n=', num2str(z)]), Conditions, trials_per_condition, 'UniformOutput', false);
legend(handles, leg_text);
legend('boxoff');