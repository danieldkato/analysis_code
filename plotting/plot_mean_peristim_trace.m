function plot_mean_peristim_trace(Conditions, frame_rate, pre_stim_frames, post_stim_frames, output_path, stim_duration, figtitle)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-25 


%% I. OVERVIEW: 
% This function plots the mean dFF trace with SEM bars over some peri-stimulus
% period for a given neuron.


%% II. REQUIREMENTS:
% 1) MATLAB v >= ???


%% III. INPUTS: 
% 1) data - c x 1 array of structs, where c is the number of
%    conditions in which the current neuron was observed. Each element
%    of Conditions has the following fields:
%           
%       Mean - 1 x p vector specifying the mean peri-stimulus dFF trace
%       for the corresponding neuron and condition, where p is the
%       duration of the peri-stimulus window in frames.
%
%       SEM - 1 x p vector specifying the SEM of the dFF trace for the
%       corresponding neuron and condition, where p is the duration of
%       the peri-stimulus window in frames.
%
%       Name - char array specifying human-readable condition name 
%
% 2) frame_rate - rate at which data was acquired, in frames per second
%
% 3) pre_stim_frames - number of frames within peri-stimulus period before
%    stimulus onset
%
% 4) post_stim_frames - number of frames within peri-stimulus period
%    including and after stimulus onset
%
% 5) conditions - c x 1 array of structs, where c is the number of
%    stimulus or trial conditions being analyzed. Each element should
%    minimimally include the following fields: 
%       Name - char array specifying human-readable condition name 
%
%       Abbreviation - char array specifying concise, human-readable
%       condition name (useful for plotting)
%
%       Color - 1 x 3 RBG vector specifying the color code for the
%       corresponding condition.
%
%     For example, an element of Conditions might be as follows:
%
%       Conditions{1}.name = 'stepper only';
%       Conditions{1}.abbreviation = 'W';
%       Conditions{1}.color = [1.00 0.00 0.00];
%       Conditions(1).params.STPRIDX = 1;
%       Conditions{1}.params.SPKRDX = 0;
%
% 6) stim_duration (optional)- stimulus duration in seconds
%
% 7) title (optional) - char vector or cell array of char vectors
%    specifying figure title


%% IV. OUTPUTS:



%% Create figure:
mean_figure = figure();
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
        
        % Retrieve mean and SEM:
        Mean = curr_condition.Mean;
        SEM = curr_condition.SEM;
        
        % Plot mean:
        hold on;
        handles(c) = plot(Mean, 'Color', color, 'LineWidth', 1.25);

        % Plot SEM:
        sem_over = [Mean fliplr(Mean + SEM)];
        sem_under = [Mean fliplr(Mean - SEM)];
        sem_over_patch = patch(sem_X, sem_over, color, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2);
        sem_under_patch = patch(sem_X, sem_under, color, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2);
        
        % Record how many trials there were of this condition:
        trials_per_condition(c) = length(curr_condition.Trials);    
end


%% Format whole figure:

% If we know the stimulus duration and it's the same for every stimulus delivery, draw a rectangle covering the stimulus period:
if ~isempty(stim_duration)
    yl = ylim;
    recY = [yl fliplr(yl)];
    %h = fill([stimPeriod(1) stimPeriod(1) stimPeriod(2) stimPeriod(2)], recY, [0.01 0.01 0.01], 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1);
    p = patch([stim_period(1) stim_period(1) stim_period(2) stim_period(2)], recY, [0.75 0.75 0.75], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
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

% Save figure:
saveas(gcf, output_path);

% Close figure;
close(mean_figure);