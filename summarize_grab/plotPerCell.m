function [meanPaths, rawPaths] = plotPerCell(params_file)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-16


%% OVERVIEW:
% This function creates peri-stimulus activity plots for every unit recorded
% from a given acquisition session. For every unit, this function creates 2
% plots:

% 1) A plot with mean activity traces plus SEM bars, and
% 2) A plot with activity traces for every individual trial

% All plots are saved to storage in an output directory specified in the
% function call.

% If the data include responses to multiple stimulus or trial conditions,
% then individual activity traces will be color-coded according to the
% mapping defined either in the data file itself or in a condition settings
% file specified in the function call (see below for detail).

% Note that this function itself is ONLY responsible for taking the mean
% and SEM of the data passed as input and plotting it. It dos NOT perform
% any other analysis, e.g. compute dF/F from raw fluorescence values.
% Thus, if the data passed as input are raw fluorescence values, then
% the values plotted will be raw fluorescence values. 


%% REQUIREMENTS:
% This function should only be used on datasets where the stimulus period
% for all trials is the same. This function draws a shaded rectangle over
% the area corresponding to the stimulus period on every plot, so the
% figures are really only sensible if the stimulus period for all trials is
% the same.


%% INPUTS:
% 1) activity - path to an HDF5 file containing the data to be plotted.
% Data must be parsed as follows: the HDF5 must include one dataset for
% each trial condition or stimulus condition to be analyzed. Each dataset
% should be an N x T x P activity matrix, where N is the number of ROIs to
% be analyzed, T is the number of samples in the peri-stimulus period to be
% plotted, and P is the number of presentations of the corresponding trial
% or stimulus condition. N and and T must be the same for all datasets, but
% P may be different for each.

% In order to properly lay out figure objects like shaded rectangles for
% stimulus epochs, the HDF5 root must have the following attributes:
    
%   a) num_samples_pre_stim: the number of frames before stimulus onset
%   included in the trace for each trial

%   b) num_samples_post_stim: the number of frames after stimulus onset
%   included in the trace for each trial

% In addition to the required attributes described above, the HDF5 may also
% have the following optional attributes for specifying the appearance of
% figures:


% 2) outputDirectory - directory where all created figures should be saved.

% 3) grabMetadata - path to a MATLAB-evaluable .txt file containing
% information about the acquisition session. This must include the
% data acquisition rate in samples per second. 

% 4) conditionSettings - path to a MATLAB-evaluable .txt file defining a 
% C x 1 cell array of structs, where C is the number of distinct trial
% conditions presented during throughout the course of the trials to be
% analyzed. Each struct must have at least the following three fields:

%   a) Name - the name of the trial condition. This must exactly match the
%   trial condition descriptions in the second column of trials, described
%   above.

%   b) Color - color code for the given condition, in any valid MATLAB
%   format for encoding color. Will be used in plotting. 

%   c) abbreviation - abbreviation for the given trial condition. Will be
%   used in creating legends for each figure.


%% OUTPUTS:
% This function saves to disk 2 plots for each ROI from the analyzed
% dataset:

% 1) A plot with mean peri-stimulus activity traces plus SEM bars for each condition, and
% 2) A plot with peri-stimulus traces for every individual trial, color-coded by condition

% In addition, this function formally returns:
% 1) meanPaths - an N x 1 cell array of strings containing full paths to
% the saved mean activity figures, and 
% 2) rawPaths - an N x 1 cell array of strings containing full paths to the
% saved individual trial activity figures. 


% TODO:
% 1) Stimulus duration is currently hard-coded into the script. This should
% also be read in dynamically from somewhere else, like trials (although
% trials currently only includes trial start frame and condition, not
% duration, so this would entail changes to how the trial matrix is
% created). In getting the trial duration from trials, one could also
% verify that all trial durations are the same (which is the only situation
% in which the shaded stimulus period rectangle makes sense).

% 2) This function requires that the condition names in trials exactly
% match the condition names in conditions. Should probably throw up an
% error if this doesn't hold. 

% 3) Maybe think about making the conditions input optional, as it's
% probably really not necessary; it could just read out the trial
% conditions from trials, then automatically assign colors and
% abbreviations to each if none are specified in conditions.

% 4) Perhaps make the outputDirectory argument optional, and have it
% default to the current working directory.


%% Get params from params file:
params = loadjson(params_file);
output_directory = params.output_directory;


%% Load condition information like color codes, etc:
C_struct = loadjson(params.conditions_path);
Conditions = C_struct.conditions;


%% Trialize data for each neuron:
N = trialize_neurons(params.rawF_path, ... 
    params.galvo_path, ...
    params.timer_path, ...
    params.ardu_path, ...
    params.conditions_path, ...
    params.grab_metadata, ...
    params.pre_sec, ...
    params.post_sec, ...
    params.show_inflection_points);


%% Get some parameters from trialized neural data that will be useful for plotting :
num_ROIs = length(N.Neurons);
frame_rate = N. frame_rate;
pre_stim_frames = N.pre_frames;
disp(['pre_stim_frames = ' num2str(pre_stim_frames)]);
post_stim_frames = N.post_frames;
disp(['post_stim_frames = ' num2str(post_stim_frames)]);
peri_stim_frames = pre_stim_frames + post_stim_frames;

% Confirm that the stimulus duration is the same for every trial, and if
% not, throw a warning and skip drawing stimulus window:
all_trial_durations = [];
for n = 1:length(N.Neurons)
    Neuron = N.Neurons(n);
    for c = 1:length(Neuron.Conditions)
        Condition = Neuron.Conditions(c);
        all_trial_durations = [all_trial_durations Condition.Trials.STIMDUR];
    end
end
equal_stimdurs = isempty(find(all_trial_durations ~= all_trial_durations(1), 1)); % test whether each reported STIMDUR matches the first
if equal_stimdurs
    stim_dur = all_trial_durations(1)/1000; % convert from milliseconds to seconds
else    
    warning('Not all stimuli have the same duration; stimulus period windows will be ommitted from plots.');
end


%% Compute some variables that will be useful for plotting:

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

% Define vectors that will be used for plotting shaded stim period rectangle (recall, stim begins at index preStim + 1):
if equal_stimdurs
    stim_period = [pre_stim_frames+1 pre_stim_frames+1+stim_dur*frame_rate];
end


%% Prepare figure windows and file paths for plotting: 

% Create cell arrays that will contain full paths to created figures; these will be returned to calling function:
meanPaths = cell(num_ROIs, 1);
rawPaths = cell(num_ROIs, 1);

% Create figure windows (these will be erased then reused between ROIs):
mean_fig = figure(); % one for mean traces
raw_fig = figure(); % one for raw traces
figures = {mean_fig, raw_fig};
titles = {'mean', 'individual'};
outputs = {meanPaths, rawPaths};

% Make sure the output directory exists, create it if it doesn't then cd into it:
status = exist(output_directory, 'dir');
if status == 0
    mkdir(output_directory);
end
old = cd(output_directory);


%% For each ROI, create 2 figures: 
% 1) a figure plotting mean traces (with SEM bars) for each condition
% 2) a figure plotting traces for all individual trials (color coded by condition)
n_figures = 1;
for n = 1:num_ROIs
    
    disp(['Plotting ROI ' num2str(n) ' out of ' num2str(num_ROIs)]);
    
    Neuron = N.Neurons(n);
    C = Neuron.Conditions;
    
    % Initialize empty vector that will state how many trials of each
    % condition were delivered:
    trials_per_condition = NaN(1,length(C));
    
    % Initialize empty vectors for storing numeric handles to plots; these
    % will be necessary for creating legends later on
    mean_plot_handles = NaN(1,length(C)); 
    raw_plot_handles = NaN(1,length(C)); 

    % For current ROI, plot data from each condition:
    for d = 1:length(C)
        
        Mean = C(d).Mean;
        SEM = C(d).SEM;

        % Get the color code for the current condition:
        cond_idx = find(cell2mat(cellfun(@(x) strcmp(C(d).Name, x.name), Conditions, 'UniformOutput', false)));
        color = Conditions{cond_idx}.color;
        
        % Plot mean of current condition for current ROI:
        figure(mean_fig);
        hold on;
        mean_plot_handles(d) = plot(Mean, 'Color', color, 'LineWidth', 1.25);

        % Plot SEM bars of current condition for current ROI:
        sem_over = [Mean fliplr(Mean + SEM)];
        sem_under = [Mean fliplr(Mean - SEM)];
        sem_over_patch = patch(sem_X, sem_over, color, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2);
        sem_under_patch = patch(sem_X, sem_under, color, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2);

        % Plot traces for every individual trial of current condition for current ROI:
        figure(raw_fig);
        hold on;
        Trials = C(d).Trials;
        for t = 1:length(Trials)
            h = plot(Trials(t).dFF, 'Color', color);
        end
        raw_plot_handles(d) = h(1);
        
        % Record how many trials there were of this condition:
        trials_per_condition(d) = length(Trials);
        
    end

    % Create legends:
    leg_text = arrayfun(@(y, z) strcat([y.Abbreviation, ', n=', num2str(z)]), C, trials_per_condition, 'UniformOutput', false);
    
    figure(mean_fig);
    legend(mean_plot_handles, leg_text);
    legend('boxoff');

    figure(raw_fig);
    legend(raw_plot_handles, leg_text);
    legend('boxoff');

    % Format and save both figures for the current ROI:
    for f = 1:length(figures)

        figure(figures{f});
        hold on;
        
        % If we know the stimulus duration and it's the same for every stimulus delivery, draw a rectangle covering the stimulus period:
        if equal_stimdurs
            yl = ylim;
            recY = [yl fliplr(yl)];
            %h = fill([stimPeriod(1) stimPeriod(1) stimPeriod(2) stimPeriod(2)], recY, [0.01 0.01 0.01], 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1);
            p = patch([stim_period(1) stim_period(1) stim_period(2) stim_period(2)], recY, [0.75 0.75 0.75], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            uistack(h, 'bottom');
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
        title(strcat(['ROI #', num2str(n), ' ', titles{f}, ' dF/F traces by condition']));

        xl = xlim;
        xlim([1 peri_stim_frames]);

        % Save figure:
        num_str = pad(num2str(n), 3, 'left', '0');
        title_str = strcat(['ROI_', num_str, '_', titles{f}, '_traces.fig']);
        fig_full_path = fullfile(output_directory, title_str);
        saveas(gcf, title_str);

        % Add the name of the figure to the appropriate outputs structure
        if mod(f,2) == 0
            rawPaths{n} = fig_full_path;
        else
            meanPaths{n} = fig_full_path;
        end

        % Clear figure for the next ROI:
        clf(gcf)
        
        % Add the full path to the figure to Metadata struct:
        Metadata.outputs(n_figures).path = fig_full_path;
        n_figures = n_figures + 1; 
    end
end


%% Write metadata:
Metadata.inputs(1).path = params.galvo_path;
Metadata.inputs(2).path = params.timer_path;
Metadata.inputs(3).path = params.ardu_path;
Metadata.inputs(4).path = params.conditions_path;
Metadata.inputs(5).path = params.grab_metadata;
Metadata.params.pre_onset_period = params.pre_sec;
Metadata.params.post_onset_period = params.post_sec;

write_metadata(Metadata, 'plot_every_cell_metadata.json');