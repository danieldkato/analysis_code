function [f, S] = compare_2_conditions(params_file, cond1_name, cond2_name, paired)
%% Load parameters and metadata:

% Load params from params file:
params = loadjson(params_file);
output_directory = params.output_directory;

% Load some useful imaging session meadata:
grab_metadata = loadjson(params.grab_metadata);
mouse = grab_metadata.mouse;
date = grab_metadata.date;

% Load condition information like name, hardware parameters, color codes, etc:
C_struct = loadjson(params.conditions_path);
Conditions = cell2mat(C_struct.conditions);

% Compute some quantities that will be needed later on:
frame_rate = grab_metadata.frame_rate;
pre_stim_frames = ceil(params.pre_sec * grab_metadata.frame_rate);


%% Do initial processing:

% Split data by trial:
T = trialize_data(params.rawF_path, ... 
    params.galvo_path, ...
    params.timer_path, ...
    params.ardu_path, ...
    params.conditions_path, ...
    params.grab_metadata, ...
    params.pre_sec, ...
    params.post_sec, ...
    params.show_inflection_points);

stim_duration = T.Trials(1).STIMDUR;
end_frame = pre_stim_frames + frame_rate * stim_duration;

Conditions = split_trials_by_condition(T.Trials, Conditions); % Split trials by condition
Conditions = get_condition_means(Conditions); % Get the mean peristimulus dF/F trace of every neuron for every condition;
Conditions = get_condition_amplitudes(Conditions, pre_stim_frames, end_frame);

% Get the mean of the trace post-stim onset:
for c = 1:length(Conditions)
	Conditions(c).distribution = mean(Conditions(c).Mean(:, pre_stim_frames+1:end), 2);
end

Condition1 = get_condition(cond1_name, Conditions);
Condition2 = get_condition(cond2_name, Conditions);


%% Run a t-test on the 2 conditions of interest and plot a histogram:
fig_title = {['Peak dF/F response for \color[rgb]{' num2str(Condition1.color) '}' Condition1.name ' \color[rgb]{0 0 0} vs \color[rgb]{' num2str(Condition2.color) '}' Condition2.name]; ['\color[rgb]{0 0 0} Mouse ' mouse ', session ' date]};
[f, p, stats] = ttest_and_histogram(Condition1, Condition2, paired, fig_title);
S.ttest.stats = stats;
S.ttest.p = p;


%% If the condition data are paired, do a regression and create a scatter plot:


if paired
    % Do a linear regression; we want to test whether the slope of the best fit
    % line for Condition1 vs Condition2 is significantly different from 1, but
    % MATLAB's fitlm function tests whether the slope of a best fit line is
    % significantly different from 0. So instead of regressing Condition1 on
    % Condition2, we'll regress (Condition1 - Condition2) on to condition
    % Condition2; if Condition1 = Condition2, then the slope of
    % (Condition1-Condition2) vs Condition2 should be 0:
    Tbl = table();
    Tbl.X = Condition2.amplitudes;
    Tbl.Y = Condition1.amplitudes - Condition2.amplitudes;
    lm = fitlm(Tbl, 'Y ~ X');    
    disp(lm);
    S.linear_model = lm;
    
    % Create a scatterplot of each cell's response to both conditions:
    %scatter_fig_title = {['Peak dF/F response to \color[rgb]{' num2str(Condition1.color) '}' Condition1.name '\color[rgb]{0 0 0} vs \color[rgb]{' num2str(Condition2.color) '}' Condition2.name]; ['\color[rgb]{0 0 0} \fontsize{10}Mouse ' mouse]; ['\fontsize{10}Session ' date]};
    scatter_handle = scatter_conditions(cond1_name, cond2_name, Conditions, 'amplitudes', fig_title); % Create scatter plot of W+T1 vs W
    lims = [xlim ylim];
    lowest = min(lims);
    highest = max(lims);
    xlim([lowest highest]);
    ylim([lowest highest]);
    hline = refline(1, 0);
    hline.Color = [0.40 0.40 0.40];
    lsline; % add least-squares line
end


%% Save figures and stats:

% Save stats:
stats_path = [output_directory filesep Condition1.name '_vs_' Condition2.name '.mat'];
save(stats_path, 'S');

% Save histogram figure:
hist_fig_path = [output_directory filesep Condition1.name '_vs_' Condition2.name '_distributions.fig'];
savefig(f, hist_fig_path);

% Save scatterplot figure:
if paired
    scatter_fig_path = [output_directory filesep Condition1.abbreviation '_v_' Condition2.abbreviation '_scatter.fig']; 
    savefig(scatter_handle, scatter_fig_path);
end


%% Save metadata: 
Metadata.inputs(1).path = params.rawF_path;
Metadata.inputs(2).path = params.galvo_path;
Metadata.inputs(3).path = params.timer_path;
Metadata.inputs(4).path = params.ardu_path;
Metadata.inputs(5).path = params.conditions_path;
Metadata.inputs(6).path = params.grab_metadata;
Metadata.params.pre_onset_period = params.pre_sec;
Metadata.params.post_onset_period = params.post_sec;
Metadata.outputs(1).path = hist_fig_path;
Metadata.outputs(2).path = stats_path;
if paired
    Metadata.outputs(3).path = scatter_fig_path;
end 

write_metadata(Metadata, [output_directory filesep Condition1.name '_vs_' Condition2.name '.json']);