function [f, S] = compare_2_conditions(params_file, cond1_name, cond2_name)
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

Conditions = split_trials_by_condition(T.Trials, Conditions); % Split trials by condition
Conditions = get_condition_means(Conditions); % Get the mean peristimulus dF/F trace of every neuron for every condition;

% Get the mean of the trace post-stim onset:
for c = 1:length(Conditions)
	Conditions(c).distribution = mean(Conditions(c).Mean(:, pre_stim_frames+1:end), 2);
end

Condition1 = get_condition(cond1_name, Conditions);
Condition2 = get_condition(cond2_name, Conditions);


%% Run a t-test on the 2 conditions of interest and plot a histogram:
fig_title = {['\color[rgb]{' num2str(Condition1.color) '}' Condition1.name ' \color[rgb]{0 0 0} vs \color[rgb]{' num2str(Condition2.color) '}' Condition2.name]; ['\color[rgb]{0 0 0} Mouse ' mouse ', session ' date]};
[f, p, stats] = ttest_and_histogram(Condition1, Condition2, true, fig_title);
S.stats = stats;
S.p = p;


%% Save figure and stats:
fig_path = [output_directory filesep Condition1.name '_vs_' Condition2.name '_distributions.fig'];
stats_path = [output_directory filesep Condition1.name '_vs_' Condition2.name '.mat'];

savefig(f, fig_path);
save(stats_path, 'S');


%% Save metadata: 
Metadata.inputs(1).path = params.rawF_path;
Metadata.inputs(2).path = params.galvo_path;
Metadata.inputs(3).path = params.timer_path;
Metadata.inputs(4).path = params.ardu_path;
Metadata.inputs(5).path = params.conditions_path;
Metadata.inputs(6).path = params.grab_metadata;
Metadata.params.pre_onset_period = params.pre_sec;
Metadata.params.post_onset_period = params.post_sec;
Metadata.outputs(1).path = fig_path;
Metadata.outputs(2).path = stats_path;

write_metadata(Metadata, [output_directory filesep 'regress_neurons_v_trial_params_within_session.json']);