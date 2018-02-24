function [Stats, hist_fig_handle, scatter_handle] = compare_2_conditions_within_session(params_file, cond1_name, cond2_name, field_name)
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


%% Perform ttest and regression on the 2 conditions:
if strcmp(field_name, 'amplitudes')
    title_root = 'Peak';
elseif strcmp(field_name, 'Mean')
    title_root = 'Mean';
end
fig_title = {[title_root ' dF/F response for \color[rgb]{' num2str(Condition1.color) '}' Condition1.name ' \color[rgb]{0 0 0} vs \color[rgb]{' num2str(Condition2.color) '}' Condition2.name]; ['\color[rgb]{0 0 0} Mouse ' mouse ', session ' date]};
[ttest_results, hist_fig_handle, regression_results, scatter_handle] = compare_2_conditions(Condition1, Condition2, field_name, true, fig_title);
Stats.ttest = ttest_results;
Stats.linear_model = regression_results;


%% Save figures and stats:

% Make output directory if it doesn't exist:
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end
output_basename = [Condition1.abbreviation '_vs_' Condition2.abbreviation];
mkdir([output_directory filesep output_basename]);

% Save stats:
stats_path = [output_directory filesep output_basename filesep output_basename '.mat'];
save(stats_path, 'Stats');

% Save histogram figure:
hist_fig_path = [output_directory filesep output_basename filesep output_basename '_distributions.fig'];
savefig(hist_fig_handle, hist_fig_path);

% Save scatterplot figure:
scatter_fig_path = [output_directory filesep output_basename filesep output_basename '_scatter.fig']; 
savefig(scatter_handle, scatter_fig_path);



%% Save metadata: 
Metadata.inputs(1).path = params.rawF_path;
Metadata.inputs(2).path = params.galvo_path;
Metadata.inputs(3).path = params.timer_path;
Metadata.inputs(4).path = params.ardu_path;
Metadata.inputs(5).path = params.conditions_path;
Metadata.inputs(6).path = params.grab_metadata;
Metadata.params.pre_onset_period = params.pre_sec;
Metadata.params.post_onset_period = params.post_sec;
Metadata.params.compared_variable = field_name;
Metadata.outputs(1).path = hist_fig_path;
Metadata.outputs(2).path = stats_path;
if paired
    Metadata.outputs(3).path = scatter_fig_path;
end 

write_metadata(Metadata, [output_basename '.json']);