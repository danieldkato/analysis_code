function Neurons = fitlm_per_neuron_within_session(params_file)
%% Load parameters and metadata:

% Load params from params file:
params = loadjson(params_file);
output_directory = params.output_directory;

% Load condition information like color codes, etc:
C_struct = loadjson(params.conditions_path);
Conditions = cell2mat(C_struct.conditions);

% Get grab metadata:
grab_metadata = loadjson(params.grab_metadata);
pre_stim_frames = ceil(params.pre_sec * grab_metadata.frame_rate);
mouse = grab_metadata.mouse;
date = grab_metadata.date;


%% Trialize data:
[Trials, Meta] = trialize_data(params.rawF_path, ... 
    params.galvo_path, ...
    params.timer_path, ...
    params.ardu_path, ...
    params.conditions_path, ...
    params.grab_metadata, ...
    params.pre_sec, ...
    params.post_sec, ...
    params.show_inflection_points);


%% Do regressions for each neuron:

% Initialize some variables that will be useful for regression:
Tbl = table();

w_idx = 1;
t1_idx = 2;
t2_idx = 3;
w_t1_idx = 4;
w_t2_idx = 5;

regressors{w_idx} = 'W';
regressors{t1_idx} = 'T1';
regressors{t2_idx} = 'T2';
regressors{w_t1_idx} = 'W:T1';
regressors{w_t2_idx} = 'W:T2';

model_spec = ['response ~ -1 + ' strjoin(regressors, ' + ')];
disp(model_spec);

% Get the number of ROIs in the dataset:
n_ROIs_each_trial = arrayfun(@(x) size(x.dFF, 1), Trials);
check_n_ROIs = circshift(n_ROIs_each_trial, 1);
if isequal(n_ROIs_each_trial, check_n_ROIs)
    num_ROIs = n_ROIs_each_trial(1);
else
    error('Not all trials include observations from the same number of ROIs. Please check integrity of input data');
end

% Regress peak dF/F response against trial parameters for each neuron:
for n = 1:num_ROIs
    
    disp(['Fitting linear model for ROI ' num2str(n) ' out of ' num2str(num_ROIs)]);
    
    % Get data just for current neuron:
    curr_neuron_trials = get_neuron_from_trials(Trials, n);
    
    % Assemble response and parameter data into table:
    Tbl.W = [curr_neuron_trials.STPRIDX]';
    Tbl.T1 = ([curr_neuron_trials.SPKRIDX] == 1)';
    Tbl.T2 = ([curr_neuron_trials.SPKRIDX] == 2)';
    trials_pre_stim_cropped = arrayfun(@(x) x.dFF(:, pre_stim_frames+1:end), curr_neuron_trials, 'UniformOutput', false)';
    trials_pre_stim_cropped = cell2mat(trials_pre_stim_cropped);
    disp(size(trials_pre_stim_cropped));
    Tbl.response = mean(vertcat(trials_pre_stim_cropped), 2);
    
    % Fit linear model:
    Neurons(n).lm = fitlm(Tbl, model_spec); 
    
    %{
    Neurons(n).w_pval = Neurons(n).lm.Coefficients{w_idx, 4};
    Neurons(n).t1_pval = Neurons(n).lm.Coefficients{t1_idx, 4};
    Neurons(n).t2_pval = Neurons(n).lm.Coefficients{t2_idx, 4};
    Neurons(n).w_t1_pval = Neurons(n).lm.Coefficients{w_t1_idx, 4};
    Neurons(n).w_t2_pval = Neurons(n).lm.Coefficients{w_t2_idx, 4};
    %}
    
end


%% Compute the percentage of cells that have significant coefficients for each regressor:
for r = 1:length(regressors)
    pvals = arrayfun(@(x) x.lm.Coefficients{r, 4}, Neurons);
    Regressors(r).num_sig_pvals = sum(pvals < 0.05);
    Regressors(r).pct_sig_pvals = (Regressors(r).num_sig_pvals/num_ROIs)*100;
    Regressors(r).regressor_name = regressors{r};
    
    % hack for plotting purposes later:
    Regressors(r).abbreviation = Conditions(r).abbreviation;
    Regressors(r).color = Conditions(r).color;
    
end

%{
pcts_table = table(R.pct_sig_pvals, 'VariableNames', regressors);
disp(pcts_table);
%}

% Write to CSV:
dat = cell(2, length(regressors));
dat(1,:) = regressors;
dat(2,:) = num2cell([Regressors.pct_sig_pvals]);
disp(dat);


%% Plot bar graph of percentages of neurons with significant coefficients for each regressor:
bar_fig = figure();
hold on;
for r = 1:length(Regressors)
    c = categorical({Regressors(r).abbreviation});
    bar(c, Regressors(r).pct_sig_pvals, 'FaceColor', Regressors(r).color, 'EdgeColor', 'none');
end
yl = ylim;
ylim([yl(1) yl(2)*1.1]);
title({'Percentage neurons with significant coefficients'; ['Mouse ' mouse ', session ' date]});


%% Save output to secondary storage:

% Create output directory if it doesn't exist and cd to it:
if ~exist(output_directory, 'dir')
    mkdir(output_directory)
end
old = cd(output_directory);

% Create dedicated sub-directory for output of this function and cd to it:
mkdir('fitlm_per_neuron');
cd('fitlm_per_neuron');

% Save stats:
fitlm_per_neuron_path = [output_directory filesep 'fitlm_per_neuron' filesep 'fitlm_per_neuron.mat'];
save(fitlm_per_neuron_path, 'Neurons', 'Regressors');

% Save figure:
bar_fig_path = [output_directory filesep 'fitlm_per_neuron' filesep 'pct_neurons_sig_coeffs.fig'];
savefig(bar_fig, bar_fig_path);


%% Save metadata:
Metadata.inputs(1).path = params.rawF_path;
Metadata.inputs(2).path = params.galvo_path;
Metadata.inputs(3).path = params.timer_path;
Metadata.inputs(4).path = params.ardu_path;
Metadata.inputs(5).path = params.conditions_path;
Metadata.inputs(6).path = params.grab_metadata;
Metadata.params.pre_onset_period = params.pre_sec;
Metadata.params.post_onset_period = params.post_sec;
Metadata.outputs(1).path = fitlm_per_neuron_path;
Metadata.outputs(2).path = bar_fig_path;
write_metadata(Metadata, 'fitlm_per_neuron_metadata.json');

cd(old);