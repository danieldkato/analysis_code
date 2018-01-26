function stim_coeffs_within_session(params_file)
%% Load parameters and metadata:

% Load params from params file:
params = loadjson(params_file);
output_directory = params.output_directory;

% Load some useful imaging session meadata:
grab_metadata = loadjson(params.grab_metadata);
mouse = grab_metadata.mouse;
date = grab_metadata.date;

% Load condition information like color codes, etc:
C_struct = loadjson(params.conditions_path);
Conditions = C_struct.conditions;
Conditions = cell2mat(Conditions); % format as an array of structs, rather than cell array

% Compute some quantities that will be needed later on:
frame_rate = grab_metadata.frame_rate;
pre_frames = ceil(frame_rate * params.pre_sec);
post_frames = ceil(frame_rate * params.post_sec);


%% Trialize data:

T = trialize_data(params.rawF_path, ... 
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

regressor{w_idx} = 'whisker';
regressor{t1_idx} = 'tone_1';
regressor{t2_idx} = 'tone_2';
regressor{w_t1_idx} = 'whisker:tone_1';
regressor{w_t2_idx} = 'whisker:tone_2';

model_spec = ['response ~ -1 + ' strjoin(regressor, ' + ')];
disp(model_spec);

% Get the number of ROIs in the dataset:
n_ROIs_each_trial = arrayfun(@(x) size(x.dFF, 1), T.Trials);
check_n_ROIs = circshift(n_ROIs_each_trial, 1);
if isequal(n_ROIs_each_trial, check_n_ROIs)
    num_ROIs = n_ROIs_each_trial(1);
else
    error('Not all trials include observations from the same number of ROIs. Please check integrity of input data');
end

% Regress peak dF/F response against trial parameters for each neuron:
for n = 1:num_ROIs
    
    % Get data just for current neuron:
    C = get_neuron_from_trials(T, n);
    curr_neuron_trials = C.Trials;
    
    % Assemble response and parameter data into table:
    Tbl.whisker = [curr_neuron_trials.STPRIDX]';
    Tbl.tone_1 = ([curr_neuron_trials.SPKRIDX] == 1)';
    Tbl.tone_2 = ([curr_neuron_trials.SPKRIDX] == 2)';
    Tbl.response = max(vertcat(curr_neuron_trials.dFF), [], 2);
    
    % Fit linear model:
    Neurons(n).lm = fitlm(Tbl, model_spec); 
end

% Save linear models to secondary storage:
neurons_lms_path = ([output_directory filesep 'neurons_lms.mat']);
save(neurons_lms_path, 'Neurons');


%% Save metadata:
Metadata.inputs(1).path = params.rawF_path;
Metadata.inputs(2).path = params.galvo_path;
Metadata.inputs(3).path = params.timer_path;
Metadata.inputs(4).path = params.ardu_path;
Metadata.inputs(5).path = params.conditions_path;
Metadata.inputs(6).path = params.grab_metadata;
Metadata.params.pre_onset_period = params.pre_sec;
Metadata.params.post_onset_period = params.post_sec;
Metadata.outputs(1).path = neurons_lms_path;

write_metadata(Metadata, 'regress_neurons_v_trial_params.json');
