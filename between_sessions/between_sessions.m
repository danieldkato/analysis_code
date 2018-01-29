C_struct = loadjson('/mnt/nas2/homes/dan/MultiSens/analysis/condition_settings.json');
Conditions = cell2mat(C_struct.conditions);

Sessions(1).params_path = 'pre2.json';
Sessions(2).params_path = 'post1.json';

for s = 1:length(Sessions)
    Sessions(s).params = loadjson(Sessions(s).params_path);
    
    T = trialize_data(Sessions(s).params.rawF_path, ... 
        Sessions(s).params.galvo_path, ...
        Sessions(s).params.timer_path, ...
        Sessions(s).params.ardu_path, ...
        Sessions(s).params.conditions_path, ...
        Sessions(s).params.grab_metadata, ...
        Sessions(s).params.pre_sec, ...
        Sessions(s).params.post_sec, ...
        Sessions(s).params.show_inflection_points);    
    
    Conditions = split_trials_by_condition(T.Trials, Conditions); % Split trials by condition
    Conditions = get_condition_means(Conditions); % Get the mean peristimulus dF/F trace of every neuron for every condition
    
    for c = 1:length(Conditions)
        Conditions(c).distribution = mean(Conditions(c).Mean(:, T.pre_frames + 1:end), 2);
    end
    
    Sessions(s).Conditions = Conditions;
    Sessions(s).Coeffs = C;
    
    % hack:
    grab_metadata = loadjson(Sessions(s).params.grab_metadata);
    mouse = grab_metadata.mouse;
end

output_directory = ['/mnt/nas2/homes/dan/MultiSens/data/' mouse '/2P/between_sessions_analysis/'];
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

t1_idx = find(arrayfun(@(x) strcmp(x.name, 'test tone only'), Conditions));
t2_idx = find(arrayfun(@(x) strcmp(x.name, 'control tone only'), Conditions));
w_t1_idx = find(arrayfun(@(x) strcmp(x.name, 'stepper and test tone'), Conditions));

% Compare distributions:

% Compare test tone:
Pre_pairing_test_tone = get_condition('test tone only', Sessions(1).Conditions);
Post_pairing_test_tone = get_condition('test tone only', Sessions(2).Conditions);
fig1_title = ['\color[rgb]{' num2str(Conditions(t1_idx).color) '}' Conditions(t1_idx).name '\color[rgb]{0 0 0} pre vs. post'];
Pre_pairing_test_tone.abbreviation = 'T1 pre';
Post_pairing_test_tone.abbreviation = 'T1 post';
Pre_pairing_test_tone.color = [0.5 0.5 0.5];
%Post_pairing_test_tone.color = 'green';
[f1, p1, stats1] = ttest_and_histogram(Pre_pairing_test_tone, Post_pairing_test_tone, false, fig1_title);
figname1 = [output_directory filesep 'T1_pre_vs_post.fig'];
savefig(f1, figname1);
S.p = p1;
S.stats = stats1;
fname1 = [output_directory filesep 'T1_pre_vs_post_stats.mat'];
save(fname1, 'S');

% Compare control tone:
Pre_pairing_control_tone = get_condition('control tone only', Sessions(1).Conditions);
Post_pairing_control_tone = get_condition('control tone only', Sessions(2).Conditions);
fig2_title = ['\color[rgb]{' num2str(Conditions(t2_idx).color) '}' Conditions(t2_idx).name '\color[rgb]{0 0 0} pre vs. post'];
Pre_pairing_control_tone.abbreviation = 'T2 pre';
Post_pairing_control_tone.abbreviation = 'T2 post';
Pre_pairing_control_tone.color = [0.5 0.5 0.5];
%Post_pairing_control_tone.color = 'green';
[f2, p2, stats2] = ttest_and_histogram(Pre_pairing_control_tone, Post_pairing_control_tone, false, fig2_title);
figname2 = [output_directory filesep 'T2_pre_vs_post.fig'];
savefig(f2, figname2);
S.p = p2;
S.stats = stats2;
fname2 = [output_directory filesep 'T2_pre_vs_post_stats.mat'];
save(fname2, 'S');

% Compare conjunction:
Pre_pairing_conjunction = get_condition('stepper and test tone', Sessions(1).Conditions);
Post_pairing_conjunction = get_condition('stepper and test tone', Sessions(2).Conditions);
fig3_title = ['\color[rgb]{' num2str(Conditions(w_t1_idx).color) '}' Conditions(w_t1_idx).name '\color[rgb]{0 0 0} pre vs. post'];
Pre_pairing_conjunction.abbreviation = 'W+T1 pre';
Post_pairing_conjunction.abbreviation = 'W+T1 post';
Pre_pairing_conjunction.color = [0.5 0.5 0.5];
%Post_pairing_conjunction.color = 'green';
[f3, p3, stats3] = ttest_and_histogram(Pre_pairing_conjunction, Post_pairing_conjunction, false, fig3_title);
figname3 = [output_directory filesep 'W_T1_pre_vs_post.fig'];
savefig(f3, figname3);
S.p = p3;
S.stats = stats3;
fname3 = [output_directory filesep 'W_T1_pre_vs_post_stats.mat'];
save(fname3, 'S');

% Compare fraction cells with significan regression coefficients:
%[C1, N1] = regress_neurons_v_trial_params();


% Write metadata:
Metadata.inputs(1).path = Sessions(1).params.rawF_path;
Metadata.inputs(2).path = Sessions(1).params.galvo_path;
Metadata.inputs(3).path = Sessions(1).params.timer_path;
Metadata.inputs(4).path = Sessions(1).params.ardu_path;
Metadata.inputs(5).path = Sessions(1).params.conditions_path;
Metadata.inputs(6).path = Sessions(1).params.grab_metadata;
Metadata.inputs(7).path = Sessions(2).params.rawF_path;
Metadata.inputs(8).path = Sessions(2).params.galvo_path;
Metadata.inputs(9).path = Sessions(2).params.timer_path;
Metadata.inputs(10).path = Sessions(2).params.ardu_path;
Metadata.inputs(11).path = Sessions(2).params.conditions_path;
Metadata.inputs(12).path = Sessions(2).params.grab_metadata;
Metadata.params.pre = Sessions(1).params.pre_sec; % hack
Metadata.params.post = Sessions(1).params.post_sec; % hack
Metadata.outputs(1).path = figname1;
Metadata.outputs(2).path = figname2;
Metadata.outputs(3).path = figname3;
Metadata.outputs(4).path = fname1;
Metadata.outputs(5).path = fname2;
Metadata.outputs(6).path = fname3;
write_metadata(Metadata, [output_directory filesep 'between_sessions.json']);