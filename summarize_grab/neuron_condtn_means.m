function M = neuron_condtn_means(N, Conditions)


M.frame_rate = N.frame_rate;
M.pre_frames = N.pre_frames;
M.post_frames = N.post_frames;
peri_stim_window = M.pre_frames + M.post_frames;

for n = 1:length(N.Neurons)
    
    Neuron = N.Neurons(n);
    Trials = Neuron.Trials;
    
    % Find which trials belong to each condition:
    Conditions = match_trials_to_conditions(Conditions, Trials);
    
    for c = 1:length(Conditions)
        curr_cond_trials = Conditions.;
    end
    
    
end
