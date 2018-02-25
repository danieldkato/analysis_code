function T = get_trial_amplitudes(T, start_frame, end_frame)


for t = 1:length(T.Trials)
    
    dFF = T.Trials(t).dFF;
    baselines = mean(dFF(:,1:start_frame-1), 2);
    
    if nargin < 2
        abs_peaks = max(dFF(:, start_frame:end_frame), [], 2);
    else
        abs_peaks = max(dFF(:, start_frame:end_frame), [], 2);
    end
    
    T.Trials(t).amplitudes = abs_peaks - baselines;
    
end