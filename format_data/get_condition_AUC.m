function Conditions = get_condition_AUC(Conditions, start, stop)

for c = 1:length(Conditions)
    
    mean_traces = Conditions(c).Mean;
    
    if nargin < 2
        AUC = sum(mean_traces(:, start:stop), 2);
    else
        AUC = sum(mean_traces(:, start:end), 2);
    end
    
    Conditions(c).AUC = AUC;
end