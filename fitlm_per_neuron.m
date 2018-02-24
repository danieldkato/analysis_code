function [Neurons, Regressors] = fitlm_per_neuron(T, modelspec)
%% Identify regressors and validate modelspec and Trials:

Trials = T.Trials;

split = strsplit(modelspec);
term_indices = cellfun(@(c) ~isempty(regexp(c, '\w*', 'ONCE')), split); % find all 'words' (i.e., elements containing alphanumeric characters) in modelspec
terms = split(term_indices);

% Do some validation of modelspec to make sure it's formatted correctly:
if term_indices(1) ~= 1 || ~strcmp(split{2}, '~') || length(find(term_indices)) < 2
    error(['Modelspec ' modelspec ' is not formatted correctly; modelspec should include one response variable on the left hand side of a tidla followed by one or more regressor terms on the right-hand side.']);
elseif term_indices(end) ~= 1
    error(['Modelspec ' modelspec ' is not formatted correctly; last character must be part of a regressor term.']);
end

% Split response term from regressor terms:
response_var = terms{1};
regressors = terms(2:end);

% Exclude any offset terms:
contains_letters = cellfun(@(c) ~isempty(regexp(c, '[a-zA-Z]*', 'ONCE')), regressors);
regressors = regressors(contains_letters);

% Get just the basic regressors (i.e., no interactive terms)
no_nonwords = cellfun(@(c) isempty(regexp(c, '\W*', 'ONCE')), regressors);
basic_regressors = terms(no_nonwords);

% Validate that each basic regressor has a corresponding field in Trials
% with the same name:
Trials_fieldnames = fieldnames(Trials);
for r = 1:length(basic_regressors)
    rname = basic_regressors{r};
    matches = contains(Trials_fieldnames, rname);
    if isempty(find(matches, 1))
        error(['Input Trials struct does not contain field ' rname ' specified in modelspec. Please make sure that modelspec is specified correctly and that input Trials struct indcludes fields corresponding to each basic regressor term in modelspec.']);
    end
end


%% Set up table with regressor information:
Tbl = table();
for r = 1:length(basic_regressors)
    rname = basic_regressors{r};
    Tbl.(rname) = [Trials.(rname)]'; 
end


%% Do regressions for each neuron:

% Validate that all trials include observations from the same number of neurons:
n_ROIs_each_trial = arrayfun(@(x) size(x.(response_var), 1), T.Trials);
if range(n_ROIs_each_trial) == 0
    n_ROIs = n_ROIs_each_trial(1);
else
    error('Not all trials include observations from the same number of ROIs. Please check integrity of input data');
end

for n = 1:n_ROIs 
    
    disp(['Fitting linear model for ROI ' num2str(n) ' out of ' num2str(num_ROIs)]);
    
    % Get data just for current neuron:
    C = get_neuron_from_trials(T, n);
    curr_neuron_trials = C.Trials;
    
    % Put response variable for current neuron into table:
    Tbl.(response_var) = [curr_neuron_trials.(response_var)]';
    
    % Fit linear model:
    Neurons(n).lm = fitlm(Tbl, modelspec);
end


%% Compute the percentage of cells that have significant coefficients for each regressor:
for r = 1:length(regressors)
    Regressors(r).regressor_name = regressors{r};
    pvals = arrayfun(@(x) x.lm.Coefficients{r, 4}, Neurons);
    Regressors(r).num_sig_pvals = sum(pvals < 0.05);
    Regressors(r).pct_sig_pvals = (Regressors(r).num_sig_pvals/num_ROIs)*100;
end