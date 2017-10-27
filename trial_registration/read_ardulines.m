% Last updated DDK 2016-09-27    

% OVERVIEW: 
% This function takes a .txt file containing serial output read from an
% Arduino running a given ArduFSM protocol for a single behavior session,
% and returns a T x 2 cell array of the trial type of each trial delivered
% during the session, where T is the number of trials delivered during the
% session.


% INPUTS: 
% 1) input - a path to a .txt file containing serial output from an Arduino
% running an ArduFSM protocol for a single behavior session. In accordance
% with the general ArduFSM framework, each line of output either
% acknowledges the receipt of instructions from the host PC, asserts
% upcoming trial parameters, reports recorded behavior parameters, or
% signals the start of a trial. More information about the ArduFSM
% framework can be found at: 

% https://github.com/cxrodgers/ArduFSM

% 2) conditionSettings - path to a .txt file contaning information about
% different trial conditions used in the current experiment. This text file
% should contain a series of MATLAB statements defining a c x 1 cell array
% of structs called Conditions, where c is the number of conditions. Each
% struct should minimally include fields corresponding to whatever trial
% parameters are necessary for defining a condition in the current
% experiment, as well as a concise condition name that will be assigned to
% each trial and that can be used later on to easily parse experiments by
% condition. Example contents of a conditionSettings.txt file might
% include:

% Conditions{1}.STPRIDX = 1;
% Conditions{1}.SPKRIDX = 0;
% Conditions{1].Name = ' stepper only ';

% Conditions{2}.STPRIDX = 0;
% Conditions{2}.SPKRIDX = 1;
% Conditions{2}.Name = ' speaker only ';


% OUTPUTS:
% 1) trialTypes - a T x 2 cell array of strings describing the trial type
% of each trial delivered over the course of a behavior session, where T is
% the number of trials. The first column contains the name of the trial
% condition, and the second column contains the trial duration. Trial types
% are listed in the order they are delivered.


% TODO: 
% This function is currently hard-coded to look for trial parameters and
% return trial types specific to the multiSens protocol. It may be
% worthwhile building more flexibility into this function. s

%%
function T = read_ardulines(txt, condition_settings)
    
    % Load data into an k x 1 cell array, where k is the number of lines:
    linesFID = fopen(txt);
    txt = textscan(linesFID, '%s', 'Delimiter', '\r\n');
    fclose(linesFID);
    txt = txt{1,1};
    txt = txt(~cellfun(@isempty, txt)); % For some reason txt includes some empty lines; not sure why, but strip these
    
    % Load a structure containing information about the conditions we're looking for:
    condsFID = fopen(condition_settings);
    content = fscanf(condsFID, '%c');
    eval(content);
    
    % (Might eventually want to include some code here to ensure that each line begins with a number)
    
    % Find all trial start lines:
    trial_start_lines = find(~cellfun(@isempty, regexp(txt, 'TRL_[0-9]*_START')));
    
    % Initialize a t x 2 cell array, where t is the number of trials, to be returned to the calling function:
    T = cell(length(trial_start_lines), 2);
    
    % Determine the condition of each trial:
    for t = 1:length(trial_start_lines)
        
        % Get all lines associated with the current trial:
        if t < length(trial_start_lines)
            this_trial_last_line = trial_start_lines(t+1) - 1;
        elseif t == length(trial_start_lines)
            this_trial_last_line = length(txt);
        end
        this_trial_lines = txt(trial_start_lines(t):this_trial_last_line);
        
        % Extract trial parameter lines for current trial:
        param_def_lines = this_trial_lines(~cellfun(@isempty, regexp(this_trial_lines, 'TRLP')));
        
        % Make a cell array of current trial parameters:
        this_trial_params = cell(length(param_def_lines), 2);
        for p = 1:length(this_trial_params)
           [trlp1, trlpEnd] = regexp(param_def_lines{p}, 'TRLP [A-Z]+');
           this_trial_params{p, 1} = param_def_lines{p}(trlp1+5:trlpEnd); % parameter name
           this_trial_params{p, 2} = str2double(param_def_lines{p}(trlpEnd+2:end)); % parameter value
        end 
        
        % For each possible condition, check if the defining parameters match the known parameters of the current trial:
        for c = 1:length(Conditions)
            
            n_params_compared = 0;
            n_params_matching = 0;
            
            field_names = fieldnames(Conditions{c});
            
            % Go through each parameter for the current possible condition under consideration...
            for n = 1:length(field_names)
                
                % ... and see if that parameter is specified for the current trial:
                matches = cellfun(@(c) strcmp(field_names{n}, c), this_trial_params(:,1));
                
                % If so:
                if any(matches) 
                    n_params_compared = n_params_compared+1; 
                    disp(eval(strcat('Conditions{',num2str(c),'}.',field_names{n})));
                    if this_trial_params{matches,2} == eval(strcat('Conditions{',num2str(c),'}.',field_names{n}))
                        n_params_matching = n_params_matching + 1;
                    end
                end
            end
            
            % If all compared parameters match:
            if n_params_matching/n_params_compared == 1
                
                % Write the name of the matched condition into T{t,1}, where t is the current trial number:
                T{t,1} = Conditions{c}.Name;
                
                % Try to get the trial duration:
                durInd = cellfun(@(c) strcmp('STIMDUR', c), this_trial_params(:,1));
                if any(durInd)
                    T{t,2} = this_trial_params{durInd, 2};
                end
                
                break
                
            end
            
        end
        
    end
end