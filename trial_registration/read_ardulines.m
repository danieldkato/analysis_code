function T = read_ardulines(ardulines_path)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-12    


%% I. OVERVIEW: 
% This function takes a .txt file containing serial output read from an
% Arduino running a given ArduFSM protocol for a single behavior session,
% and returns a T x 1 array of structs describing each trial delivered
% during the session, where T is the number of trials delivered during the
% session.


%% II. REQUIREMENTS:
% 1) The MATLAB toolbox JSONlab, available at https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files


%% III. INPUTS: 
% 1) txt - a path to a .txt file containing serial output from an Arduino
% running an ArduFSM protocol for a single behavior session. In accordance
% with the general ArduFSM framework, each line of output either
% acknowledges the receipt of instructions from the host PC, asserts
% upcoming trial parameters, reports recorded behavior parameters, or
% signals the start of a trial. More information about the ArduFSM
% framework can be found at https://github.com/cxrodgers/ArduFSM. 


%% IV. OUTPUTS:
% 1) Trials - a T x 1 array of structs describing the trial type of each
% trial delivered over the course of a behavior session, where T is the
% number of trials. Each element includes fields corresponding to various
% trial parameters as well as a concise, human-readable description of the
% trial condition.


%%

    
    % Load data into an k x 1 cell array, where k is the number of lines:
    linesFID = fopen(ardulines_path);
    txt = textscan(linesFID, '%s', 'Delimiter', '\r\n');
    fclose(linesFID);
    txt = txt{1,1};
    txt = txt(~cellfun(@isempty, txt)); % For some reason txt includes some empty lines; not sure why, but strip these
    
    % (Might eventually want to include some code here to ensure that each line begins with a number)
    
    % Find all trial start lines:
    trial_start_lines = find(~cellfun(@isempty, regexp(txt, 'TRL_[0-9]*_START')));
    
    % Determine the condition of each trial:
    for t = 1:length(trial_start_lines)
        
        % Get all lines associated with the current trial:
        if t < length(trial_start_lines)
            this_trial_last_line = trial_start_lines(t+1) - 1;
        elseif t == length(trial_start_lines)
            this_trial_last_line = length(txt);
        end
        this_trial_lines = txt(trial_start_lines(t):this_trial_last_line);
        
        % Extract the lines that define the current trial parameters:
        trial_param_lines = this_trial_lines(~cellfun(@isempty, regexp(this_trial_lines, 'TRLP')));
        
        % For each current trial parameter, get its name and value:
        for tp = 1:length(trial_param_lines)
            [start_idx, pName_end_idx] = regexp(trial_param_lines{tp}, 'TRLP [A-Z]+'); % get index of where parameter name ends
            param_name = trial_param_lines{tp}(start_idx+5:pName_end_idx); % get parameter name
            param_val = trial_param_lines{tp}(pName_end_idx+2:end); % get parameter value
            T(t).(param_name) = str2double(param_val);
        end
        
    end
end