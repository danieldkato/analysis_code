function Condition = get_condition(name, Conditions)

% Verify that struct array Conditions has the field 'name':
fields = fieldnames(Conditions);
has_names = find(cellfun(@(c) strcmp(c, 'name'), fields),1);
if isempty(has_names)
    error('Input struct array Conditions has no ''name'' field. Please ensure that input struct array Conditions is formatted correctly.');
end

% Try to find the index of the unique element of Conditions whose 'name' field matches the specified input name: 
idx = find(arrayfun(@(x) strcmp(x.name, name), Conditions), 1);

if length(idx) == 1
    Condition = Conditions(idx);
elseif isempty(idx)
    Condition = [];
    warning(['Could not find condition named ' name '; please ensure that requested condition name was specified correctly and that data includes trial from requested condition.']);
elseif length(idx) > 1
    error(['More than one condition named ' name ' found in data; please ensure that each condition has a unique name.'])
end

