function h = scatter_conditions(Condition_Y, Condition_X, field_name, fig_title)


% Get field name
if nargin > 3
    field_name = field_name;
else
    field_name = 'distribution';
end

% Get the number of ROIs for each condition, and validate that they are the same:
if size(Condition_X.(field_name), 1) == size(Condition_Y.(field_name), 1)
    num_ROIs = size(Condition_X.(field_name), 1);
else
    error(['Condition data for ' Condition_Y ' and ' Condition_X ' include observations from different numbers of ROIs. Please verify that data are formatted correctly.']);
end

% Define figure size:
width = 400;
height = width; % Both axes measure the same quantity so I think it's clearest if the figure is square

% Create scatterplot:
h = figure;
scatter(Condition_X.(field_name), Condition_Y.(field_name), 80, '.k');

% Add axis labels:
if strcmp(field_name, 'amplitudes')
    axis_description = ' peak dF/F (a.u.)';
elseif strcmp(field_name, 'Mean')
    axis_description = ' mean dF/F (a.u.)';
else
    axis_description = field_name;
end
Condition_X.label = xlabel([Condition_X.abbreviation axis_description], 'FontSize', 11);
Condition_Y.label = ylabel([Condition_Y.abbreviation axis_description], 'FontSize', 11);

% Add annotation with number of ROIs:
annotation('textbox', [0.6 0.2 0.2 0.2], 'String', ['n = ' num2str(num_ROIs)], 'LineStyle', 'none');

% Add unity line for reference:
hline = refline(1, 0);
hline.Color = [0.50 0.50 0.50];

% Adjust figure size:
curr_pos = get(gcf, 'position');
set(gcf, 'units', 'points', 'position', [curr_pos(1) curr_pos(2) width height]);

% Color-code axis labels if possible:
X_fields = fieldnames(Condition_X);
defines_color = find(cellfun(@(c) strcmp(c, 'color'), X_fields));
if defines_color
    Condition_X.label.Color = Condition_X.color;
end
Y_fields = fieldnames(Condition_Y);
defines_color = find(cellfun(@(c) strcmp(c, 'color'), Y_fields));
if defines_color
    Condition_Y.label.Color = Condition_Y.color;
end

% Add title: 
if nargin > 3
    title(fig_title, 'FontWeight', 'normal');
end