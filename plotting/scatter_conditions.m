function h = scatter_conditions(y_cond_name, x_cond_name, Conditions, field_name, fig_title)


% Get field name
if nargin > 3
    field_name = field_name;
else
    field_name = 'distribution';
end

% Get requested condition info:
X = get_condition(x_cond_name, Conditions);
Y = get_condition(y_cond_name, Conditions);

% Get the number of ROIs for each condition, and validate that they are the same:
if size(X.distribution, 1) == size(Y.(field_name), 1)
    num_ROIs = size(X.(field_name), 1);
else
    error(['Condition data for ' y_cond_name ' and ' x_cond_name ' include observations from different numbers of ROIs. Please verify that data are formatted correctly.']);
end

% Define figure size:
width = 400;
height = width; % Both axes measure the same quantity so I think it's clearest if the figure is square

% Create scatterplot:
h = figure;
scatter(X.(field_name), Y.(field_name), 80, '.k');

% Add axis labels:
X.label = xlabel([X.abbreviation ' peak dF/F (a.u.)'], 'FontSize', 11);
Y.label = ylabel([Y.abbreviation ' peak dF/F (a.u.)'], 'FontSize', 11);

% Add annotation with number of ROIs:
annotation('textbox', [0.6 0.2 0.2 0.2], 'String', ['n = ' num2str(num_ROIs)], 'LineStyle', 'none');

% Add unity line for reference:
hline = refline(1, 0);
hline.Color = [0.50 0.50 0.50];

% Adjust figure size:
curr_pos = get(gcf, 'position');
set(gcf, 'units', 'points', 'position', [curr_pos(1) curr_pos(2) width height]);

% Color-code axis labels if possible:
X_fields = fieldnames(X);
defines_color = find(cellfun(@(c) strcmp(c, 'color'), X_fields));
if defines_color
    X.label.Color = X.color;
end
Y_fields = fieldnames(Y);
defines_color = find(cellfun(@(c) strcmp(c, 'color'), Y_fields));
if defines_color
    Y.label.Color = Y.color;
end

% Add title: 
if nargin > 4
    title(fig_title, 'FontWeight', 'normal');
end