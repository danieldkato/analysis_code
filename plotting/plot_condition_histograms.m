function hist_figure = plot_condition_histograms(Conditions, field_name, fig_title)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-29 


%% I. OVERVIEW: 
% This function plots histograms of some distribution associated with each
% trial/stimulus condition within a single imaging session. 


%% II. REQUIREMENTS:
% 1) MATLAB v >= ???


%% III. INPUTS: 
% 1) Conditions - c x 1 array of structs, where c is the number of
%    conditions being analyzed. Each element of Conditions should minimally
%    include the following fields:
%
%       Name - char array specifying human-readable condition name 
%
%       Color - 1 x 3 RGB vector specifying the color code for the
%       corresponding condition.
%
%       In addition, each element of `Conditions` should have a field
%       consisting of a d x 1 vector to be plotted in a histogram, where d
%       is the number of observations in the distribution. If the
%       `field_name` input argument is specified, then the name of the
%       field to be plotted should be specified by `field_name`. If
%       `field_name` is not specified, then this function will look for a
%       field simply called `distribution`.
%
% 2) field_name (optional)- char vector speciying the name of the field for
%    which a histogram should be plotted for each condition. If this input
%    argument is not specified, this function will look for a field simply
%    called `distribution`.
%
% 3) fig_title (optional)- char vector or cell specifying the figure title.


%% IV. OUTPUTS:
% 1) mean_figure - handle to created figure


%% Get the name of the field to plot:
if nargin > 1
    field_name = field_name;
elseif nargin == 1 || isempty(field_name)
    field_name = distribution;
end


%% Create figure:
hist_figure = figure();
num_conditions = length(Conditions);


%% Plot each histogram:
for c = 1:num_conditions
    subplot(num_conditions, 1, c);
    h(c).handle = histogram(Conditions(c).(field_name), 'FaceColor', Conditions(c).color, 'EdgeColor', 'none');
    h(c).axes = gca;
    h(c).xlim = xlim;
    h(c).ylim = ylim;
    h(c).nbins = h(c).handle.NumBins;
    h(c).bin_edges = h(c).handle.BinEdges;
    pos = get(gca, 'Position');
    h(c).width = pos(3);
    h(c).height = pos(4);
    
    % Create title and axis labels:
    %title(Conditions(c).abbreviation, 'FontWeight', 'normal', 'FontSize', 10, 'Color', Conditions(c).color);
    ylabel('Count', 'FontSize', 10);
    if c == num_conditions
        xlabel('dF/F (a.u.)', 'FontSize', 10);
    end
end


%% Make axes, bin edges consistent:

% Get axis limits in data units:
max_x = max([h.xlim]);
min_x = min([h.xlim]);
max_y = max([h.ylim]);
min_y = min([h.ylim]);

% Get bin edges of histogram with greatest number of bins:
[m, I] = max([h(c).nbins]);
bin_edges = h(I).bin_edges;

% Get max height and width of all axis objects:
max_width = max([h.width]);
max_height = max([h.height]);

cumulative_y_shift = 0;
for c = 1:length(Conditions)
    axes(h(c).axes);
    xlim([min_x max_x]);
    ylim([min_y max_y]);
    h(c).handle.BinEdges = bin_edges;
    pos = get(gca, 'Position');
    cumulative_y_shift = cumulative_y_shift + max_height - pos(4);
    set(gca, 'Position', [pos(1), pos(2)-cumulative_y_shift, max_width, max_height]);
    
    % Add annotation with condition abbreviation and number of observations:
    num_observations = size(Conditions(c).(field_name), 1);
    pos = get(gca, 'Position');
    n_annotation_width = 0.1;
    n_annotation_height = 0.1;
    color_str = num2str(Conditions(c).color);
    cond_str = ['\color[rgb]{' color_str '}' Conditions(c).abbreviation '\color[rgb]{0 0 0}'];
    n_str = ['n = ' num2str(num_observations)];
    total_str = {cond_str; n_str};
    b = annotation('textbox', [pos(1)+pos(3)-n_annotation_width, pos(2)+pos(4)-n_annotation_height, n_annotation_width, n_annotation_height], 'String', total_str, 'LineStyle', 'none');
end

% Plot overall title if specified as input argument:
if nargin > 1
    text_box_width = 0.2; % normalized units
    text_box_height = 0.1;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2)-cumulative_y_shift pos(3) pos(4)*1.2 + cumulative_y_shift]);
    a = annotation('textbox', [(1-text_box_width)/2, 1-text_box_height , text_box_width, text_box_height], 'String', fig_title, 'HorizontalAlignment', 'center', 'FontSize', 11, 'LineStyle', 'none', 'FitBoxToText', 'on');
end
