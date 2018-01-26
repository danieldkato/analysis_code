function hist_figure = plot_condition_histograms(Conditions, fig_title)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-25 


%% I. OVERVIEW: 
% This function plots histograms of some distribution associated with each
% trial/stimulus condition within a single imaging session. 


%% II. REQUIREMENTS:
% 1) MATLAB v >= ???


%% III. INPUTS: 
% 1) Conditions - c x 1 array of structs, where c is the number of
%    conditions in which the current neuron was observed. Each element
%    of Conditions should minimally include the following fields:
%           
%       Distribution - 1 x d vector specifying some distributuon associated
%       with the corresponding condition, where d is the number of
%       observations in the distribution.
%
%       Name - char array specifying human-readable condition name 
%
%       Color - 1 x 3 RGB vector specifying the color code for the
%       corresponding condition.


%% IV. OUTPUTS:
% 1) mean_figure - handle to created figure


%% Create figure:
hist_figure = figure();
num_conditions = length(Conditions);


%% Plot each histogram:
for c = 1:num_conditions
    subplot(num_conditions, 1, c);
    h(c).handle = histogram(Conditions(c).distribution, 'FaceColor', Conditions(c).color, 'EdgeColor', 'none');
    h(c).axes = gca;
    h(c).xlim = xlim;
    h(c).ylim = ylim;
    h(c).nbins = h(c).handle.NumBins;
    h(c).bin_edges = h(c).handle.BinEdges;
    title(Conditions(c).abbreviation, 'FontWeight', 'normal', 'FontSize', 10, 'Color', Conditions(c).color);
    ylabel('Count', 'FontSize', 10);
    if c == num_conditions
        xlabel('dF/F (a.u.)', 'FontSize', 10);
    end
end


%% Make axes, bin edges consistent:
max_x = max([h.xlim]);
min_x = min([h.xlim]);

max_y = max([h.ylim]);
min_y = min([h.ylim]);

[m, I] = max([h(c).nbins]);
bin_edges = h(I).bin_edges;

for c = 1:length(h)
    axes(h(c).axes);
    xlim([min_x max_x]);
    ylim([min_y max_y]);
    h(c).handle.BinEdges = bin_edges;
end

% Plot overall title if specified as input argument:
if nargin > 1
    text_box_width = 0.2; % normalized units
    text_box_height = 0.1;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) pos(3) pos(4)*1.2]);
    a = annotation('textbox', [(1-text_box_width)/2, 1-text_box_height , text_box_width, text_box_height], 'String', fig_title, 'HorizontalAlignment', 'center', 'FontSize', 11, 'LineStyle', 'none', 'FitBoxToText', 'on');
end
