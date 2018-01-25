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
title(fig_title);
num_conditions = length(Conditions);


%% Plot each histogram:
for c = 1:num_conditions
    subplot(num_conditions, 1, c);
    h(c).handle = histogram(Conditions(c).distribution, 'FaceColor', Conditions(c).color);
    h(c).axes = gca;
    h(c).xlim = xlim;
    h(c).ylim = ylim;
    h(c).nbins = h(c).handle.NumBins;
    h(c).bin_edges = h(c).handle.BinEdges;
    title(Conditions(c).name);
    ylabel('Count');
    if c == num_conditions
        xlabel('dF/F (a.u.)');
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

