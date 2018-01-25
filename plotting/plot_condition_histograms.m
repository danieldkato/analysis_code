function hist_figure = plot_condition_histograms(Conditions, fig_title)

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

