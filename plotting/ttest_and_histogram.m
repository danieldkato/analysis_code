function [f, p, stats] = ttest_and_histogram(Condition1, Condition2, paired, fig_title)

Conditions(1) = Condition1;
Conditions(2) = Condition2;


%% Do t-test:
if paired
    [h, p, ci, stats] = ttest(Conditions(1).distribution, Conditions(2).distribution);
    disp(p);
else
    [h, p, ci, stats] = ttest2(Conditions(1).distribution, Conditions(2).distribution);
end

%% Plot histograms:
f = figure();
hold on;
for c = 1:length(Conditions)
    H(c).handle = histogram(Conditions(c).distribution, 'FaceColor', Conditions(c).color, 'EdgeColor', 'none');
    H(c).num_bins = H(c).handle.NumBins;
    H(c).bin_edges = H(c).handle.BinEdges;
end

% Go back and make bin edges consistent and plot means:
[m, I] = max(H(c).num_bins);
bin_edges = H(I).bin_edges;
yl = ylim;
for c = 1:length(H)
    H(c).handle.BinEdges = bin_edges;
    condition_mean = mean(Conditions(c).distribution);
    line([condition_mean condition_mean], [yl(1) yl(2)], 'Color', Conditions(c).color);
end


%% Create title, labels, and annotation:
%title_str = {'\color[rgb]{' num2str(Conditions(1).color) '}' Conditions(1).name '\color[rgb]{0 0 0} vs \color[rgb]{' num2str() '}'};
ylabel('count');
xlabel('dF/F response');
pos = get(gca, 'Position');
%t = text('textbox', [pos(1)+pos(3)*0.8,  pos(2)+pos(4)*0.8, 0.1, 0.1], 'String', {['p = ' num2str(p)]}, 'LineStyle', 'none');
legend(Conditions(1).abbreviation, Conditions(2).abbreviation);
legend('boxoff');

if nargin > 3
    title(fig_title);
end