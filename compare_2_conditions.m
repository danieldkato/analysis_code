function [ttest_results, hist_fig_handle, regression_results, scatter_handle] = compare_2_conditions(Condition1, Condition2, field, paired, fig_title)

Conditions(1) = Condition1;
Conditions(2) = Condition2;


%% Do t-test:
if paired
    [h, p, ci, stats] = ttest(Conditions(1).(field), Conditions(2).(field));
    disp(p);
else
    [h, p, ci, stats] = ttest2(Conditions(1).(field), Conditions(2).(field));
end

ttest_results.h = h;
ttest_results.p = p;
ttest_results.ci = ci;
ttest_results.stats = stats;


%% Plot histograms:
hist_fig_handle = figure();
hold on;
for c = 1:length(Conditions)
    H(c).handle = histogram(Conditions(c).(field), 'FaceColor', Conditions(c).color, 'EdgeColor', 'none');
    H(c).num_bins = H(c).handle.NumBins;
    H(c).bin_edges = H(c).handle.BinEdges;
end

% Go back through histograms and a) make bin edges consistent and b) plot means:
[m, I] = max(H(c).num_bins);
bin_edges = H(I).bin_edges;
yl = ylim;
for c = 1:length(H)
    H(c).handle.BinEdges = bin_edges;
    condition_mean = mean(Conditions(c).(field));
    line([condition_mean condition_mean], [yl(1) yl(2)], 'Color', Conditions(c).color);
end

% Create title, labels, and annotation:
pos = get(gca, 'Position');
legend(Conditions(1).abbreviation, Conditions(2).abbreviation);
legend('boxoff');

if nargin > 3
    title(fig_title, 'FontWeight', 'normal');
end


%% If the condition data are paired, do a regression and create a scatter plot:

if paired
    % Do a linear regression to get the slope:
    Tbl = table();
    Tbl.X = Condition2.(field);
    Tbl.Y = Condition1.(field) - Condition2.(field);
    regression_results.lm = fitlm(Tbl, 'Y ~ X');    
    
    % Do a separate linear regression to get the p-value; we want to test
    % whether the slope of the best fit line for Condition1 vs Condition2
    % is significantly different from 1, but MATLAB's fitlm function tests
    % whether the slope of a best fit line is significantly different from
    % 0. So instead of regressing Condition1 on Condition2, we'll regress
    % (Condition1 - Condition2) on to condition Condition2; if Condition1 =
    % Condition2, then the slope of (Condition1-Condition2) vs Condition2
    % should be 0:
    Tbl2 = table();
    Tbl2.X = Condition2.(field);
    Tbl2.Y = Condition1.(field) - Condition2.(field);
    lm2 = fitlm(Tbl2, 'Y ~ X');    
    regression_results.p = lm2.Coefficients{2, 4};
    
    % Create a scatterplot of each cell's response to both conditions:
    %scatter_fig_title = {['Peak dF/F response to \color[rgb]{' num2str(Condition1.color) '}' Condition1.name '\color[rgb]{0 0 0} vs \color[rgb]{' num2str(Condition2.color) '}' Condition2.name]; ['\color[rgb]{0 0 0} \fontsize{10}Mouse ' mouse]; ['\fontsize{10}Session ' date]};
    scatter_handle = scatter_conditions(cond1_name, cond2_name, Conditions, field, fig_title); % Create scatter plot of W+T1 vs W
    lims = [xlim ylim];
    lowest = min(lims);
    highest = max(lims);
    xlim([lowest highest]);
    ylim([lowest highest]);
    hline = refline(1, 0);
    hline.Color = [0.40 0.40 0.40];
    lsline; % add least-squares line
end