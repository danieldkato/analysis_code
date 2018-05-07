function [ttest_results, hist_fig_handle, regression_results, scatter_handle] = compare_2_conditions(Condition1, Condition2, field_name, paired, fig_title)

Conditions(1) = Condition1;
Conditions(2) = Condition2;


%% Compute means:
for c = 1:2
    Conditions(c).mean = mean(Conditions(c).(field_name));
end


%% Do t-test:
if paired
    [h, p, ci, stats] = ttest(Conditions(1).(field_name), Conditions(2).(field_name));
    disp(p);
else
    [h, p, ci, stats] = ttest2(Conditions(1).(field_name), Conditions(2).(field_name));
end

ttest_results.means = [Conditions.mean];
ttest_results.h = h;
ttest_results.p = p;
ttest_results.ci = ci;
ttest_results.stats = stats;


%% Plot histograms:
hist_fig_handle = figure();
hold on;
for c = 1:length(Conditions)
    H(c).handle = histogram(Conditions(c).(field_name), 'FaceColor', Conditions(c).color, 'EdgeColor', 'none');
    H(c).num_bins = H(c).handle.NumBins;
    H(c).bin_edges = H(c).handle.BinEdges;
end

% Go back through histograms and make bin edges consistent:
disp('[H.num_bins] = ');
disp([H.num_bins]);
[m, I] = max([H.num_bins]);
disp('I');
disp(I);
bin_edges = H(I).bin_edges;
for c = 1:length(H)
    H(c).handle.BinEdges = bin_edges;
end

pause(2);

% Go back through histograms and plot means:
yl = ylim;
disp('ylim 1');
disp(ylim)
for c = 1:length(H)
    condition_mean = mean(Conditions(c).(field_name));
    line([condition_mean condition_mean], [yl(1) yl(2)], 'Color', Conditions(c).color);
end

disp('ylim 2');
disp(ylim)
% Create title, labels, and annotation:
ylabel('Count');
if strcmp(field_name, 'amplitudes')
    xlabel('Peak dF/F (a.u.)');
elseif strcmp(field_name, 'Mean')
    xlabel('Mean dF/F (a.u.)');
end
pos = get(gca, 'Position');

disp('ylim 3');
disp(ylim)

legend(Conditions(1).abbreviation, Conditions(2).abbreviation);
legend('boxoff');


if nargin > 3
    title(fig_title, 'FontWeight', 'normal');
end


%% If the condition data are paired, do a regression and create a scatter plot:

if paired
    % Do a linear regression to get the slope:
    Tbl = table();
    Tbl.X = Condition2.(field_name);
    Tbl.Y = Condition1.(field_name) - Condition2.(field_name);
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
    Tbl2.X = Condition2.(field_name);
    Tbl2.Y = Condition1.(field_name) - Condition2.(field_name);
    lm2 = fitlm(Tbl2, 'Y ~ X');    
    regression_results.p = lm2.Coefficients{2, 4};
    
    % Create a scatterplot of each cell's response to both conditions:
    %scatter_fig_title = {['Peak dF/F response to \color[rgb]{' num2str(Condition1.color) '}' Condition1.name '\color[rgb]{0 0 0} vs \color[rgb]{' num2str(Condition2.color) '}' Condition2.name]; ['\color[rgb]{0 0 0} \fontsize{10}Mouse ' mouse]; ['\fontsize{10}Session ' date]};
    scatter_handle = scatter_conditions(Condition1, Condition2, field_name, fig_title); % Create scatter plot of W+T1 vs W
    lims = [xlim ylim];
    lowest = min(lims);
    highest = max(lims);
    xlim([lowest highest]);
    ylim([lowest highest]);
    hline = refline(1, 0);
    hline.Color = [0.40 0.40 0.40];
    lsline; % add least-squares line
end