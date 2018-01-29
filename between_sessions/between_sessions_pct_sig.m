C_struct = loadjson('/mnt/nas2/homes/dan/MultiSens/analysis/condition_settings.json');
Conditions = cell2mat(C_struct.conditions);

path1 = '/mnt/nas2/homes/dan/MultiSens/data/5037-1/2P/2017-12-18/site2/grab01/2P/segmentation/manual/manual3/plots4/within_session_regression.mat';
path2 = '/mnt/nas2/homes/dan/MultiSens/data/5037-1/2P/2017-12-31/site2/grab01/2P/segmentation/manual/manual1/csv_merge/within_session_analysis/within_session_regression.mat';


t1_idx = find(arrayfun(@(x) strcmp(x.name, 'test tone only'), Conditions));
t2_idx = find(arrayfun(@(x) strcmp(x.name, 'control tone only'), Conditions));
w_t1_idx = find(arrayfun(@(x) strcmp(x.name, 'control tone only'), Conditions));

C1 = load(path1);
C2 = load(path2);

Session(1).Conditions = C1.within_session_regressions.Conditions;
Session(2).Conditions = C2.within_session_regressions.Conditions;


% compare fraction with sig responses to tone 1
Pre_T1 = get_condition('test tone only', Session(1).Conditions);
Post_T1 = get_condition('test tone only', Session(2).Conditions);
Pre_T1_pvals = Pre_T1.pvals;
Post_T1_pvals = Post_T1.pvals;
x1 = sum(Pre_T1_pvals < 0.05);
x2 = sum(Post_T1_pvals < 0.05);
n1 = length(Pre_T1_pvals);
n2 = length(Post_T1_pvals);
p_hat = (x1 + x2)/(n1 + n2);
s = sqrt(p_hat*(1-p_hat)*(1/n1+1/n2));
z = ((x1/n1) + (x2/n2)) / s;
pvals(1) = normcdf(z)
pre_color = [162 180 232]/255;
post_color = [33 71 239]/255;
bar(categorical({'T1 pre'}), (x1/n1)*100, 'FaceColor', pre_color, 'EdgeColor', 'none');
hold on;
bar(categorical({'T1 post'}), (x2/n2)*100, 'FaceColor', post_color, 'EdgeColor', 'none');


% compare fraction with sig responses to tone 2
Pre_T2 = get_condition('control tone only', Session(1).Conditions);
Post_T2 = get_condition('control tone only', Session(2).Conditions);
Pre_T2_pvals = Pre_T2.pvals;
Post_T2_pvals = Post_T2.pvals;
x1 = sum(Pre_T2_pvals < 0.05);
x2 = sum(Post_T2_pvals < 0.05);
n1 = length(Pre_T2_pvals);
n2 = length(Post_T2_pvals);
p_hat = (x1 + x2)/(n1 + n2);
s = sqrt(p_hat*(1-p_hat)*(1/n1+1/n2));
z = ((x1/n1) + (x2/n2)) / s;
pvals(2) = normcdf(z)
pre_color = [177 226 182]/255;
post_color = [24 165 24]/255;
bar(categorical({'T2 pre'}), (x1/n1)*100, 'FaceColor', pre_color, 'EdgeColor', 'none');
hold on;
bar(categorical({'T2 post'}), (x2/n2)*100, 'FaceColor', post_color, 'EdgeColor', 'none');


% compare fraction with sig responses to conjunction
Pre_T3 = get_condition('stepper and test tone', Session(1).Conditions);
Post_T3 = get_condition('stepper and test tone', Session(2).Conditions);
Pre_T3_pvals = Pre_T3.pvals;
Post_T3_pvals = Post_T3.pvals;
x1 = sum(Pre_T3_pvals < 0.05);
x2 = sum(Post_T3_pvals < 0.05);
n1 = length(Pre_T3_pvals);
n2 = length(Post_T3_pvals);
p_hat = (x1 + x2)/(n1 + n2);
s = sqrt(p_hat*(1-p_hat)*(1/n1+1/n2));
z = ((x1/n1) + (x2/n2)) / s;
pvals(3) = normcdf(z)
pre_color = [160 141 196]/255;
post_color = [97 56 193]/255;
bar(categorical({'W+T1 pre'}), (x1/n1)*100, 'FaceColor', pre_color, 'EdgeColor', 'none');
hold on;
bar(categorical({'W+T1 post'}), (x2/n2)*100, 'FaceColor', post_color, 'EdgeColor', 'none');

title('Percentage neurons with significant responses, pre vs. post');
ylabel('% neurons with p-val < 0.05');