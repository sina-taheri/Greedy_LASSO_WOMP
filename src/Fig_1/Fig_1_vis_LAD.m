close all;
clear;
clc;

addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\data\Fig_1
addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\src\utils

WOMP_type = 'ladlassol1';
switch WOMP_type
    case 'ladlassol0'
        load('Fig_1_13_ladlassol0_N_300_s_10_m_150_n_iter_30');
        title_txt = '$\ell^0$-based LAD-LASSO-WOMP';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_1\\Fig_1_13_ladlassol0_N_300_s_10_m_150_n_iter_30';
    case 'ladlassol1'
        load('Fig_1_23_ladlassol1_N_300_s_10_m_150_n_iter_30');
        title_txt = '$\ell^1$-based LAD-LASSO-WOMP';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_1\\Fig_1_23_ladlassol1_N_300_s_10_m_150_n_iter_30';
end

%% Display - error vs. lambda
color_profile = [0 0 0; 0 0 1; 1 0 0; 0 1 0.5; 0 1 1; 1 0 1; 1 1 0];

% data structure
data.Y = lambda_err;
data.X = lambda_vals;
data.nX = n_lambda;
data.PLOT = beta;
data.nPLOT = n_beta;

% params structure
for i = 1:data.nPLOT
    legend_txt{i} = sprintf('$K = %.2f m$', beta(i));
end

params = set_plot_params(1, title_txt, legend_txt);

plot_error_lambda(data, color_profile, params)

% print(im_save_txt, '-depsc', '-opengl');
% saveas(gcf, im_save_txt, 'epsc');

%{
figure;
color_profile = [0 0 0; 0 0 1; 1 0 0; 0 1 0.5; 0 1 1; 1 0 1; 1 1 0];
% ax = axes('Position', [0.08 0.093 0.90 0.85]);
ax = axes;
for i = 1:n_beta
%     lambda_err1 = reshape(lambda_err(:, i, :), n_trial, n_lambda);
%     figure;
%     boxplot(lambda_err1, lambda_vals, 'color', color_profile(i, :));
    hbp = boxplot(lambda_err{i}, lambda_vals, 'color', color_profile(i, :));
    h = findobj(hbp, 'tag', 'Outliers');
    for iH = 1:length(h)
        h(iH).MarkerEdgeColor = color_profile(i, :);
    end
    hAx = gca;
    xtk = hAx.XTick; 
    hold on;
    lgd_txt = sprintf('%d %% corruption', 100*beta(i));
%     plot(xtk, median(lambda_err1, 1), 'Marker', '*', 'MarkerSize', 8, 'color', color_profile(i, :), 'LineWidth', 1.5, 'DisplayName', lgd_txt);
    plot(xtk, median(lambda_err{i}, 1), 'Marker', '*', 'MarkerSize', 8, 'color', color_profile(i, :), 'LineWidth', 1.5, 'DisplayName', lgd_txt);
    hold on;
end

str = sprintf('%.2e', lambda_vals);
str_mat = reshape(str, [8 n_lambda])';
num_parts = str2double(split(string(str_mat), 'e'));
xtk = cell(1, n_lambda);
for i = 1:n_lambda
    a = num_parts(i, 1);
    if a == 1
        b = num_parts(i, 2);
        xtk{i} = sprintf('$10^{%d}$', b);
    end
end
ax.XTickLabels = xtk;
ax.TickLabelInterpreter = 'latex';
set(gca, 'YScale', 'log', 'FontSize', 16);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 1]);
xlabel('$\lambda$', 'Fontsize', 18, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('$||x - x^\sharp||_2/||x||_2$', 'Fontsize', 20, 'FontWeight', 'bold', 'Interpreter', 'latex');
lgd = legend('show', 'Location', 'northwest');
lgd.FontSize = 20;
lgd.FontName = 'Times';
lgd.Interpreter = 'latex';
grid on;
grid minor;
title(title_txt, 'FontName', 'Times', 'Fontsize', 20, 'Interpreter', 'latex');
axis([1 n_lambda 10^-4.5 10]);
%}