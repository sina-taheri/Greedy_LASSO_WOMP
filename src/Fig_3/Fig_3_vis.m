close all;
clear;
clc;

addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\data\Fig_3

experiment = 'srlassol1_low';
switch experiment
    case 'lassol1_low'
        load('Fig_3_11_smaller_lassol1_N_500_s_10_m_40_n_iter_30');
        title_txt = '$\ell^1$-based LASSO-WOMP, $m = 40$';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_3\\Fig_3_11_smaller_lassol1_N_500_s_10_m_40_n_iter_30';
    case 'srlassol1_low'
        load('Fig_3_12_smaller_srlassol1_N_500_s_10_m_35_n_iter_30');
        title_txt = '$\ell^1$-based SR-LASSO-WOMP, $m = 35$';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_3\\Fig_3_12_smaller_srlassol1_N_500_s_10_m_35_n_iter_30';
    case 'ladlassol1_low'
        load('Fig_3_13_smaller_ladlassol1_N_500_s_10_m_60_n_iter_30');
        title_txt = '$\ell^1$-based LAD-LASSO-WOMP, $m = 60$';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_3\\Fig_3_13_smaller_ladlassol1_N_500_s_10_m_60_n_iter_30';
    case 'lassol1_high'
        load('Fig_3_21_smaller_lassol1_N_500_s_10_m_100_n_iter_30');
        title_txt = '$\ell^1$-based LASSO-WOMP, $m = 100$';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_3\\Fig_3_21_smaller_lassol1_N_500_s_10_m_100_n_iter_30';
    case 'srlassol1_high'
        load('Fig_3_22_smaller_srlassol1_N_500_s_10_m_100_n_iter_30');
        title_txt = '$\ell^1$-based SR-LASSO-WOMP, $m = 100$';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_3\\Fig_3_22_smaller_srlassol1_N_500_s_10_m_100_n_iter_30';
    case 'ladlassol1_high'
        load('Fig_3_23_smaller_ladlassol1_N_500_s_10_m_120_n_iter_30');
        title_txt = '$\ell^1$-based LAD-LASSO-WOMP, $m = 120$';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_3\\Fig_3_23_smaller_ladlassol1_N_500_s_10_m_120_n_iter_30';
end

%% Display - error vs. lambda
color_profile = [0 0 1; 1 0 0; 0 1 0.5; 0 0 0; 0 1 1; 1 0 1; 1 1 0];

% data structure
data.Y = lambda_err;
data.X = lambda_vals;
data.nX = n_lambda;
data.PLOT = w_vals;
data.nPLOT = n_w;

% params structure
for i = 1:data.nPLOT
    str = sprintf('%.2e', w_vals(i));
    num_parts = str2double(split(string(str), 'e'));
    if num_parts(1) == 1
        if num_parts(2) == 0
            legend_txt{i} = '$w_0 = 1.00$';
        else
            legend_txt{i} = sprintf('$w_0 = 10^{%d}$', num_parts(2));
        end
    else
        legend_txt{i} = sprintf('$w = %.2f \\times 10^{%d}$', num_parts(1), num_parts(2));
    end
end

params = set_plot_params(3, title_txt, legend_txt);

plot_error_lambda(data, color_profile, params)
% print(im_save_txt, '-depsc', '-opengl');

%{
figure;
color_profile = [0 0 1; 1 0 0; 0 1 0.5; 0 0 0; 0 1 1; 1 0 1; 1 1 0];
ax = axes;
for i = 1:n_w
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
    str = sprintf('%.2e', w_vals(i));
    num_parts = str2double(split(string(str), 'e'));
    if num_parts(1) == 1
        lgd_txt = sprintf('$w = 10^{%d}$', num_parts(2));
    else
        lgd_txt = sprintf('$w = %.2f \\times 10^{%d}$', num_parts(1), num_parts(2));
    end
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
axis([1 n_lambda 10^-3.5 10^2]);
%}
