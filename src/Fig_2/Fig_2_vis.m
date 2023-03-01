close all;
clear;
clc;

addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\data\Fig_2
addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\src\utils

WOMP_type = 'lassol0';
switch WOMP_type
    case 'lassol0'
        load('Fig_2_11_lassol0_N_426_d_5_m_200_n_iter_50');
        title_txt = '$\ell^0$-based LASSO-WOMP';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_2\\Fig_2_11_lassol0_N_426_d_5_m_200_n_iter_50';
    case 'srlassol0'
        load('Fig_2_12_srlassol0_N_426_d_5_m_200_n_iter_50');
        title_txt = '$\ell^0$-based SR-LASSO-WOMP';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_2\\Fig_2_12_srlassol0_N_426_d_5_m_200_n_iter_50';
    case 'lassol1'
        load('Fig_2_21_lassol1_N_426_s_5_m_200_n_iter_50');
        title_txt = '$\ell^1$-based LASSO-WOMP';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_2\\Fig_2_21_lassol1_N_426_d_5_m_200_n_iter_50';
    case 'srlassol1'
        load('Fig_2_22_srlassol1_N_426_d_5_m_200_n_iter_50');
        title_txt = '$\ell^1$-based SR-LASSO-WOMP';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_2\\Fig_2_22_srlassol1_N_426_d_5_m_200_n_iter_50';
end

%% Display - error vs. lambda
color_profile = [0 0 1; 1 0 0; 0 1 0.5; 0 0 0; 0 1 1; 1 0 1; 1 1 0];

% data structure
data.Y = lambda_err;
data.X = lambda_vals;
data.nX = n_lambda;
data.PLOT = eta;
data.nPLOT = n_eta;

% params structure
for i = 1:data.nPLOT
    str = sprintf('%.2e', eta(i));
    num_parts = str2double(split(string(str), 'e'));
    if num_parts(1) == 1
        legend_txt{i} = sprintf('$\\eta = 10^{%d}$', num_parts(2));
    elseif num_parts(1) == 0
        legend_txt{i} = '$\eta = 0$';
    else
        legend_txt{i} = sprintf('$\\eta = %.2f \\times 10^{%d}$', num_parts(1), num_parts(2));
    end
end

params = set_plot_params(2, title_txt, legend_txt);

plot_error_lambda(data, color_profile, params)

% print(im_save_txt, '-depsc', '-opengl');

% saveas(gcf, im_save_txt, 'epsc');
%{
figure;
color_profile = [0 0 1; 1 0 0; 0 1 0.5; 0 0 0; 0 1 1; 1 0 1; 1 1 0];
% ax = axes('Position', [0.08 0.093 0.90 0.85]);
ax = axes;
for i = 1:n_eta
%     lambda_err1 = reshape(lambda_err(:, i, :), n_trial, n_lambda);
%     figure;
%     boxplot(lambda_err1, lambda_vals, 'color', color_profile(i, :));
    hbp = boxplot(lambda_err{i}, lambda_vals, 'color', color_profile(i, :));
    h = findobj(hbp, 'tag', 'Outliers');
    for iH = 1:length(h)
        h(iH).MarkerEdgeColor = color_profile(i, :);
    end
    % this block of code removes half of the boxplots to make the final
    % plot less messy
    %{
    h = findobj(hbp);
    ih = [];
    il = 1:length(h);
    for iH = 0:2:n_lambda - 1
        ih = [ih, il(7*iH+1:7*iH + 7)];
    end
    delete(h(ih));
    %}
    hAx = gca;
    xtk = hAx.XTick; 
    hold on;
    str = sprintf('%.2e', eta(i));
    num_parts = str2double(split(string(str), 'e'));
    if num_parts(1) == 1
        lgd_txt = sprintf('$\\eta = 10^{%d}$', num_parts(2));
    elseif num_parts(1) == 0
        lgd_txt = '$\eta = 0$';
    else
        lgd_txt = sprintf('$\\eta = %.2f \\times 10^{%d}$', num_parts(1), num_parts(2));
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
ylabel('$ Relative-L_\varrho^2(\mathcal{U})\ error $', 'Fontsize', 20, 'FontWeight', 'bold', 'Interpreter', 'latex');
lgd = legend('show', 'Location', 'northwest');
lgd.FontSize = 20;
lgd.FontName = 'Times';
lgd.Interpreter = 'latex';
grid on;
grid minor;
title(title_txt, 'FontName', 'Times', 'Fontsize', 20, 'Interpreter', 'latex');
axis([1 n_lambda 10^-4.5 10^0.5]);
%}