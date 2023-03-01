close all;
clear;
clc;

addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\data\Fig_4
addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\src\utils

womp_type = 'srlassol1';
switch womp_type
    case 'lassol1'
        load('Fig_4_11_lassol1_N_200_s_15_m_100_n_iter_150');
        title_txt = '$\ell^1$-based LASSO-WOMP';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_4\\Fig_4_11_lassol1_N_200_s_15_m_100_n_iter_150';
        % extract a portion of lambdas (if required)
        [~, idx] = intersect(lambda_vals, [0, 10^-1, 10^-3, 10^-5], 'stable');
        lambda_vals = lambda_vals([idx', end]);
        lambda_err = lambda_err(:, [idx', end], :);
        n_plot = numel(lambda_vals);
    case 'srlassol1'
        load('Fig_4_12_srlassol1_N_200_s_15_m_100_n_iter_150');
        title_txt = '$\ell^1$-based SR-LASSO-WOMP';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_4\\Fig_4_12_srlassol1_N_200_s_15_m_100_n_iter_150';
        % extract a portion of lambdas (if required)
%         [~, idx] = intersect(lambda_vals, [0, 10^-1, 1, 10], 'stable');
%         lambda_vals = lambda_vals([idx', end]);
%         lambda_err = lambda_err(:, [idx', end], :);
%         n_plot = numel(lambda_vals);
    case 'ladlassol1'
        load('Fig_4_13_ladlassol1_N_200_s_15_m_100_n_iter_150');
        title_txt = '$\ell^1$-based LAD-LASSO-WOMP';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_4\\Fig_4_13_ladlassol1_N_200_s_15_m_100_n_iter_150';
        % extract a portion of lambdas (if required)
        [~, idx] = intersect(lambda_vals, [10, 1, 0, 10^-1], 'stable');
        lambda_vals = lambda_vals([idx', end]);
        lambda_err = lambda_err(:, [idx', end], :);
        n_plot = numel(lambda_vals);
end

%% Display - error vs. iterations
colors = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 0 0 0; 1 0 1; 1 1 0];

% data structure
data.Y = lambda_err;
data.nX = n_iter;
data.PLOT = lambda_vals;
data.nPLOT = n_plot;

% params structure
params = set_plot_params(4, title_txt);

plot_error_iteration(data, colors, params)
% print(im_save_txt, '-depsc', '-opengl');

%{
std_err = std(log10(lambda_err), 0, 1);
std_err = reshape(std_err, n_plot, n_iter + 1);
mean_err = mean(log10(lambda_err), 1);
mean_err = reshape(mean_err, n_plot, n_iter + 1);
upper_err = mean_err + std_err;
lower_err = mean_err - std_err;

% colors = color_profile(n_plot);
% colors = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];
figure;
ax = axes;
% ax.ActivePositionProperty = 'outerposition';
for i = 1:n_plot
    fig = fill(ax, [(1:n_iter + 1)'; (n_iter + 1:-1:1)'], [10.^lower_err(i, :), 10.^rot90(upper_err(i, :), 2)],...
        colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    set(get(get(fig,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on;
    string = sprintf('%.2e', lambda_vals(i)); % Convert number to a scientific notation string with 2 decimal places of precision
    stringParts = strsplit(string, 'e'); % Split the string where 'e' is
    a = str2double(stringParts(1));
    b = str2double(stringParts(2));
    if a == 1 && b == 0
        lgd_txt = '$\lambda = 1$';
    elseif a == 1 && b == 1
        lgd_txt = '$\lambda = 10$';
    elseif b == 0
        lgd_txt = sprintf('$\\lambda = %.2f$', a);
    elseif a == 1
        lgd_txt = sprintf('$\\lambda = 10^{%d}$', b);
    else
        lgd_txt = sprintf('$\\lambda = %.2f \\times 10^{%d}$', a, b);
    end
%     lgd_txt = sprintf('%2.0e', lambda_vals(i));
    plot(ax, 10.^mean_err(i, :), 'LineWidth', 1.7, 'color', colors(i, :), 'DisplayName', lgd_txt);
end
set(gca, 'YScale', 'log', 'FontSize', 16);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 1]);
xlabel('Iteration Number', 'FontName', 'Times', 'Fontsize', 20); xlim([1 n_iter + 1]);
ylabel('$||x - x^\sharp||_2/||x||_2$', 'Interpreter', 'Latex', 'Fontsize', 20, 'FontWeight', 'bold'); % ylim([0.03 3]);
lgd = legend('show', 'Location', 'northwest');
lgd.FontSize = 20;
lgd.FontName = 'Times';
lgd.Interpreter = 'latex';
grid on;
grid minor;
title(title_txt, 'FontName', 'Times', 'Fontsize', 20, 'Interpreter', 'latex');
axis([1 n_iter 10^-4.5 10]);

% This function creates a color pool that are discriminative enough for
% plots.
function colors = color_profile(n)
    p = n^(1/3);
    if mod(p, 1) == 0
        n = n + 1;
        p = n^(1/3);
    end
    lp = floor(p);
    up = ceil(p);
    
    a = up;
    b = lp;
    c = lp;
    if a*b*c < n
        b = up;
        if a*b*c < n
            c = up;
        end
    end
        
    colors = [];
    for i = 1:a
        for j = 1:b
            for k = 1:c
                colors = [colors; (i - 1)/a, (j - 1)/b, (k - 1)/c];
            end
        end
    end
    colors = colors(1:n, :);
end
%}