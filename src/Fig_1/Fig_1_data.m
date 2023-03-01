close all;
clear;
clc;

addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\src\utils
addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\data

%% Experiment hyperparameters
m = 150;
N = 300;
s = 10;

womp_type = 'srlassol1';
switch womp_type
    case 'lassol0'
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_1\\Fig_1_11_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', womp_type, N, s, m, n_iter);
    case 'srlassol0'
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_1\\Fig_1_12_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', womp_type, N, s, m, n_iter);
    case 'lassol1'
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_1\\Fig_1_21_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', womp_type, N, s, m, n_iter);
    case 'srlassol1'
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_1\\Fig_1_22_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', womp_type, N, s, m, n_iter);
end

n_iter = 3*s;

lambda_vals = [0, 10.^(-8:0.25:2)];
n_lambda = numel(lambda_vals);

eta = 10.^(-3:-1);
n_eta = numel(eta);

w = ones(N, 1);

%% Main routine
n_trial = 30;
lambda_err = cell(1, n_eta);

for i = 1:n_trial    
    for i_eta = 1:n_eta
        [x, A, noise, corruption] = cs_synth_model(m, N, s, eta(i_eta));
        y_n = A*x + noise;
%         e = randn(m, 1);
%         y_n = y + eta(i_eta)*e/norm(e);
        
        %% Signal recovery
        for k = 1:n_lambda
            lambda = lambda_vals(k);
            x_vals = womp(A, y_n, w, lambda, n_iter, womp_type);
            lambda_err{i_eta}(i, k) = norm(x_vals(:, n_iter + 1) - x)/norm(x);
        end
    end
end

save(save_txt, 'lambda_err', 'lambda_vals', 'n_lambda', 'eta', 'n_eta', 'womp_type');

% clear all;
% load('Fig_11_500_10_150_30.mat');

%% Display - error vs. lambda
figure;
color_profile = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];
% ax = axes('Position', [0.08 0.093 0.90 0.85]);
ax = axes;
for i = 1:n_eta
%     lambda_err1 = reshape(lambda_err(:, i, :), n_trial, n_lambda);
%     figure;
%     boxplot(lambda_err1, lambda_vals, 'color', color_profile(i, :));
    hbp = boxplot(lambda_err{i}(:, 1:end), lambda_vals(1:end), 'color', color_profile(i, :));
    h = findobj(hbp, 'tag', 'Outliers');
    for iH = 1:length(h)
        h(iH).MarkerEdgeColor = color_profile(i, :);
    end
    hAx = gca;
    xtk = hAx.XTick; 
    hold on;
    str = sprintf('%.2e', eta(i));
    num_parts = str2double(split(string(str), 'e'));
    lgd_txt = sprintf('$\\eta = %.2f \\times 10^{%d}$', num_parts(1), num_parts(2));
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
lgd.FontSize = 16;
lgd.FontName = 'Times';
lgd.Interpreter = 'latex';
grid on;
grid minor;
title(sprintf('%s WOMP', womp_type), 'FontName', 'Times', 'Fontsize', 20, 'Interpreter', 'latex');
axis([1 n_lambda 1e-5 10]);
% txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\\Fig_1\\Fig_12_%s_N_%d_s_%d_m_%d_n_iter_%d.eps', WOMP_type, N, s, m, n_iter);
% saveas(gcf, txt);
