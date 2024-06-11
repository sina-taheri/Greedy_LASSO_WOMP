close all;
clear;
clc;

addpath '../src/utils'

%% Experiment hyperparameters
m = 100;
N = 200;
s = 15;
n_iter = 150;

% Experiment-specific macros
experiment = 'srlasso';
switch experiment
    case 'lasso'
        womp_type = 'lassol1';
        lambda_omp_vals = [0, 10.^[-8, -6:1]];
        lambda_cvx_vals = 1.78e-4;
%         lambda_cvx_vals = [0, 10.^[-6:0.5:-1]];
        prg_cvx = 'wlasso';
        Lev_c = 0;
        save_txt = sprintf('C:...\\data\\Fig_4\\Fig_4_11_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', womp_type, N, s, m, n_iter);
    case 'srlasso'
        womp_type = 'srlassol1';
%         lambda_omp_vals = 1.58e-1;
        lambda_omp_vals = [0, 0.158*(10.^(-0.25:0.25:0.25))];
        lambda_cvx_vals = 1.26e-1;
%         lambda_cvx_vals = [0, 10.^[-1.5:0.2:-0.5]];
        prg_cvx = 'wsrlasso';
        Lev_c = 0;
        save_txt = sprintf('C:...\\data\\Fig_4\\Fig_4_12_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', womp_type, N, s, m, n_iter);
    case 'ladlasso'
        womp_type = 'ladlassol1';
        lambda_cvx_vals = 1;
        prg_cvx = 'wladlasso';
        Lev_c = 0.05;
        save_txt = sprintf('C:...\\data\\Fig_4\\Fig_4_13_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', womp_type, N, s, m, n_iter);
end

eta = 1e-3;

K = round(Lev_c*m);
Mc = 100;

n_lambda_omp = length(lambda_omp_vals);
n_lambda_cvx = length(lambda_cvx_vals);

w = ones(N, 1);
n_plot = n_lambda_omp + n_lambda_cvx;

%% Main routine
n_trial = 30;
lambda_err = zeros(n_trial, n_plot, n_iter + 1);

%% Running the experiment
for i = 1:n_trial
    disp(i)
    % Building the compressed sensing model
    [x, A, noise, corruption] = cs_synth_model(m, N, s, eta, K, Mc);
    y_n = A*x + corruption + noise;
    
    err = zeros(n_plot, n_iter + 1);
    
    % Signal recovery via WOMP
    for j = 1:n_lambda_omp
        lambda = lambda_omp_vals(j);
        x_vals = womp(A, y_n, w, lambda, n_iter, womp_type);
        err(j, :) = vecnorm(repmat(x, [1, n_iter + 1]) - x_vals)'/norm(x);
    end
    
    % Signal recovery via CVX
    for k = 1:n_lambda_cvx
        lambda_cvx = lambda_cvx_vals(k);
        x_vals = womp(A, y_n, w, lambda_cvx, n_iter, prg_cvx);
        err(n_lambda_omp + k, :) = vecnorm(repmat(x, [1, n_iter + 1]) - x_vals)'/norm(x);
    end
    
    lambda_err(i, :, :) = err;
end

%% Save the results
lambda_vals = [lambda_omp_vals, lambda_cvx_vals];
% data = struct('lambda_err', lambda_err, 'lambda_vals', lambda_vals, 'n_lambda_omp', n_lambda_omp, ...
%     'n_lambda_cvx', n_lambda_cvx, 'n_plot', n_plot, 'eta', eta, 'K', K, 'womp_type', womp_type, 'n_iter', n_iter);
% save(save_txt, 'data');
save(save_txt, 'lambda_err', 'lambda_vals', 'n_lambda_omp', 'n_lambda_cvx', 'n_plot', 'eta', 'K', 'womp_type', 'n_iter');

%% Display - error vs. iterations
std_err = std(log10(lambda_err), 0, 1);
std_err = reshape(std_err, n_plot, n_iter + 1);
mean_err = mean(log10(lambda_err), 1);
mean_err = reshape(mean_err, n_plot, n_iter + 1);
upper_err = mean_err + std_err;
lower_err = mean_err - std_err;

colors = color_profile(n_plot);
% colors = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0; 0.4 0.4 0.4; 0.8 0.4 0.8];
ax = axes('Position', [0.096 0.093 0.88 0.89]);
ax.ActivePositionProperty = 'outerposition';
for i = 1:n_plot
    fig = fill(ax, [(1:n_iter + 1)'; (n_iter + 1:-1:1)'], [10.^lower_err(i, :), 10.^rot90(upper_err(i, :), 2)],...
        colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    set(get(get(fig,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on;
    string = sprintf('%.2e', lambda_vals(i)); % Convert number to a scientific notation string with 2 decimal places of precision
    stringParts = strsplit(string, 'e'); % Split the string where 'e' is
    a = str2double(stringParts(1));
    b = str2double(stringParts(2));
    if a == 1
        lgd_txt = sprintf('$\\lambda = 10^{%d}$', b);
    elseif b == 0
        lgd_txt = sprintf('$\\lambda = %.2f$', a);
    else
        lgd_txt = sprintf('$\\lambda = %.2f \\times 10^{%d}$', a, b);
    end
%     lgd_txt = sprintf('%2.0e', lambda_vals(i));
    plot(ax, 10.^mean_err(i, :), 'LineWidth', 1.7, 'color', colors(i, :), 'DisplayName', lgd_txt);
end
set(gca, 'YScale', 'log', 'FontSize', 12);
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.7 0.95]);
xlabel('Iteration Number', 'FontName', 'Times', 'Fontsize', 18, 'FontWeight', 'bold'); xlim([1 n_iter + 1]);
ylabel('$||x - x^\sharp||_2/||x||_2$', 'Interpreter', 'Latex', 'Fontsize', 18, 'FontWeight', 'bold'); % ylim([0.03 3]);
% title('$\ell_1$-based LAD-LASSO WOMP', 'Interpreter', 'Latex');
lgd = legend('show', 'Location', 'northwest');
lgd.FontSize = 16;
lgd.FontName = 'Times';
lgd.Interpreter = 'latex';
grid on;
grid minor;

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
