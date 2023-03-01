% close all;
clear;
clc;

addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\src\utils
addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\data

%% Function approximation setup
d = 5;               % dimension of the domain
func = @iso_exp;     % isotropic exponential function (See book, f1 in Appendix A)
index_type = 'HC';   % use hyperbolic cross index set
err_grid_ratio = 10; % ratio between maximum index set size and error grid size
N_max = 500;        % Maximum size of ambient space
poly_setting = 2;    % 1 = Leg + Unif, 2 = Cheby + Cheby

switch poly_setting
    case 1
        poly_type = 'legendre'; % use Legendre polynomials
        samp_type = 'uniform';  % use samples drawn randomly from the uniform distribution
    case 2
        poly_type = 'chebyshev'; % use Chebyshev polynomials
        samp_type = 'chebyshev'; % use samples drawn randomly from the Chebyshev distribution
    otherwise
        error('invalid column number');
end

%%% Construct index set and define N %%%

n = find_order(index_type,d,N_max); % find the maximum polynomial order for the given index_type
I = generate_index_set(index_type,d,n); % compute index set
N = size(I,2); % number of basis functions

u = generate_intrinsic_weights(poly_type,I); % generate the intrinsic weights

%%% Construct error grid, error matrix and error vector %%%

M = ceil(err_grid_ratio*N);
err_grid = generate_sampling_grid(samp_type,d,M);
A_err_grid = generate_measurement_matrix(poly_type,I,err_grid);
b_err_grid = func(err_grid)/sqrt(M);

%%% Compute a reference solution via least squares on the error grid %%%

c_ref = A_err_grid \ b_err_grid; %% reference solution (used to compute the error)

%% Experiments hyperparameters
m = 200;

womp_type = 'srlassol1';
switch womp_type
    case 'lassol0'
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_2\\Fig_2_Fig_11_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', womp_type, N, d, m, n_iter);
    case 'srlassol0'
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_2\\Fig_2_Fig_12_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', womp_type, N, d, m, n_iter);
    case 'lassol1'
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_2\\Fig_2_Fig_21_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', womp_type, N, d, m, n_iter);
    case 'srlassol1'
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_2\\Fig_2_Fig_22_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', womp_type, N, d, m, n_iter);
end

n_iter = 30;

lambda_vals = [10.^(-8:0.25:0)];
n_lambda = numel(lambda_vals);

eta = [10.^(-3:-1), 0];
n_eta = numel(eta);

%% Main routine
n_trial = 30;
lambda_err = cell(1, n_eta);

for i = 1:n_trial
    %% Generating the compressed sensing model for function approximation
    y_grid = generate_sampling_grid(samp_type,d,m); % generate sample points
    A = generate_measurement_matrix(poly_type,I,y_grid); % generate measurement matrix
    b = func(y_grid)/sqrt(m); % generate measurement vector

    %A = normc(randn(m, N));
    %ix = randi(N, 1, s);
    %x = zeros(N, 1);
    %x(ix) = randn(s, 1);
    %x = x / norm(x); %% Normalization of the signal 
    %y = A*x;

    %% Adding noise to the signal
    %     err = zeros(n_snr, n_lambda);
    for j = 1:n_eta

        %% Add normally distributed and rescaled noise of l2-norm eta
        rand_vec = randn(m,1);
        noise    = eta(j) * rand_vec / norm(rand_vec);
        b_noisy  = b + noise;

        %% Signal recovery
        for k = 1:n_lambda
            lambda   = lambda_vals(k);
            norms    = sqrt(sum(abs(A).^2,1));
            A_norm   = A * diag(1 ./ norms);
            x_vals_1 = womp(A_norm, b_noisy, u, lambda, n_iter, womp_type);
            x_vals   = diag(1 ./ norms) * x_vals_1;
            c        = x_vals(:, n_iter + 1); %% computed vector of coefficients
            lambda_err{j}(i, k) = norm(A_err_grid * c - b_err_grid)/norm(b_err_grid); % compute L^2_rho-norm error
            % lambda_err{j}(i, k) = norm(A_err_grid*c - b_err_grid,Inf)/norm(b_err_grid,Inf); % compute L^inf-norm error
        end
    end
    %     lambda_err(i, : , :) = err;
end

txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_2\\Fig_2_21_%s_N_%d_d_%d_m_%d_n_iter_%d.mat', womp_type, N, d, m, n_iter);
save(txt, 'lambda_err', 'lambda_vals', 'n_lambda', 'eta', 'n_eta', 'womp_type');

% clear all;
% load('Fig_11_500_10_150_30.mat');

%% Display - error vs. lambda
figure;
color_profile = [0 0 1; 1 0 0; 0 1 0.5; 0 0 0; 0 1 1; 1 0 1; 1 1 0];
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
% txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\\Fig_2\\Fig_12_%s_N_%d_s_%d_m_%d_n_iter_%d.eps', womp_type, N, d, m, n_iter);
% saveas(gcf, txt);
