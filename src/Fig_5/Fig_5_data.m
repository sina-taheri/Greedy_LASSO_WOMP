close all;
clear;
clc;

addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\src\utils

%% Function approximation initialization
d = 10;               % dimension of the domain
func = @iso_exp;     % isotropic exponential function (See book, f1 in Appendix A)
index_type = 'HC';   % use hyperbolic cross index set
err_grid_ratio = 10; % ratio between maximum index set size and error grid size
N_max = 500;        % Maximum size of ambient space

poly_setting = 1;    % 1 = Leg + Unif, 2 = Cheby + Cheby
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

% Construct index set and define N %%%
n = find_order(index_type,d,N_max); % find the maximum polynomial order for the given index_type
I = generate_index_set(index_type,d,n); % compute index set
N = size(I,2); % number of basis functions

u = generate_intrinsic_weights(poly_type,I); % generate the intrinsic weights

% Construct error grid, error matrix and error vector %%%
Mc = ceil(err_grid_ratio*N);
err_grid = generate_sampling_grid(samp_type,d,Mc);
A_err_grid = generate_measurement_matrix(poly_type,I,err_grid);
b_err_grid = func(err_grid)/sqrt(Mc);

% Compute a reference solution via least squares on the error grid %%%
c_ref = A_err_grid \ b_err_grid; %% reference solution (used to compute the error)

%% Experiment hyperparameters
m = 200;
% lambda_omp_vals = [];
n_iter = 250;
experiment = 'srlasso';
switch experiment
    case 'lasso'
        womp_type = 'lassol1';
        lambda_omp_vals = [0, (3.16e-5)*(10.^(-3:3:3))];
%         lambda_omp_vals = [0, 10.^[-8, -6:0.5:-4, -3, -2, -1, 0]];
        lambda_cvx_vals = 3.16e-5;
%         lambda_cvx_vals = [0, 10.^[-6:0.5:-1]];
        prg_cvx = 'wlasso';
        lev_c = 0;
    case 'srlasso'
        womp_type = 'srlassol1_revised';
        lambda_omp_vals = [0, (1.78e-2)*(10.^(-0.5:0.5:0.5))];
        lambda_cvx_vals = 1.78e-2;
%         lambda_cvx_vals = [0, 10.^[-6:0.25:0]];
        prg_cvx = 'wsrlasso';
        lev_c = 0;
    case 'ladlasso'
        womp_type = 'ladlassol1';
        lambda_omp_vals = [0, 10.^[-8, -6:0.5:-4, -3, -2, -1, 0]];
        lambda_cvx_vals = 1;
        prg_cvx = 'wladlasso';
        lev_c = 0.05;
end

n_lambda_omp = length(lambda_omp_vals);
n_lambda_cvx = length(lambda_cvx_vals);
eta = 1e-3;
Kc = round(lev_c*m);
Mc = 50;
n_plot = n_lambda_omp + n_lambda_cvx;

%% Main routine
n_trial = 30;
lambda_err = zeros(n_trial, n_plot, n_iter + 1);

for i = 1:n_trial
    disp(i)
    % Sampling the function
    y_grid = generate_sampling_grid(samp_type, d, m); % generate sample points
    A = generate_measurement_matrix(poly_type, I, y_grid); % generate measurement matrix
    norms    = sqrt(sum(abs(A).^2,1));
    A_norm   = A * diag(1 ./ norms);
    b = func(y_grid)/sqrt(m); % generate measurement vector
    
    % Add normally distributed and rescaled noise of l2-norm eta, and or high
    % magnitude sparse corruptions
    rand_vec = randn(m,1);
    noise    = eta * rand_vec / norm(rand_vec);
    
    ik = randi(m, 1, Kc);
    corruption = zeros(m, 1);
    corruption(ik) = Mc*randn(Kc, 1);
    
    b_noisy  = b + noise + corruption;
    
    err = zeros(n_plot, n_iter + 1);
    
    % Signal recovery via WOMP
    for j = 1:n_lambda_omp
        lambda   = lambda_omp_vals(j);
        x_vals_1 = womp(A_norm, b_noisy, u, lambda, n_iter, womp_type);
        x_vals   = diag(1 ./ norms) * x_vals_1;
        err(j, :) = vecnorm(A_err_grid * x_vals - repmat(b_err_grid, [1, n_iter + 1]))/norm(b_err_grid); % compute L^2_rho-norm error
    end
    
    % Signal recovery via CVX
    for k = 1:n_lambda_cvx
        lambda_cvx = lambda_cvx_vals(k);
        x_vals_1 = womp(A_norm, b_noisy, u, lambda_cvx, n_iter, prg_cvx);
        x_vals   = diag(1 ./ norms) * x_vals_1;
        err(n_lambda_omp + k, :) = vecnorm(A_err_grid * x_vals - repmat(b_err_grid, [1, n_iter + 1]))/norm(b_err_grid); % compute L^2_rho-norm error
    end
    
    lambda_err(i, :, :) = err;
end

%% Save the results
lambda_vals = [lambda_omp_vals, lambda_cvx_vals];
switch experiment
    case 'lasso'
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_5\\Fig_5_13_%s_d_%d_N_%d_m_%d_n_iter_%d.mat', womp_type, d, N, m, n_iter);
    case 'srlasso'
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_5\\Fig_5_23_%s_d_%d_N_%d_m_%d_n_iter_%d.mat', womp_type, d, N, m, n_iter);
    case 'ladlasso'
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_5\\Fig_5_33_%s_d_%d_N_%d_m_%d_n_iter_%d.mat', womp_type, d, N, m, n_iter);
end
save(save_txt, 'lambda_err', 'lambda_vals', 'n_lambda_omp', 'n_lambda_cvx', 'n_plot', 'eta', 'Kc', 'womp_type', 'n_iter');

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
