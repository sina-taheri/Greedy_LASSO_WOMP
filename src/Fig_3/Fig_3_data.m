clear
clc
% close all

addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\src\utils
addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\data


%% Experiment hyperparameters
N = 500;
s = 10;

% oracle_type = 'exact';   % knows s out of s indices
oracle_type = 'smaller'; % knows s/2 out of s indices
%oracle_type = 'larger';  % knows 2s indices contains the s real ones
%oracle_type = 'hybrid';  % knows s/2 out of s indices, + s/5 wrong indices

experiment = 'lassol1_low';
switch experiment
    case 'lassol1_low'
        womp_type = 'lassol1';
        m = 40;
        corruption_lev = 0;
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_3\\Fig_3_11_smaller_%s_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', oracle_type, womp_type, N, s, m, n_iter);
    case 'lassol1_high'
        womp_type = 'lassol1';
        m = 100;
        corruption_lev = 0;
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_3\\Fig_3_21_smaller_%s_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', oracle_type, womp_type, N, s, m, n_iter);
    case 'srlassol1_low'
        womp_type = 'srlassol1_revised';
        m = 35;
        corruption_lev = 0;
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_3\\Fig_3_12_smaller_%s_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', oracle_type, womp_type, N, s, m, n_iter);
    case 'srlassol1_high'
        womp_type = 'srlassol1_revised';
        m = 100;
        corruption_lev = 0;
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_3\\Fig_3_22_smaller_%s_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', oracle_type, womp_type, N, s, m, n_iter);
    case 'ladlassol1_low'
        womp_type = 'ladlassol1';
        m = 80;
        corruption_lev = 0.1;
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_3\\Fig_3_13_smaller_%s_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', oracle_type, womp_type, N, s, m, n_iter);
    case 'ladlassol1_high'
        womp_type = 'ladlassol1';
        m = 100;
        corruption_lev = 0.1;
        save_txt = sprintf('C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\data\\Fig_3\\Fig_3_23_smaller_%s_%s_N_%d_s_%d_m_%d_n_iter_%d.mat', oracle_type, womp_type, N, s, m, n_iter);
end

n_corruption = round(corruption_lev*m);
Mc = 100;
eta = 1e-3;

n_iter = 3*s;

lambda_vals = 10.^(-6:0.2:2);
n_lambda = length(lambda_vals);

w_vals = 10.^[0, -3]; % weights value in support
n_w = numel(w_vals);

%% Main routine
n_trial = 30;

%% Main routine
for i = 1:n_trial
    disp(i);
    % random Gaussian matrix
    A = randn(m,N)/sqrt(m);
    %imagesc(abs(A))

    % generat sparse random vector
    rand_perm = randperm(N);
    support = rand_perm(1:s);
    x = zeros(N,1);
    x(support) = randn(s,1);
    %figure
    %stem(x)
    corruption = zeros(m, 1);
    corruption(randperm(m, n_corruption)) = Mc*randn(n_corruption, 1);
    
    b = A*x + eta * randn(m,1) + corruption;

    for i_lambda = 1:n_lambda

        for i_w = 1:length(w_vals)

            w = ones(N,1);

            switch oracle_type
                case 'exact'
                    oracle_support = support;
                case 'smaller'
                    oracle_support = rand_perm(1:round(s/2)); % half of the support known
                case 'larger'
                    oracle_support = rand_perm(1:2*s); % twice the support known
                case 'hybrid'
                    oracle_support = rand_perm([1:round(s/2), s+1:s+round(s/5)]); % knowns s/2 good indices and s/5 bad indices
            end
            w(oracle_support) = w_vals(i_w); % set weight on known support

            lambda   = lambda_vals(i_lambda);
            norms    = sqrt(sum(abs(A).^2,1));
            A_norm   = A * diag(1 ./ norms);
            x_vals_1 = womp(A_norm, b, w, lambda, n_iter, womp_type);
            x_vals   = diag(1 ./ norms) * x_vals_1;
            c        = x_vals(:, n_iter + 1); %% computed vector of coefficients
            lambda_err{i_w}(i, i_lambda) = norm(x - c)/norm(x); % compute L^2_rho-norm error
            % lambda_err{j}(i, k) = norm(A_err_grid*c - b_err_grid,Inf)/norm(b_err_grid,Inf); % compute L^inf-norm error
        end
    end
end

save(save_txt, 'lambda_err', 'lambda_vals', 'n_lambda', 'w_vals', 'n_w', 'womp_type');

%% Display - error vs. lambda
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
    lgd_txt = sprintf('$w = %.2f \\times 10^{%d}$', num_parts(1), num_parts(2));
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