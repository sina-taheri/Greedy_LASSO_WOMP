close all;
clear;
clc;

addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\data\Fig_1
addpath C:\Users\sinat\Dropbox\Shared_Sina_Simone\Numerics\WOMP_Sina_numerics\WOMP_paper_numerics\src\utils

WOMP_type = 'lassol1';
switch WOMP_type
    case 'lassol0'
        load('Fig_1_11_lassol0_N_300_s_10_m_150_n_iter_30');
        title_txt = '$\ell^0$-based LASSO-WOMP';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_1\\Fig_1_11_lassol0_N_300_s_10_m_150_n_iter_30';
    case 'srlassol0'
        load('Fig_1_12_srlassol0_N_300_s_10_m_150_n_iter_30');
        title_txt = '$\ell^0$-based SR-LASSO-WOMP';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_1\\Fig_1_12_srlassol0_N_300_s_10_m_150_n_iter_30';
    case 'lassol1'
        load('Fig_1_21_lassol1_N_300_s_10_m_150_n_iter_30');
        title_txt = '$\ell^1$-based LASSO-WOMP';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_1\\Fig_1_21_lassol1_N_300_s_10_m_150_n_iter_30';
    case 'srlassol1'
        load('Fig_1_22_srlassol1_N_300_s_10_m_150_n_iter_30');
        title_txt = '$\ell^1$-based SR-LASSO-WOMP';
        im_save_txt = 'C:\\Users\\sinat\\Dropbox\\Shared_Sina_Simone\\Numerics\\WOMP_Sina_numerics\\WOMP_paper_numerics\\figs\Fig_1\\Fig_1_22_srlassol1_N_300_s_10_m_150_n_iter_30';
end

%% Display - error vs. lambda
color_profile = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];

% data structure
data.Y = lambda_err;
data.X = lambda_vals;
data.nX = n_lambda;
data.PLOT = eta;
data.nPLOT = n_eta;

for i = 1:data.nPLOT
    str = sprintf('%.2e', eta(i));
    num_parts = str2double(split(string(str), 'e'));
    if num_parts(1) == 1
        legend_txt{i} = sprintf('$\\eta = 10^{%d}$', num_parts(2));
    else
        legend_txt{i} = sprintf('$\\eta = %.2f \\times 10^{%d}$', num_parts(1), num_parts(2));
    end
end

params = set_plot_params(1, title_txt, legend_txt);

plot_error_lambda(data, color_profile, params)

% print(im_save_txt, '-depsc', '-opengl');
% saveas(gcf, im_save_txt, 'epsc');   % saveas is problematic with figure size, and I'm saving the images by hand. I have to find a better alternative.