function x = lad_solver(A, y, type)
[m, N] = size(A);
y_n = y;


%% Solving LAD

if N==1
    % Direct solution in the one dimensional case (N = 1)
    a = A;
    enum = 1:m;
    ind_nz = enum(abs(a) > 1e-10);

    [t, ind] = sort((y_n(ind_nz)./a(ind_nz)), 'ascend');
    a1 = a(ind);
    K = find(cumsum(abs(a1)/norm(a1, 1)) >= 1/2, 1);
    p = ind_nz(ind(K)); % The true index, i.e x = y(p)/A(p)
    x = t(K);
else
    switch type
        case 'linprog'
            % Solution via Linear Programming
            f = [zeros(N, 1); zeros(N, 1); ones(m, 1)];
            B = [-A, A, -eye(m); A, -A, -eye(m)];
            b = [-y_n; y_n];
            lb = zeros(2*N + m, 1);
            ub = inf(2*N + m, 1);

            options = optimoptions('linprog', 'Display', 'off');
            x_nz = linprog(f, B, b, [], [], lb, ub, options);
            x = x_nz(1:N, :) - x_nz(N + 1:2*N, :);
            
        case 'cvx'
        % Solution via CVX (Solves an LP)
        cvx_begin quiet
            variable x(N)
            minimize(norm(A*x - y_n, 1))
        cvx_end
        
        otherwise
            error('LAD solver type not supported');
    end
    
end

%{
close all;
clear;
clc;

%% Generating the model
m = 512;

%{
a = randn(m, 1);
a = a/norm(a, 1);

x0 = randn(1);

y = x0*a;
%}
N = 10;

x0 = randn(N, 1);
A = randn(m, N);
y = A*x0;

%% Adding noise to the signal
eta = 10^-4;
e = randn(m, 1);
noise = eta*norm(y)*e/norm(e);

%% Adding corruption
k = 50;
ik = randi(m, 1, k);
c = zeros(m, 1);
c(ik) = randn(k, 1);

y_n = y + noise + c;
%}

