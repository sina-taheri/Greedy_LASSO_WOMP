function [x, A, noise, c] = cs_synth_model(m, N, s, varargin)
ix = randperm(N, s);
x = zeros(N, 1);
x(ix) = randn(s, 1);

A = randn(m, N);
A = A./repmat(sqrt(sum(A.^2)), [m, 1]);
% A = A/sqrt(m);

y = A*x;

%% Adding noise/corruption to the signal
switch nargin
    case 3
        eta = 0;
        K = 0;
        M = 1;
    case 4
        eta = varargin{1};
        K = 0;
        M = 1;
    case 5
        eta = varargin{1};
        K = varargin{2};
        M = 1;
    case 6
        eta = varargin{1};
        K = varargin{2};
        M = varargin{3};
    otherwise
        error('Too much input arguments!');
end

e = randn(m, 1);
noise = eta*e/norm(e);

ik = randi(m, 1, K);
c = zeros(m, 1);
c(ik) = M*randn(K, 1);

