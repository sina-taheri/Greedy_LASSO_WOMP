close all;
clear;
clc;

m = 20;
n = 40;
w = 2;

A = zeros(m, w);
% A(9:10, :) = zeros(2, n);
% A(33:34, :) = zeros(2, n);
% A(21:22, :) = zeros(2, n);
% A(:, 3:4) = zeros(m, 2);

f = imagesc(A);

for i = 0:m/w - 1
    line(linspace(0, n + 1, 100)', (i*w + 0.5)*ones(100, 1), 'LineWidth', 2, 'color', 'white');
end
for i = 0:n/w - 1
    line((i*w + 0.5)*ones(100, 1), linspace(0, m + 1, 100)' , 'LineWidth', 2, 'color', 'white');
end
        