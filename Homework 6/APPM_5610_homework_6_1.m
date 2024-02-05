% Problem 1
clc;
clear;

% Initial conditions
y_0 = [1/2; 0];
a = 1e-15;
b = 3*pi;

f = @(t) [-t/8; ...
          (-3/8) - (1/2) + t.^(2)/16];

tol = 1e-11;

max_rows = 10;
h = (b - a);

R = zeros(max_rows + 1);

for j = (max_rows + 1)
    y = trap_rule(f, a, b, h/2^(j-1), y_0);
    R(j, 1) = y(1);
end

for k = 2:(max_rows + 1)
    for j = k:(max_rows + 1)
        R(j, k) = (4^(k-1)*R(j, k-1) - R(j-1, k-1))/(4^(k-1) - 1);
    end
end

function [y] = trap_rule(f, a, b, h, y_0)
    t_0 = a;

    while t_0 <= b
       y_n = y_0 + (h/2)*(f(t_0) + f(t_0 + h));
       y_0 = y_n;
       t_0 = t_0 + h;
       y = y_n;
    end
end