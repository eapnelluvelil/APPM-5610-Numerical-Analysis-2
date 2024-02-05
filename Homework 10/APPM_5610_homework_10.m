clc;
clear;
close all;

PI = 2*asin(1);

m = 7;
n = 2^(m) - 1;

a = 0;
b = 1;
x = linspace(a, b, n)'; y = x;
[X, Y] = meshgrid(x, y);

f = @(x, y) -cos(PI*x).*sin(PI*y);

f_eval = f(X, Y);
f_eval_dst_yx = dst(dst(f_eval).').';

% Compute the sine coefficients of the solution
p = (1:n)'; q = p;
[P, Q] = meshgrid(p, q);
mat = P.^(2) + Q.^(2);
u_hat = f_eval_dst_yx./(-PI^(2)*mat);

% Apply inverse DST to recover solution at gridpoints
u = idst(idst(u_hat).').';

figure;
surf(X, Y, u);
colorbar;