clc;
clear;

tspan = [1e-15 3*pi];
y0    = [1/2; 0];

[t, y] = ode45(@bessel, tspan, y0);

% plot(t, t.*y(:, 1), '-o', 'DisplayName', 'Numerical solution');
% hold on;
% plot(t, besselj(1, t), '-o', 'DisplayName', 'J_{1}');
% legend;

A = [1 1 1  0  0  0; ...
     0 1 2 -1 -1 -1; ...
     0 1 4  0 -2 -4; ...
     0 1 8  0 -3 -12; ...
     0 1 16 0 -4 -32];
%      0 1 32 0 -5 -80];

b = [-1; -3; -9; -27; -81];

x = A\b;
disp((x(1) + x(3)));
disp(x(2));

A = [A, b];
disp(rref(A));

function dydt = bessel(t, y)
    dydt = [y(2); (-3./t).*y(2) - y(1)];
end