clc;
clear;

t_0 = 0;
t_f = 1;
h = (t_f - t_0)/10;
y_0 = [1; 999/10];

y = method(@odefun, t_0, t_f, h, y_0);

V = [1, 1; ...
     0, 999/10];
D = [-100, 0; ...
     0, -1/10];

y_actual = @(t, y_0) V*expm(t*D)*inv(V)*y_0;

disp(y);
disp(y_actual(t_f, y_0));

function [f] = odefun(t, y)
    f = [-100, 1; ...
         0, (-1/10)]*y;
end

function [y] = method(odefun, t_0, t_f, h, y_0)
    y0 = y_0;

    % Compute y_1 using one iteration of forward Euler
    y1 = y0 + h*odefun(t_0, y0);

    while t_0 <= t_f
        y2 = 4*y1 - 3*y0 - 2*h*odefun(t_0, y0);
        y0 = y1;
        y1 = y2;
        t_0 = t_0 + h;
    end

    y = y1;
end