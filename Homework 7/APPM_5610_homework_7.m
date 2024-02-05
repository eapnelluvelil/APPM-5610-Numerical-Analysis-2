clc;
clear;
close;

a = -5;
b = 5;
n = 400;
x = linspace(a, b, n);
[X, Y] = meshgrid(x, x);
Z = X + 1i*Y;
inds = zeros(size(Z));

% 2-step BDF stability polynomial

% 3-step BDF stability polynomial
BDF_3_step_poly = @(r) [1 - (6/11)*r; ...
                        -18/11; ...
                        9/11; ...
                        -2/11];


for i = 1:n
    for j = 1:n
        % Compute roots
        r = roots(BDF_3_step_poly(Z(i, j)));
        if (abs(r) < ones(size(r)))
            inds(i, j) = 1;
        end
    end
end

figure;
plot(real(Z(inds == 1)), imag(Z(inds == 1)));