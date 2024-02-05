clc;
clear;

N = 3;
h = 1/(N+1);
K = 1;

I = eye(N);
T = diag(ones(N-1, 1), -1) + diag(-(4 + h^(2)*K)*ones(N, 1)) + ...
    diag(ones(N-1, 1), 1);

A = zeros(N*N);

A(1:3, 1:3) = T;
A(1:3, 4:6) = I;

A(4:6, 1:3) = I;
A(4:6, 4:6) = T;
A(4:6, 7:9) = I;

A(7:9, 4:6) = I;
A(7:9, 7:9) = T;

A = (1/h^(2))*A;