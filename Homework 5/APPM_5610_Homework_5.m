clc;
clear;

%%%
% Question 1
%%%
a = 2;
b = 1;
c = 1;

C = [a b; b c];

iters = 1;
C_k = chol_alg(C, iters);

disp(C_k);

%%%
% Problem 2
%%%
eps = 1e-10;
B = [2 eps; eps 1];

% Compute one step of the QR algorithm with and without a shift
iters = 1;
B_k_1 = qr_alg(B, 0, iters);
B_k_2 = qr_alg(B, 1, iters);

% disp(B_k_1);
% disp(B_k_2);

%%%
% Question 3
%%%
n = 99;
sub_diag  = ones(n, 1);
main_diag = 2*ones(n+1, 1);
A = diag(main_diag) + diag(sub_diag, -1) + diag(sub_diag, 1);

iters = 1000;
A_k = qr_alg(A, 0, iters);


% Compare the error between the eigenvalues returned by the QR iteration
% and MATLAB's eigs function
E_computed = diag(A_k);
E_actual   = eigs(A, n+1);

disp(abs(E_computed - E_actual));

function [A_prev] = chol_alg(A, iters)
    A_prev = A;
    for k = 1:iters
        A_prev = chol(A_prev); A_prev = A_prev*A_prev';
    end
end

function [A_prev] = qr_alg(A, mu, iters)
    A_prev = A;
    for k = 1:iters
        [Q_prev, R_prev] = qr(A_prev - mu*eye(size(A)));
        A_prev = R_prev*Q_prev + mu*eye(size(A));
    end
end