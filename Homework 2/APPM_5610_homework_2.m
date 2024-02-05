clc;
clear;

n_max = 16;
max_iters = 100;

actual_evals = zeros(n_max-1, 1); 
dom_evals    = zeros(n_max-1, 1);
dom_evecs    = zeros(n_max-1, n_max-1);

for n = 2:n_max
    H = hilb(n);
     
    [V, D] = eigs(H);
    
    z_0 = zeros(length(H), 1); z_0(1) = 1;
    [lambda, v] = power_method(H, z_0, max_iters);

    actual_evals(n-1) = D(1, 1); 
    dom_evals(n-1)    = lambda;
    dom_evecs(1:n, n-1) = v;

end

disp("Actual eigenvalues and computed eigenvalues of Hilbert matrix");
fprintf("\t\tn\t\t\t\tActual eigenvalues\t Computed eigenvalues\n");
disp([(2:n_max)', actual_evals, dom_evals]);
fprintf("\n\nCorresponding computed eigenvectors\n");
disp(dom_evecs);

% Run the modified power iteration to obtain the smallest eigenvalues of
% the Hilbert matrix

actual_smallest_evals = zeros(n_max-1, 1);
smallest_evals        = zeros(n_max-1, 1);

disp("Actual smallest eigenvalues and computed smallest eigenvalues of Hilbert matrix");
fprintf("\t\tn\t\t\t\tActual smallest eigenvalues\t Computed smallest eigenvalues\n");

for n = 2:n_max
    H = hilb(n);
    [V, D] = eigs(H, n);

    z_0 = zeros(length(H), 1); z_0(1) = 1;
    [lambda_2, ~] = power_method_inv(H, z_0, max_iters);
    
    disp([n, lambda_2, min(diag(D))]);
end

% 
% A = diag([1, -1, 0.5]);
% [l_1, v_1] = power_method(A^2, max_iters);
% [l_2, v_2] = power_method(A^2, max_iters);
% 
% disp([v_1, v_2]);

% A = rand(9);
% A = triu(A) - diag(diag(A)) + diag([1, 1, 1, 1, 1, 1, 1, 1, 1e-2]);
% 
% V = zeros(length(A)-1, 1);
% D = zeros(length(A), 8);
% I = eye(size(A));
% 
% for i = 1:8
%     z_0 = I(:, i);
%     [lambda, v] = power_method(A, z_0, max_iters + 1000);
%     V(i) = lambda;
%     D(:, i) = v;
% end
% 
% 
% disp([A*D(:, 2), V(2)*D(:, 2)]);


function [lambda, v] = power_method(A, z_0, max_iters)

    for k = 1:max_iters
        z_k = A*z_0;
        q_k = z_k/norm(z_k);
        lambda = dot(q_k, A*q_k);
        z_0 = z_k;
    end

    v = q_k;
end

function [lambda, v] = power_method_inv(A, z_0, max_iters)
    [largest_e_val, ~] = power_method(A, z_0, max_iters);

    A_tilde = A - largest_e_val*eye(size(A));

    for k = 1:max_iters
        z_k = A_tilde*z_0;
        q_k = z_k/norm(z_k);
        lambda = dot(q_k, A_tilde*q_k);
        z_0 = z_k;
    end

    lambda = lambda + largest_e_val;
    v = q_k;
end