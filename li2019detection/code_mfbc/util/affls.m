% Solve the following minimization problem:
%
% min_X = norm(B - A * X)
%
% subject to colsums(X) = 1
%
%
% Input:
%  
% A - coefficient matrix
% B - matrix of right hand sides
%
% Output:
%
% Solution matrix X.
%

function [X] = affls(A, B)

k = size(A, 2);
nrhs = size(B, 2);

C = [A'*A ones(k , 1); ones(1, k) 0];
R = [A' * B; ones(1, nrhs)];

Xtilde = C \ R;
X = Xtilde(1:end-1, :);

end

