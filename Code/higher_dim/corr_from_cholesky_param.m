function [R, L] = corr_from_cholesky_param(h, d)
% Build a d x d correlation matrix R via the Huang-Ye-Wang (2025)
% Cholesky-based parameterisation.
%
% The mapping is
%   h_{i,j} in R   -->   theta_{i,j} = pi / (1 + exp(-h_{i,j})) in (0, pi)
%                  -->   S_{i,j} = sin(theta_{i,j}),  C_{i,j} = cos(theta_{i,j})
%                  -->   L (upper triangular, paper eq. 5)
%                  -->   R = L' * L
% R is symmetric, positive definite with unit diagonal for any h in R^{d(d-1)/2}.
%
% INPUT
%   h : (d*(d-1)/2 x 1) vector of unconstrained real parameters, stored
%       column-by-column for the upper-triangular off-diagonals:
%         d=3 -> h = [h_{1,2}; h_{1,3}; h_{2,3}]
%         d=4 -> h = [h_{1,2}; h_{1,3}; h_{2,3}; h_{1,4}; h_{2,4}; h_{3,4}]
%       In paper notation, these are the optimisation variables h_{i,j}.
%   d : matrix dimension (>= 2).
%
% OUTPUT
%   R : (d x d) correlation matrix.
%   L : (d x d) upper-triangular Cholesky factor, R = L' * L.
%
% Reference
%   Huang, Ye, Wang (2025), J. Choice Modelling 57, 100580, Section 2.2.

    h = h(:);
    expected = d*(d-1)/2;
    if numel(h) ~= expected
        error('corr_from_cholesky_param:badSize', ...
              'h must have d*(d-1)/2 = %d elements (got %d).', ...
              expected, numel(h));
    end

    % Step 1: unconstrained h -> bounded angles in (0, pi)
    theta_vec = pi ./ (1 + exp(-h));

    % Step 2: stash sin/cos in upper-triangular arrays
    S = zeros(d, d);
    C = zeros(d, d);
    idx = 0;
    for j = 2:d
        for k = 1:j-1
            idx      = idx + 1;
            S(k, j)  = sin(theta_vec(idx));
            C(k, j)  = cos(theta_vec(idx));
        end
    end

    % Step 3: build upper-triangular L via paper eq. (5)
    L = zeros(d, d);
    L(1, 1) = 1;
    for j = 2:d
        L(j, j) = S(1, j);                              % diagonal
        L(1, j) = prod(C(1:j-1, j));                    % top row
        for i = 2:j-1                                   % interior 1 < i < j
            L(i, j) = prod(C(1:j-i, j)) * S(j-i+1, j);
        end
    end

    % Step 4: correlation matrix
    R = L' * L;

    % Strict symmetry and exact unit diagonal (clean roundoff)
    R = 0.5 * (R + R');
    R(1:d+1:end) = 1;
end