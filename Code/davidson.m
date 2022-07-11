function [V, E, FLAG] = davidson(A, varargin)
%DAVIDSON Davidson iterative method for a few of the eigenvalues and 
%eigenvectors of real symmetric and diagonally dominant matrices
% E = DAVIDSON(A) returns the lowest eigenvalue of a real symmetric and
% diagonally dominant matrix A.
% E = DAVIDSON(A, I) returns the lowest few eigenvalues specified by index
% array I.
% E = DAVIDSON(..., 'largest') returns the largest (few) eigenvalue(s).
% E = DAVIDSON(..., 'rough') uses a less strict convergence criterion.
% [V, E] = DAVIDSON(...) returns the eigenpairs.
% [V, E, FLAG] = DAVIDSON(...) also returns a convergence flag. FLAG is the
% number of eigenvalues that are not converged.

% Reference:
% [1] Davidson, E.R., "The iterative calculation of a few of the lowest 
% eigenvalues and corresponding eigenvectors of large real-symmetric 
% matrices", J. Comput. Phys. 17, 87-94 (1975)

debug = 0;

if ~ismatrix(A) || size(A,1) ~= size(A,2) || ~isreal(A)
    error('Input matrix is invalid.')
end

if any(~isfinite(A(:)))
    error('Input matrix contains NaN or Inf.')
end

% check whether input matrix is symmetric
if any(abs(A-A') > 10*eps)
    error('Input matrix should be symmetric.')
end

opt_conv = 'strict';
if ~isempty(varargin) && any(strcmp(varargin, 'rough'))
    opt_conv = 'rough';
    varargin(strcmp(varargin, 'rough')) = [];
end

opt_eval = 'smallest';
if ~isempty(varargin) && any(strcmp(varargin, 'largest'))
    opt_eval = 'largest';
    varargin(strcmp(varargin, 'largest')) = [];
    A = -A;
end

if ~isempty(varargin)
    I = varargin{1};
    varargin(1) = [];
else
    I = 1;
end

if ~isempty(varargin)
    error('Option is not recognized.');
end

% check whether the input index array is valid
if ~isvector(I) || any(mod(I, 1)) || any(I < 1)
    error('Index array is invalid.')
end

sz = size(A,1);
if max(I) > sz
    error('Index exceeds matrix size.')
end

[I_sort, idx_sort_I] = sort(I);
if I_sort(end) > nthroot(sz, 3)
    warning('Computing general eigenvalues is inefficient.')
end

nvals = length(I); % number of eigenvalues to compute
E = zeros(nvals, 1); % eigenvalues
V = zeros(sz, nvals); % eigenvectors

% convergence criteria
crit1 = @(q) norm(q) < 1e-8;
crit2 = @(q, alphaM_end) norm(q) < 1e-4 && abs(alphaM_end) < 1e-12;

% max number of iterations for each eigenvalue
iter_count_max = 1000;
if debug
    hist_conv_crit = zeros(nvals, 1);
    hist_iter_count = zeros(nvals, 1);
    hist_qnorm = zeros(iter_count_max, nvals);
    hist_alpha = zeros(iter_count_max, nvals);
end

% subspace size limits for the k-th eigenvalue
sz_subspace_min = @(k) k + 5;
sz_subspace_max = @(k) k + 20;

if sz_subspace_max(I_sort(end)) > sz
    warning('Matrix size is too small. Use built-in dense algorithm instead.')
    [V_all, E_all] = eig(A, 'vector');
    [E_all, idx_sort_E_all] = sort(E_all);
    V_all = V_all(:, idx_sort_E_all);
    
    if strcmp(opt_eval, 'largest')
        E_all = -E_all;
    end
    
    E = E_all(I);
    V = V_all(I);
    
    if nargout == 1
        V = E;
    end
    
    FLAG = sum(isnan(E));
    return
end

% preallocate space for b and Ab
b = zeros(sz, sz_subspace_min(I_sort(end)));
Ab = zeros(sz, sz_subspace_max(I_sort(end)));

for ival = 1 : nvals
    
    k = I_sort(ival); % seek the k-th eigenvalue
    
    % subspace size limits for the k-th eigenvalue
    sz_min = sz_subspace_min(k);
    sz_max = sz_subspace_max(k);
    
    % guess space
    if ival == 1
        % generate from scratch
        b(:, 1:sz_min) = guess_init(diag(A), sz_min);
        sz_b = sz_min;
    else
        sz_alpha = size(alphaM, 1);
        b(:, 1:sz_alpha) = b(:, 1:sz_alpha) * alphaM;
        sz_b = sz_alpha;
        
        if sz_b < sz_min
            b(:, 1:sz_min) = [b(:, 1:sz_b), guess_init(diag(A), sz_min-sz_b)];
            sz_b = sz_min;
        end
        [b(:, 1:sz_b), ~] = qr(b(:, 1:sz_b), 0);
    end
    
    Ab(:, 1:sz_b) = A * b(:, 1:sz_b);
    
    iter_count = 0;
    while true
        if iter_count > iter_count_max
            E(ival) = nan;
            V(:, ival) = nan;
            if debug
                hist_iter_count(ival) = iter_count;
            end
            warning('Eigenvalue not converged.')
            break
        end
        Asub = b(:, 1:sz_b)' * Ab(:, 1:sz_b);
        Asub = (Asub + Asub') / 2;
        
        [alphaM, lambda] = eig(Asub, 'vector');
        [lambda, idx_sort_lambda] = sort(lambda);
        alphaM = alphaM(:, idx_sort_lambda);

        if sz_b == sz_max
            b(:, 1:sz_min) = b(:, 1:sz_max) * alphaM(:, 1:sz_min);
            sz_b = sz_min;
            [b(:, 1:sz_b), ~] = qr(b(:, 1:sz_b), 0);
            Ab(:, 1:sz_b) = A * b(:, 1:sz_b);
            continue
        end
        
        % residual vector
        q = ( Ab(:, 1:sz_b) - lambda(k)*b(:, 1:sz_b) ) * alphaM(:, k);
        
        is_conv = false;
        if strcmp(opt_conv, 'strict')
            is_conv = crit1(q);
        elseif strcmp(opt_conv, 'rough')
            is_conv = crit1(q) || crit2(q, alphaM(end, k));
        end
        
        if is_conv
            E(ival) = lambda(k);
            V(:, ival) = b(:, 1:sz_b) * alphaM(:, k);
            if debug
                hist_conv_crit(ival) = 1*crit1(q) + 2*crit2(q, alphaM(end, k));
                hist_iter_count(ival) = iter_count;
            end
            break
        end
        
        d = q ./ (lambda(k) - diag(A));
        if any(~isfinite(d))
            d(isfinite(d)) = 0;
            d(~isfinite(d)) = 1;
        end
        [b_new, ~] = qr([b(:, 1:sz_b), d], 0);
        b(:, sz_b+1) = b_new(:, end);
        Ab(:, sz_b+1) = A*b_new(:, end);
        sz_b = sz_b + 1;
        
        iter_count = iter_count + 1;
        
        if debug
            hist_qnorm(iter_count, ival) = norm(q);
            hist_alpha(iter_count, ival) = abs(alphaM(end, k));
        end
    end
end

if strcmp(opt_eval, 'largest')
    E = -E;
end

[~, idx_revsort_I] = sort(idx_sort_I);
E = E(idx_revsort_I);
V = V(:, idx_revsort_I);

if nargout == 1
    V = E;
end

FLAG = sum(isnan(E));

end

function [ b ] = guess_init(A_diag, sz_sub)
% generate initial guess space based on the the diagonal elements
sz = length(A_diag);
[~, idx_small] = mink(A_diag, sz_sub);
b = 0.001 * (rand(sz, sz_sub) - 0.5);
%b = 0.001 * ones(sz, sz_sub);
b(idx_small + (0:sz_sub-1)'*sz) = 1;
[b, ~] = qr(b, 0);
end

