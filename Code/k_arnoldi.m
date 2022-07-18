function [lam,res] = k_arnoldi(Afun,n,m,k,shift)

if nargin < 5, shift = inf; end

%% initialize
V = zeros(n,m+1);
V(:,1) = ones(n,1)/sqrt(n);
H = zeros(m+1,m);

%% Arnoldi iteration
for j = 1:k
    if n > 1000, fprintf('arnoldi iteration %i\n',j); end
    [V,H] = arnoldi_step(Afun,V,H,j);
end

%% eigenvalues and residual
if nargout < 2
    lam = eig(H(1:m,1:m));
    if ~isinf(shift), lam = shift + 1./lam; end
    lam = sort(lam);
else
    [S,lam] = eig(H(1:m,1:m),'vector');
    if ~isinf(shift), lam = shift + 1./lam; end
    [lam,idx] = sort(lam);
    res = abs(H(m+1,m)*S(end,idx))'; %%check for res computations before calling it/ 
end

end


function [V,H] = arnoldi_step(Afun,V,H,j)

%% last vector
w = V(:,j);

%% shift-and-invert
w = Afun(w);
delta = norm(w);

%% orthogonalize
h1 = V(:,1:j)'*w;
w = w - V(:,1:j)*h1;

%% reorthogonalize
h2 = 0;
if norm(w) < 0.5*delta
    h2 = V(:,1:j)'*w;
    w = w - V(:,1:j)*h2;
end

%% update H
H(1:j,j) = h1 + h2;
H(j+1,j) = norm(w);

%% expand subspace
V(:,j+1) = w/H(j+1,j);

end
