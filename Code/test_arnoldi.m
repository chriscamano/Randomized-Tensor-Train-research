% test_arnoldi

%% parameters
n = 400;
maxit = 120;

%% generate matrix
rng(100);
A = randn(n);

%% eigenvalues
lam = eig(A);

%% arnoldi
[ritz,res] = arnoldi(@(x) A*x,n,maxit);

%% output
idx = res < 1e-2;
[ritz(idx),res(idx)]

%% figure
figure;
plot(real(lam),imag(lam),'*'); hold on
plot(real(ritz(idx)),imag(ritz(idx)),'o');
legend('lambda','ritz');
