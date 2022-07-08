function H = HamHeis(n,Jx,Jy,Jz,hx,hy,hz)
% Heisenberg Hamiltonian
% See https://en.wikipedia.org/wiki/Quantum_Heisenberg_model

%construct tt_HamHeis right away as a tensor train

%% default arguments
if nargin < 1, n = 2; end     % number of spins
if nargin < 2, Jx = 0.1; end  % parameter of Hxx
if nargin < 3, Jy = 0.2; end  % parameter of Hyy
if nargin < 4, Jz = 0.3; end  % parameter of Hzz
if nargin < 5, hx = 0.5; end  % parameter of Hx
if nargin < 6, hy = 0.7; end  % parameter of Hy
if nargin < 7, hz = 0.9; end  % parameter of Hz

%% pauli matrices
sigmax = [0 1; 1 0];
sigmay = [0 -1i; 1i 0];
sigmaz = [1 0; 0 -1];

%% Hamiltonian: Hxx
%Sxx = zeros(2^n);
Sxx=tt_zeros(2^n,1)

for i = 1:n-1
    %add tensors to cell array as arg of mkron 
    k=tt_matrix(kron(sigmax,sigmax))
    ia=tt_eye(2^(i-1))
    ib=tt_eye(2^(n-i-1))
    size(ib)
    ic=kron(ia,k);
    c=kron(ic,ib)
    Sxx=Sxx+c
    
   % Sxx = Sxx + kron3(eye(2^(i-1)),kron(sigmax,sigmax),eye(2^(n-i-1)));
    %%write this as a tt of rank 1 and sum over rank one tensor train oseldets and how to 
    %construct TT by iterative svd core construction 
    %compare over original script for accuracy. Use the TT_tensor
    %constructor. Use version of TT_tensor where you pass cell array of
    %cores or copy source code of TT_tensor sontrsutor so that you do not
    %have to form the full matrix or vector directly. 
end
Hxx = -Jx*Sxx;

%% Hamiltonian: Hyy
Syy = zeros(2^n);
for i = 1:n-1
    Syy = Syy + kron3(eye(2^(i-1)),kron(sigmay,sigmay),eye(2^(n-i-1)));
end
Hyy = -Jy*Syy;

%% Hamiltonian: Hzz
Szz = zeros(2^n);
for i = 1:n-1
    Szz = Szz + kron3(eye(2^(i-1)),kron(sigmaz,sigmaz),eye(2^(n-i-1)));
end
Hzz = -Jz*Szz;

%% Hamiltonian: Hx
Sx = zeros(2^n); 
for i = 1:n
    Sx = Sx + kron3(eye(2^(i-1)),sigmax,eye(2^(n-i)));
end
Hx = -hx*Sx;

%% Hamiltonian: Hy
Sy = zeros(2^n); 
for i = 1:n
    Sy = Sy + kron3(eye(2^(i-1)),sigmay,eye(2^(n-i)));
end
Hy = -hy*Sy;

%% Hamiltonian: Hz
Sz = zeros(2^n); 
for i = 1:n
    Sz = Sz + kron3(eye(2^(i-1)),sigmaz,eye(2^(n-i)));
end
Hz = -hz*Sz;

%% full Hamiltonian
H = Hxx + Hyy + Hzz + Hx + Hy + Hz;

end


function K = kron3(A,B,C)
%KRON3  Kronecker tensor product.

K = kron(A,kron(B,C));

end