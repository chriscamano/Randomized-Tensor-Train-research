function H = HamHeis(n,Jx,Jy,Jz,hx,hy,hz)
% Heisenberg Hamiltonian
% Author: Roel Van Beumeen
% See https://en.wikipedia.org/wiki/Quantum_Heisenberg_model

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
Sxx = zeros(2^n);
for i = 1:n-1
    Sxx = Sxx + kron3(eye(2^(i-1)),kron(sigmax,sigmax),eye(2^(n-i-1)));
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
