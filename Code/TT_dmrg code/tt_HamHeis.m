function H = tt_HamHeis(n,Jx,Jy,Jz,hx,hy,hz)
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
Sxx=ttm_zeros(2,n);
for i = 1:n-1
    
    if( i==1)
        %special handling for first core contraction
        Sxx=Sxx+tkron(tt_matrix(kron(sigmax,sigmax),eps,[2 2],[2 2]),tt_eye(2,(n-i-1)));
        
    elseif (i==n-1)
        Sxx=Sxx+tkron(tt_eye(2,(i-1)),tt_matrix(kron(sigmax,sigmax),eps,[2 2],[2 2]));
        
    else
        %special handling for Last core contraction
        Sxx=Sxx+tkron(tkron(tt_eye(2,(i-1)),tt_matrix(kron(sigmax,sigmax),eps,[2 2],[2 2])),tt_eye(2,(n-i-1)));
    end
    %Sxx = Sxx + kron3(eye(2^(i-1)),kron(sigmax,sigmax),eye(2^(n-i-1))); 
end
Hxx = -Jx*Sxx;

%% Hamiltonian: Hyy
Syy = ttm_zeros(2,n);
for i = 1:n-1
    if( i==1)
        %special handling for first core contraction
        Syy=Syy+tkron(tt_matrix(kron(sigmay,sigmay),eps,[2 2],[2 2]),tt_eye(2,(n-i-1)));
        
    elseif (i==n-1)
        Syy=Syy+tkron(tt_eye(2,(i-1)),tt_matrix(kron(sigmay,sigmay),eps,[2 2],[2 2]));
        
    else
        %special handling for Last core contraction
        Syy=Syy+tkron(tkron(tt_eye(2,(i-1)),tt_matrix(kron(sigmay,sigmay),eps,[2 2],[2 2])),tt_eye(2,(n-i-1)));
    end
    
    
    % Syy = Syy + kron3(eye(2^(i-1)),kron(sigmay,sigmay),eye(2^(n-i-1)));
end
Hyy = -Jy*Syy;

%% Hamiltonian: Hzz
Szz = ttm_zeros(2,n);
for i = 1:n-1
    if( i==1)
        %special handling for first core contraction
        Szz=Szz+tkron(tt_matrix(kron(sigmaz,sigmaz),eps,[2 2],[2 2]),tt_eye(2,(n-i-1)));
        
    elseif (i==n-1)
        Szz=Szz+tkron(tt_eye(2,(i-1)),tt_matrix(kron(sigmaz,sigmaz),eps,[2 2],[2 2]));
        
    else
        %special handling for Last core contraction
        Szz=Szz+tkron(tkron(tt_eye(2,(i-1)),tt_matrix(kron(sigmaz,sigmaz),eps,[2 2],[2 2])),tt_eye(2,(n-i-1)));
    end
   
    %Szz = Szz + kron3(eye(2^(i-1)),kron(sigmaz,sigmaz),eye(2^(n-i-1)));
end
Hzz = -Jz*Szz;

%% Hamiltonian: Hx
Sx = ttm_zeros(2,n);
for i = 1:n
    if( i==1)
        %special handling for first core contraction
        Sx=Sx+tkron(tt_matrix(sigmax,eps,[2 2]),tt_eye(2,(n-i)));
        
    elseif (i==n)
        Sx=Sx+tkron(tt_eye(2,(i-1)),tt_matrix(sigmax));
        
    else
        %special handling for Last core contraction
        Sx=Sx+tkron(tkron(tt_eye(2,(i-1)),tt_matrix(sigmax,eps,[2 2])),tt_eye(2,(n-i)));
    end
   
    %Sx = Sx + kron3(eye(2^(i-1)),sigmax,eye(2^(n-i)));
end
Hx = -hx*Sx;

%% Hamiltonian: Hy
Sy = ttm_zeros(2,n);
for i = 1:n
    
    if( i==1)
        %special handling for first core contraction
        Sy=Sy+tkron(tt_matrix(sigmay,eps,[2 2]),tt_eye(2,(n-i)));
    elseif (i==n)
        Sy=Sy+tkron(tt_eye(2,(i-1)),tt_matrix(sigmay));
        
    else
        %special handling for Last core contraction
        Sy=Sy+tkron(tkron(tt_eye(2,(i-1)),tt_matrix(sigmay,eps,[2 2])),tt_eye(2,(n-i)));
    end
    
    %Sy = Sy + kron3(eye(2^(i-1)),sigmay,eye(2^(n-i)));
end
Hy = -hy*Sy;

%% Hamiltonian: Hz
Sz = ttm_zeros(2,n);
for i = 1:n
    if( i==1)
        %special handling for first core contraction
        Sz=Sz+tkron(tt_matrix(sigmaz,eps,[2 2]),tt_eye(2,(n-i)));
    elseif (i==n)
        Sz=Sz+tkron(tt_eye(2,(i-1)),tt_matrix(sigmaz));
    else
        %special handling for Last core contraction
        Sz=Sz+tkron(tkron(tt_eye(2,(i-1)),tt_matrix(sigmaz,eps,[2 2])),tt_eye(2,(n-i)));
    end
    
    %Sz = Sz + kron3(eye(2^(i-1)),sigmaz,eye(2^(n-i)));
end
Hz = -hz*Sz;

%% full Hamiltonian

H=Hxx+Hyy+Hzz+Hx+Hy+Hz;
%H = Hxx + Hyy + Hzz + Hx + Hy + Hz;

end
