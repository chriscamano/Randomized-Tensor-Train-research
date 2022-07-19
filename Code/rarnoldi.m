function [X,Lambda] = rarnoldi(A,k,tau)
% Implementation of the sketched rayleigh-ritz method described by authors
% Yuji Nakatsukasa and Joel Tropp in Feburary 2022 preprint:
% "Fast & Accurate Randomized Algorithms For Linear Systems and Eigenvalue Problems. 
% implementation by Chris Camano Lawrence Berkeley National Laboratory. 
%% Algorithm details
% Pre: Square matrix and vector , basis dimension d, nuymber k of vector for 
% partial orthogonalization, stability tolerance tol =O(u^-1), convergence tolerance 
% \tau.
% 
% Post: Approximate eigenpairs (x_i,\lambda_i) such that Ax_i \approx \lambda_ix_i. 
% Estimated residual norms \hat{r}_{est,i}



if nargin <3 tau=1e-10;
% shift for mgs **
% A
% k
%
%


d=5;
n=size(A,1);
s=4*d;                                        % target embedding dimension                         
%% Line 2                                     Create subsampled random fourier transform embedding (SRFT)
S=SRFT(s,n);                                
%% Line 3                                     % init random starting vector    
q1=rand(n,1);
q1 = q1/norm(q1);

%% Line 6                                    d-truncated Arnoldi iteration
d=6;


B = zeros(n,d);                               % pre-allocate krylov subspace B
B(:,1) = q1;                                  % Store first Arnoldi vector 
H = zeros(d,d-1);                             % pre-allocate H
z=zeros(n,1);

for i=1:d-1  
    z = A*B(:,i);
    delta = norm(z);
    %% orthogonalize
    h1 = B(:,1:i)'*z;
    z = z - B(:,1:i)*h1;
    %% reorthogonalize
    h2 = 0;
    if norm(z) < 0.5*delta
        h2 = B(:,1:i)'*z;
        z = z - B(:,1:i)*h2;
    end
    %% update H
    H(1:i,i) = h1 + h2;
    H(i+1,i) = norm(z);
    %% expand subspace
    B(:,i+1) = z/H(i+1,i);
end
%% Line 7 
 %M=A*B; small improvement about .001

 %% Line 8
C=S*B;                                      % Sketch basis C=S[b_1,...,b_dmax]  
D=S*(A*B);                                  % Sketch D=S[m_1,...m_dmax]
%% Line 9
[U,T]=qr(C,0);                              %Compute thin QR of C  
%[U,T,p]=qr(C,0);                           %Compute thin QR of C with pivoting 
%% Line 11
if(cond(T)>1/tau )
    B=B/T;                                  %whiten B  
end
%% Line 12                                  solve eigenproblem T^-1U^*Dy_i=\lambda_iy_i for i=1-d
                    
                                            
Mhat=T\(U'*D);                              %Form minimizer M hat via triangular substitution
[y,Lambda]=eig(Mhat,'vector');              %invoke QR algorithm with vector output. ;  

% [~,ii]=sort(abs(lambda));                 %sort eigenvalues by magnitidue 
% lambda=lambda(ii);
% y=y(:,ii);

%Lambda=diag(Lambda);                            
X=B*y;
X=X./vecnorm(X);




%helper fucntion for S and also CHECK 

% compute residuals of oredered ritz pairs and update the dimension of the
% krylov subspace in the event that there is a bad res on one of the
% eigenvalues. 


% only project back into the orginal subspace once all of the ritz pairs
% have convereged 
end

