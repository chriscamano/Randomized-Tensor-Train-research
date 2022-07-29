function [X,Lambda] = staticrarnoldi(A,k,tau,p)
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


if nargin <3 
    tau=1e-10;
end
n=size(A,1);
d=k;                                           %starting size for dimension of krylov subsspace      

%% Line 3                                       % init random starting vector    
%% Line 6                                       d-truncated Arnoldi iteration      

%% normal basis construction mode. 
q1=randn(n,1);
B(:,1) = q1/norm(q1);                                 %    Store first Arnoldi vector  
[B,d]=expandArnoldi(A,B,d);
for j=2:p
   [B,d]=expandArnoldi(A,B,d);
    s=4*d;                                      % target embedding dimension    
    %% Line 2                                     Create subsampled random fourier transform embedding (SRFT)
    S=SRFT(s,n);  
    C=S*B;                                      % Sketch basis C=S[b_1,...,b_dmax] 
    D=S*(A*B);                                 % Sketch D=S[m_1,...m_dmax]
    
    %% Line 9
    [U,T]=qr(C,0);                              %Compute thin QR of C  
    
    %% Line 11
    if(cond(T)>1/tau )
        B=B/T;                                  %whiten B  
    end
    
    %% Line 12                                  solve eigenproblem T^-1U^*Dy_i=\lambda_iy_i for i=1-d
    Mhat=T\(U'*D);                      %Form minimizer M hat via triangular substitution
    [y,Lambda]=eig(Mhat,'vector');              %invoke QR algorithm with vector output. ;                           
    %% check norms and recallibrate if krylov space is not large enough
    X=B*y;
    X=X./vecnorm(X);                            %Consider X=X(:,(1:k))./vecnorm(X(:,(1:k))); for 58/59
                  
    R=A*X(:,(1:k))-X(:,(1:k))*diag(Lambda(1:k));
    res=vecnorm(R);
end
[~,ii]=sort(abs(Lambda));             %sort eigenvalues by magnitidue 
ii=flip(ii);
Lambda=Lambda(ii); 

X=X(:,ii);
X=X(:,(1:k));
%Lambda=Lambda((1:k));
end
