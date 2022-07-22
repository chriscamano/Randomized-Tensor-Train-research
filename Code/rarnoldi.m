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
n=size(A,1);
d=10;                                         %starting size for dimension of krylov subsspace                             
%% Line 3                                     % init random starting vector    
q1=rand(n,1);
q1 = q1/norm(q1);
%% Line 6                                    d-truncated Arnoldi iteration
B = zeros(n,d);                               % pre-allocate krylov subspace B
B(:,1) = q1;                                  % Store first Arnoldi vector 
H = zeros(d,d-1);                             % pre-allocate H
%%add to seperate routine where you add additional vectors
%% add seperate function that expands arnoldi subspace adding d to the size of the subspace. 
% parameter for addded subspace size f( small subspace, no new vectors) 

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

for j=1:100
    s=4*d;                                        % target embedding dimension    
    %% Line 2                                     Create subsampled random fourier transform embedding (SRFT)
    S=SRFT(s,n);  
    C=S*B;                                      % Sketch basis C=S[b_1,...,b_dmax]  
    %D=S*(A*B);                                  % Sketch D=S[m_1,...m_dmax]
    %% Line 9
    [U,T]=qr(C,0);                              %Compute thin QR of C  
    %% Line 11
    if(cond(T)>1/tau )
        B=B/T;                                  %whiten B  
    end
    %% Line 12                                  solve eigenproblem T^-1U^*Dy_i=\lambda_iy_i for i=1-d
    Mhat=T\(U'*(S*(A*B)));                              %Form minimizer M hat via triangular substitution
    [y,Lambda]=eig(Mhat,'vector');              %invoke QR algorithm with vector output. ;  


    %experiment with other sorting.
%     [~,ii]=sort(abs(Lambda));                 %sort eigenvalues by magnitidue 
%     Lambda=Lambda(ii); 
%     y=y(:,ii);

                          
    %% check norms and recallibrate if krylov space is not large enough
    X=B*y;
    X=X./vecnorm(X);

    R=A*X(:,(1:k))-X(:,(1:k))*diag(Lambda(1:k));
    res=vecnorm(R);
    
    if(any(res>tau))
        [B,d]=expandArnoldi(A,B,d,d);
    else
       break; 
    end
end



%truncate final output. 
X=X(:,(1:k));
Lambda=Lambda((1:k));




%helper fucntion for S and also ``CHECK 

% compute residuals of oredered ritz pairs and update the dimension of the
% krylov subspace in the event that there is a bad res on one of the
% eigenvalues. 


% only project back into the orginal subspace once all of the ritz pairs
% have convereged 
end

