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


if nargin <3 
    tau=1e-10;
end
n=size(A,1);
d=k;                                           %starting size for dimension of krylov subsspace      
%bsize=50;
%;



%% Line 3                                       % init random starting vector    
%% Line 6                                       d-truncated Arnoldi iteration      

%% normal basis construction mode. 
q1=randn(n,1);
B(:,1) = q1/norm(q1);                                 %    Store first Arnoldi vector  
[B,d]=expandArnoldi(A,B,d);

%% Random Block Subspace mode 
% B=zeros(n,bsize);
% b_1=randn(n,bsize);
% B(:,1:bsize)=orth(b_1);
% [B d]=rBlockKrylov(A,B,d,bsize);
%% Expand subspace until tolerance on k eigen pairs is met. 
for j=1:100
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
   
    %X=X(:,(1:k));
    %Lambda=Lambda((1:k));
    
                     
    R=A*X(:,(1:k))-X(:,(1:k))*diag(Lambda(1:k));
    res=vecnorm(R);
    

%     semilogy(res);hold on
%     pause
    
    if(any(res>tau))
        %[B d ]=rBlockKrylov(A,B,d,bsize);
        [B,d]=expandArnoldi(A,B,d);
    else
        
       break; 
    end
end
[~,ii]=sort(abs(Lambda));             %sort eigenvalues by magnitidue 
ii=flip(ii);
Lambda=Lambda(ii); 

X=X(:,ii);
X=X(:,(1:k));
%Lambda=Lambda((1:k));
end
%Lambda=Lambda((1:k));

% figure();
% plot(abs(Lambda));hold on
% plot(abs(l));
% figure();
% 
% plot(l(1:20),'o','MarkerFaceColor','black'); hold on 
% plot(Lambda(1:20),'x','MarkerFaceColor','cyan'); hold on 



% for i = 1:200
%   plot(l(1:i),'s','MarkerFaceColor','green'); hold on
%   plot(Lambda(1:i),'s','MarkerFaceColor','blue'); hold on
% 
% end

%plot(l(1),'s','MarkerFaceColor','magenta'); hold on
% plot(abs(l),'o','MarkerFaceColor','pink');hold on 
% plot(abs(Lambda),'x','MarkerFaceColor','blue'); hold on
% plot(abs(Lambda(1)),'s','MarkerFaceColor','red'); hold on
% plot(abs(l(1)),'s','MarkerFaceColor','red'); hold onend
