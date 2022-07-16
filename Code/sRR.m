% Implementation of the sketched rayleigh-ritz method described by authors
% Yuji Nakatsukasa and Joel Tropp in Feburary 2022 preprint:
% "Fast & Accurate Randomized Algorithms For Linear Systems and Eigenvalue Problems. 
% implementation by Chris Camano Lawrence Berkeley National Laboratory. 


% %Algorithm details
% Pre: Square matrix and vector , basis dimension d, nuymber k of vector for 
% partial orthogonalization, stability tolerance tol =O(u^-1), convergence tolerance 
% \tau.
% 
% Post: Approximate eigenpairs (x_i,\lambda_i) such that Ax_i \approx \lambda_ix_i. 
% Estimated residual norms \hat{r}_{est,i}

d=5;
n=500;
A=rand(n,n);

mssR(A,d,n)
%mssR(A,b,d,k,u,tau)



function [norms, x,lambda]= mssR(A,d,n)
s=4*d;                                      %target embedding dimension
w=zeros(n,d);                               %init w vectors
B=zeros(n,d);                               %init Basis
AB=zeros(n,d);                               %init matrix AB
%% Line 2                                   Draw subspace embedding S with s= 4d 
                                            
D=zeros(s,n);                               %create diagonal projector matrix onto s coordinates 
for i=1:n
   D(randi([1,s]),i)=1; 
end                                         
F=dftmtx(n);                                %create unitary discrete fourier transform 
                                            %create diagonal steinhaus matrix 
stein=diag(exp(1i*2*pi*rand(n,n)));
E=zeros(n,n);
E(1:(n+1):end)=stein;

S=sqrt(d/s)*D*F*E;                          %form subsampled random fourier transform
%% Line 3
w(:,1)=randn(n,1);                          %init starting vector
%% Line 4
B(:,1)=(w(:,1)/norm(w(:,1)));               %Normalize and init basis vector
AB(:,1)=A*B(:,1);                           %init first vector of AB based on guess 

%% Line 6                                   %k Truncated Arnoldi with k=4;
k=4;
for j=2:d
    t=eye(n);                               %temp starting matrix 
    for i=1:k                   
     if(j-i<=0)                             %per the paper when i>=0 break 
         break;
     end
     t=t-B(:,j-i)*ctranspose(B(:,j-i));     %Subtract the computed value from I 
    end
   w(:,j)=t*AB(:,j-1);                      %multiply by AB
   
%% Line 7
   B(:,j)=(w(:,j)/norm(w(:,j)));            %Normalize 
   AB(:,j)=A*B(:,j);                        % update AB
end

%% Line 8
C=S*B;                                      %Sketch basis C=S[b_1,...,b_dmax]  
D=S*AB;                                     %sketch D=S[m_1,...m_dmax]

%% Line 9
[U,T]=qr(C,0);                              %Compute thin QR of C  
%[U,T,p]=qr(C,0);                           %Compute thin QR of C with pivoting 

%% Line 11
if(cond(T)>eps^-1 )
    B=B*T^-1;                               %whiten B  
end

%% Line 12                                  solve eigenproblem T^-1U^*Dy_i=\lambda_iy_i for i=1-d

B=speye(size(T));           
tinv = B(:,:)/T;                            %compute t inverse via triangular substitution

Mhat=tinv*(ctranspose(U)*D);                %Form minimizer M hat

for i=1 : 20                                %invoke QR algorithm 
    [Q,R]=qr(Mhat);
    Mhat=R*Q;
end
theta=diag(Mhat);
eigs(Mhat)
%% Line 13
%Form residual estimates ||Dy_i-\lambda_iCy_i||_2/||Cy_i||_2

%% Line 14
%Identify set I of indicies i where res is at most tol 

%% Line 15
%Compute x_i = By_i and normalize x_i = x_i/||x_i||_2 for i \ in I and
%output x(x_i, \lambda _ui) 
end

