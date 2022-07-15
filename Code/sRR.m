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

d=5
n=500
A=rand(n,n);

mssR(A,d,n)
%mssR(A,b,d,k,u,tau)
function [norms, x,lambda]= mssR(A,d,n)

%target embedding dimension
s=4*d;
%init w vectors
w=zeros(n,d);
B=zeros(n,d);
%init matrix AB
AB=zeros(n,d)

%% Line 2
% Draw subspace embedding S with s= 4d 

%create diagonal projector matrix onto s coordinates 
D=zeros(s,n);
for i=1:n
   D(randi([1,s]),i)=1; 
end

%create unitary discrete fourier transform 
F=dftmtx(n);

%create diagonal steinhaus matrix 
stein=diag(exp(1i*2*pi*rand(n,n)));
E=zeros(n,n);
E(1:(n+1):end)=stein;

S=sqrt(d/s)*D*F*E;
%% Line 3
%Starting vector w_1=randn(n,1)
w(:,1)=randn(n,1);


%% Line 4
%Normalize basis vector b_1=w_1/||w_1||_2 and apply m=AB_1
B(:,1)=(w(:,1)/norm(w(:,1)));
AB(:,1)=A*B(:,1);

%% Line 6
%Truncated Arnoldi w_j=,,,,
k=4;
for j=2:d
    t=eye(n);
    for i=1:k
     if(j-i<=0)
         break;
     end
     t=t-B(:,j-i)*ctranspose(B(:,j-i));
    end
   w(:,j)=t*AB(:,j-1);
   
%% Line 7
%Normalize and uppdate AB
   B(:,j)=(w(:,j)/norm(w(:,j)));
   AB(:,j)=A*B(:,j);
end

%% Line 8
%Sketch basis C=S[b_1,...,b_dmax] and reduce D=S[m_1,...m_dmnax]
C=S*B;
D=S*AB;

%% Line 9
%Compute thin QR factorization

[U T]=qr(C,"econ",0);

%% Line 9
%Thin QR over C=UT
%[U T]=QR(,,,,0)

if(k_2(T)>tol )
    %Line 11
    %Either Whiten B<-BT^-1 or stabilkize and solve 
end



%% Line 12
%solve eigenproblem T^-1U^*Dy_i=\lambda_iy_i for i=1-d

%% Line 13
%Form residual estimates ||Dy_i-\lambda_iCy_i||_2/||Cy_i||_2

%% Line 14
%Identify set I of indicies i where res is at most tol 

%% Line 15
%Compute x_i = By_i and normalize x_i = x_i/||x_i||_2 for i \ in I and
%output x(x_i, \lambda _ui) 
end
