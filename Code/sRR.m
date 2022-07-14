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
A=HamHeis(5);




function [norms, x,lambda]= mssR(A,b,d,k,u=1e-10,tau)

%Line 2
%Draw subspace embedding S with s= 4d 

%Line 3
%Starting vector w_1=randn(n,1)
w1=randn(n,1);

%Line 
%Line 4
%Normalize basis vector b_1=w_1/||w_1||_2 and apply m=AB_1
for j=2:d
    %Line 6
    %Truncated Arnoldi w_j=,,,,
    
    %Line 7
    %Normalize basis vector b_j=w_j/||w_j||_2 and apply m=AB_j
    
end

%Line 8
%Sketch basis C=S[b_1,...,b_dmax] and reduce D=S[m_1,...m_dmnax]

%Line 9
%Thin QR over C=UT
[U T]=QR(,,,,0)

if(k_2(T)>tol )
    %Line 11
    %Either Whiten B<-BT^-1 or stabilkize and solve 
end



%Line 12
%solve eigenproblem T^-1U^*Dy_i=\lambda_iy_i for i=1-d

%Line 13
%Form residual estimates ||Dy_i-\lambda_iCy_i||_2/||Cy_i||_2

%Line 14
%Identify set I of indicies i where res is at most tol 

%Line 15
%Compute x_i = By_i and normalize x_i = x_i/||x_i||_2 for i \ in I and
%output x(x_i, \lambda _ui) 
end

