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

%% Driver Program


n=1000;
rng(101);
A=HamHeis(10);
%A=randn(n,n);
k=5; %number eigenpairs
tic;
[x,lambda]=rarnoldi(A,k);
toc

for i =1:k
    m=norm(A*x(:,i)-lambda(i)*x(:,i));
    if(m>1e-10)
       disp(m,i) 
    end
end

tic;
[V,D]=eigs(A,k);
toc
l=diag(D);
figure();
plot(l(1:size(lambda,1)),'o','MarkerFaceColor','yellow'); hold on
plot((lambda),'x','MarkerFaceColor','red'); hold on
title('Plot of eigenvalues computer by eig(A) vs eigenvalues computed by sketched RR');
xlabel('Real(\lambda_i)');
ylabel('imag(\lambda_i)');
hold off;


figure();
plot(abs(l(1:size(lambda,1))),'MarkerFaceColor','yellow'); hold on
plot(abs((lambda)),'MarkerFaceColor','red'); hold on
title('absolute value of predicted eigenvalues versus expected. ');
xlabel('\lambda_i');
ylabel('Magnitude');
hold off;
% figure();
% plot(l,'o','MarkerFaceColor','black'); hold on 
% plot(lambda,'x','MarkerFaceColor','cyan'); hold on 

%hold off;
%plot(real(diag(D)),imag(real(diag(D))),'*'); hold on 
%plot(real(diag(lambda)),imag(real (diag(lambda))),'o');
%[V,D]=eig(A);
 %plot(real(diag(D)),imag(real(diag(D))),'*'); hold on 
 %plot(real(diag(Lambda)),imag(real (diag(Lambda))),'o')
% ;
%
%% old prototype 
% function [x,lambda]= mssR(A,d,n)%tau,m
% tau=1e-10;
% s=4*d;                                      %target embedding dimension
% w=zeros(n,d);                               %init w vectors
% B=zeros(n,d);                               %init Basis
% AB=zeros(n,d);                              %init matrix AB
% S=zeros(s,n);
% 
% %% Line 2                                   Draw subspace embedding S with s= 4d                                      
% % D=zeros(s,n);                               %create diagonal projector matrix onto s coordinates 
% % C=randi([1,s],1,d);
% % 
% % for i=1:n
% %    D(randi([1,s]),i)=1; 
% % end
% 
% %F=fft(eye(n));
% %F=dftmtx(n);                                %create unitary discrete fourier transform
% %E=diag(exp(1i*2*pi*rand(n,1)));             %create diagonal steinhaus matrix 
% % S=sqrt(d/s)*D*F*E;                         %form subsampled random fourier transform  
% S= zeros(s,n);
% for  i=1:s
%     z=zeros(n,1);
%     z(randi([1,n]))=exp(1i*2*pi*rand);
%     S(i,:)=fft(z);
%     
% %     c=randi([1,n]);
% %     S(i,:)=fft(E(:,c));
% end
% %% Line 3
% w(:,1)=randn(n,1);                          %init starting vector
% 
% %% Line 4
% B(:,1)=(w(:,1)/norm(w(:,1)));               %Normalize and init basis vector
% AB(:,1)=A*B(:,1);                           %init first vector of AB based on guess 
% 
% %% Line 6                                   %k Truncated Arnoldi with m=10
% 
% ctrans=zeros(n);                            %preallocate for optimization
% a=zeros(n,1);                               %preallocate a 
% In=eye(n);                                  %preallocate identity to avoid reformation. 
% 
%                                             %   Arnoldi iteration with N-by-N matrix A and starting vector q1
%                                             %   (which need not have unit 2-norm).  For M < N it produces
%                                             %   an N-by-(M+1) matrix Q with orthonormal columns and an
%                                             %   (M+1)-by-M upper Hessenberg matrix H such that
%                                             %   A*Q(:,1:M) = Q(:,1:M)*H(1:M,1:M) + H(M+1,M)*Q(:,M+1)*E_M',
%                                             %   where E_M is the M'th column of the M-by-M identity matrix.
% m=5;
% 
% q1=rand(n,1);
% q1 = q1/norm(q1);
% 
% Q = zeros(n,m); Q(:,1) = q1;
% H = zeros(m,m-1);
% 
% %% arnoldi method for orthogonalization of B
% for i=1:m-1
%     z = A*Q(:,i);
%     delta = norm(z);
%     %% orthogonalize
%     h1 = Q(:,1:i)'*z;
%     z = z - Q(:,1:i)*h1;
%     %% reorthogonalize
%     h2 = 0;
%     if norm(z) < 0.5*delta
%         h2 = Q(:,1:i)'*z;
%         z = z - Q(:,1:i)*h2;
%     end
%     %% update H
%     H(1:i,i) = h1 + h2;
%     H(i+1,i) = norm(z);
%     %% expand subspace
%     Q(:,i+1) = z/H(i+1,i);
% end
% %% Line 7 
% B=Q;
% M=A*B;
% 
% %% Line 8
% C=S*B;                                      %Sketch basis C=S[b_1,...,b_dmax]  
% D=S*M;                                     %sketch D=S[m_1,...m_dmax]
% 
% %% Line 9
% [U,T]=qr(C,0);                              %Compute thin QR of C  
% %[U,T,p]=qr(C,0);                           %Compute thin QR of C with pivoting 
% 
% %% Line 11
% if(cond(T)>1/tau )
%     B=B/T;                               %whiten B  
% end
% 
% %% Line 12                                  solve eigenproblem T^-1U^*Dy_i=\lambda_iy_i for i=1-d
%                     
% Mhat=T\(U'*D);                              %compute t inverse via triangular substitution
% %Mhat=tinv*(ctranspose(U)*D);                %Form minimizer M hat
% [y,lambda]=eig(Mhat,'vector');                       %invoke QR algorithm with vector output. ;  
% % [~,ii]=sort(abs(lambda));                   %sort eigenvalues by magnitidue 
% % lambda=lambda(ii);
% % y=y(:,ii);
% 
% lambda=diag(lambda);
% 
% %% Line 13                               Form residual estimates ||Dy_i-\lambda_iCy_i||_2/||Cy_i||_2
%  res=zeros(d,1);
%  for i=1:d %make this k for k eigen pairs  
%     res(i)=norm(D*y(:,i)-lambda(i)*C*y(:,i))/norm(C*y(:,i));
%  end
%  % Line 14
%  I=res<=tau;                               %Identify set I of indices i where residual is at most tau
%  
% % %% Line 15                                %Compute x_i = By_i
% % 
% %                                           %normalize x_i = x_i/||x_i||_2
% %                                           %for i \ in I andoutput x(x_i,
% %                                           %\lambda _ui) for those with res
% %                                           %issues
% x=B*y;
% x=x./vecnorm(x);
% % for i=1:size(I)
% %     x(:,i)=x(:,i)/norm(x(:,1));
% % end
% end

