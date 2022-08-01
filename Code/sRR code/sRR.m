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
A=HamHeis(11);                                  %generate a Heisenberg Hamiltonian size 2^n
%A=randn(n,n);

k=3;                                            %number eigenpairs


%% basic time test

%srr test
% tic;
% [x,lambda]=rarnoldi(A,k);
% toc

% for i =1:k
%    norm(A*x(:,i)-lambda(i)*x(:,i))
% end
% 
% %normal eigs test
% tic;
% [V,D]=eigs(A); 
% toc
%%time testing plot for dense

% c=zeros(3,1);
% d=zeros(3,1);
% ax=zeros(3,1);
% for i = 400:400:40000 
%     disp(i)
%     A=randn(n,n);
%     tic;
%     [x,lambda]=rarnoldi(A,k);
%     c(i/400)=toc;
%     tic;
%     [V,D]=eig(A);
%     d(i/400)=toc;
%     ax(i/400)=i;
% end
% 
% plot(ax,c,'-o'); hold on
% plot(ax,d,'-x'); hold on
% legend('srr','eigs()');
% title('srr vs eigs() on hamiltonian data');
% xlabel('order of hamiltonian. 2^n');
% ylabel('total compilation time (s)');



% %  [x,lambda]=staticrarnoldi(A,k,1000);

format short e	
%
%% time testing plot for hamiltonians
c=[0 0 0 0 0];
d=[0 0 0 0 0];

for m = 5:13 
    disp(m)
    A=HamHeis(m);
    tic;
    [x,lambda]=rarnoldi(A,k);
    c(m)=toc;
    tic;
    [V,D]=eigs(A,k);
    d(m)=toc;
end

ytickformat('%.1f')
plot(c,'-o','DisplayName','srr','LineWidth',5); hold on
plot(d,'-x','DisplayName','eigs','LineWidth',5); hold on
hold off
legend('srr','eigs()','FontSize',33);
title('srr vs eigs() on Hamiltonian of Size 2^n','FontSize',33);
xlabel('Order of Hamiltonian. 2^n','FontSize',33);
ylabel('Total Compilation Time (s)','FontSize',33);



% 
