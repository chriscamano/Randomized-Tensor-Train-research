
function [B,d] = rBlockKrylov(A,B,del_p)
%EXPANDARNOLDI Summary of this function goes here
%   A : current krylov subspace
%   m : amount to expand Krylov subspace by

bsize=50;
p=size(B,2)/bsize; %krylov depth
H=zeros(p*bsize+del_p*bsize,p*bsize);%compute size of current Krylov space

for i=p:(p+del_p)-1
    q=i*bsize;
    
    %% special case of i=1
    if(i==1)
        z = A*B(:,1:q);
        delta = z./norm(z);
         %% orthogonalize
        h1 = B(:,1:q)'*z; %h1
        z = z - B(:,1:q)*h1;
        
        %% reorthogonalize
        h2 =zeros(p,bsize);
        if norm(z) < 0.5*delta
            h2 = B(:,1:q)'*z;
            z = z - B(:,1:q)*h2;
        end
        
        %% update H
        H(1:q,1:q) = h1 + h2;
        H(i+bsize,1:q) = norm(z);
        
        %% expand subspace
        B(:,q+1:q+bsize) = z/H(i+1,i);
    else
        z = A*B(:,(i-1)*bsize:q-1);
        delta = norm(z); %
        
        %% orthogonalize
        h1 = B(:,1:q)'*z; %h1
        z = z - B(:,1:q)*h1;
        
        %% reorthogonalize
        h2 =zeros(q,bsize);
        if norm(z) < 0.5*delta
            h2 = B(:,1:q)'*z;
            z = z - B(:,1:q)*h2;
        end
        
        %% update H
        H(1:q,(i-1)*bsize:q-1) = h1 + h2;
        H(i+bsize,1:q) = norm(z);
        
        %% expand subspace
        B(:,q+1:q+bsize) = z/H(i+1,i);
    end
end
d=p+del_p;                              %update size in main program

end

% function [B d AB] = rBlockKrylov(A,B,k,m)%bsize,p
% %RBLOCKKRYLOV Summary of this function goes here
% %   Detailed explanation goes here
% %bsize: block size
% %p : krylov depth
% % k number of vectors for reorthogonalization
% %d amount to expand by
%
% d=size(B,2);
% n=size(A,1);
% k=2;
%
%
%
% AB=zeros(n,d);
% AB(:,1)=A*B(:,1);
% %truncated arnoldi
% I=eye(n);
%
%
% for j = d+1:(d+m)-1
%     w_j=0;
%     for i=1:k
%         if((j-i)<=0)
%         break;
%         end
%         w_j= w_j-B(:,(j-i))* B(:,(j-i))';
%     end
%
%     w_j=(I+w_j)*(A*B(:,(j-1)));
%     b_j=w_j/norm(w_j);
%
%     B(:,j)= b_j ;
%     AB(:,j)=A*b_j;
% end
%
% %GS consider double.
% for j=1: d+m-1
%     v=B(:,j);
%     for i=1:j-1
%         R(i,j)=Q(:,i)'*B(:,j);
%         v=v-R(i,j)*Q(:,i);
%     end
%     R(j,j)=norm(v);
%     Q(:,j)=v/R(j,j);
% end
% d=d+m;
% B=Q;
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
% %b: block size
% %if vargin <3 p=4; % krylov depth
% %k=2 %iterations
% % n=size(A,1);
% % %% Line 3                                       % init random starting vector
% % w_1=randn(n,bsize);
% % w_1 = w_1./norm(w_1);
% %
% % %% Line 6                                       d-truncated Arnoldi iteration                            % pre-allocate krylov subspace B
% % b_1 = w_1;
% % B_temp={b_1};
% % bstack=0;
% % for j=2:p
% %     w_j=zeros(n,bsize);
% %     %% Form identity
% %     I=eye(n);
% %     for i=1:k
% %         if(j-i<=0)                             %per the paper when i>=0 break
% %           continue;
% %         else
% %             bstack=bstack-B_temp{j-i}*B_temp{j-i}' ;
% %         end
% %     end
% %
% %     w_j=(I-bstack)*(A*B_temp{j-1});
% %     v=w_j;
% %     for i=1:bsize
% %
% %        r_ii=norm(v(:,i));
% %        q_i=v(:,i)/r_ii;
% %
% %        for j=(i+1):n
% %
% %            r_ij=q_i'*v(:,j);
% %            v(:,j)=v(:,j)-r_ij*q_i;
% %
% %        end
% %
% %     end
% %
% %
% %     B_temp=[B_temp;v];
% % end
% % B=B_temp{1};
% % for i=2:p
% %     B=B_temp{i};
% % end
% % % %% Line 6                                   %k Truncated Arnoldi with k=4;
% % % k=4;
% % % for j=2:d
% % %     t=eye(n);                               %temp starting matrix
% % %     for i=1:k
% % %      if(j-i<=0)                             %per the paper when i>=0 break
% % %          break;
% % %      end
% % %      t=t-B(:,j-i)*ctranspose(B(:,j-i));     %Subtract the computed value from I
% % %     end
% % %    w(:,j)=t*AB(:,j-1);                      %multiply by AB
% % %
% % % %% Line 7
% % %    B(:,j)=(w(:,j)/norm(w(:,j)));            %Normalize
% % %    AB(:,j)=A*B(:,j);                        % update AB
% % % end
% end