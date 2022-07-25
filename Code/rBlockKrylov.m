function [B] = rBlockKrylov(A,bsize,p)
%RBLOCKKRYLOV Summary of this function goes here
%   Detailed explanation goes here

%b: block size
%if vargin <3 p=4; % krylov depth
k=2 %iterations 
n=size(A,1);
%% Line 3                                       % init random starting vector    
w_1=randn(n,bsize);
w_1 = w_1./norm(w_1);

%% Line 6                                       d-truncated Arnoldi iteration                            % pre-allocate krylov subspace B
b_1 = w_1; 
B_temp={b_1};
bstack=0;
for j=2:p
    w_j=zeros(n,bsize);
    %% Form identity 
    I=eye(n);
    for i=1:k
        if(j-i<=0)                             %per the paper when i>=0 break 
          continue;
        else
            bstack=bstack-B_temp{j-i}*B_temp{j-i}' ;
        end
    end
    
    w_j=(I-bstack)*(A*B_temp{j-1});
    v=w_j;
    for i=1:bsize
        
       r_ii=norm(v(:,i));
       q_i=v(:,i)/r_ii;
       
       for j=(i+1):n
           
           r_ij=q_i'*v(:,j);
           v(:,j)=v(:,j)-r_ij*q_i;
           
       end
       
    end
    
    
    B_temp=[B_temp;v];
end
B=B_temp{1};
for i=2:p
    B=B_temp{i};
end
% %% Line 6                                   %k Truncated Arnoldi with k=4;
% k=4;
% for j=2:d
%     t=eye(n);                               %temp starting matrix 
%     for i=1:k                   
%      if(j-i<=0)                             %per the paper when i>=0 break 
%          break;
%      end
%      t=t-B(:,j-i)*ctranspose(B(:,j-i));     %Subtract the computed value from I 
%     end
%    w(:,j)=t*AB(:,j-1);                      %multiply by AB
%    
% %% Line 7
%    B(:,j)=(w(:,j)/norm(w(:,j)));            %Normalize 
%    AB(:,j)=A*B(:,j);                        % update AB
% end
end

