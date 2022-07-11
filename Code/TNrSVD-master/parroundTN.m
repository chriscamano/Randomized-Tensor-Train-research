function a=parroundTN(a)

d=size(a.n,1);
% right to left
for i=d:-1:2
    [B,T]=depar(reshape(a.core{i},[a.n(i,1),prod(a.n(i,2:end))])');
    a.n(i,1)=size(B,2);
    a.core{i}=reshape(B',a.n(i,:));
    a.core{i-1}=reshape(reshape(a.core{i-1},[prod(a.n(i-1,1:end-1)),a.n(i-1,end)])*T',[a.n(i-1,1:end-1) size(B,2)]);
    a.n(i-1,end)=size(B,2);        
end

% left to right
% for i=1:d-1
%     [B,T]=depar(reshape(a.core{i},[prod(a.n(i,1:3)),a.n(i,end)]));
%     a.n(i,end)=size(B,2);
%     a.core{i}=reshape(B,a.n(i,:));
%     a.core{i+1}=reshape(T*reshape(a.core{i+1},[a.n(i+1,1),prod(a.n(i+1,2:end))]),[size(B,2),a.n(i+1,2:end)]);
%     a.n(i+1,1)=size(B,2);    
% end
end