% Author:Christian Camano 
% Implementation of Als Algorithm from page 9 of:
% "The Alternating Linear Scheme for Tensor Optimization in the TT format.
% Sebastian Holtz, Thorsten Rohwedder, et al. 


%examples will be tested on an order 4 tensor train and implemented with
%TT_toolbox. 

%Initialize tensor with dimensions 3x4x5
disp('------------------------------------------------')
A=randn(8,2);
T=reshape(A,[2 2 2 2]);

%construct an order 4 tensor train
tt_y=tt_tensor(T);
d=size(core(tt_y));
d=d(1);

disp('Size of TT')
size(tt_y)

%Line 2.
tt_x=tt_tensor(tt_y);
ps=tt_x.ps;

%Line 3:
disp('╔═════════════════════════╗')
disp( "  Right to left Orthogonalization Algorithm")
disp('╚═════════════════════════╝')

for n=d:-1:2
    %---------------------------------------
    fprintf(' ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
    fprintf(' Current Core index: %d\n',n)
    fprintf(' current core dimensions: %s\n',mat2str(size(core(tt_x,n))))
    fprintf(' ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
    %---------------------------------------
    
    %Line 4
    H = unfold_H(core(tt_x,n));
    [Q,R] = qr(H',0); %0= economy mode tall skinny qr
    
    % update value of tensor array with orthogonal basis from QR
    % formula for array manip taken from tt_toolbox docs
    tt_x.core(ps(n) : ps(n+1)-1) = reshape(Q',1,[]); 
  
    %line 5
    V = unfold_V(core(tt_y,n-1))*R'; %V* R^T
    tt_x.core(ps(n-1) : ps(n)-1) = reshape(V,1,[]);
    
    %---------------------------------------
    fprintf(' loop iteration %d complete \n',n)
    %---------------------------------------

end
%space











disp('╔═════════════════════════╗')
disp( "           Error analysis step")
disp('╚═════════════════════════╝')



for n=1:4
 %---------------------------------------
 fprintf(' ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
 fprintf('   Orthogonality testing, Core:%d\n',n)
 fprintf(' ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
 %---------------------------------------
    A=unfold_H(core(tt_x,n))'*unfold_H(core(tt_x,n))
    %heatmap(unfold_H(core(tt_x,n)),'Colormap',bone)
    %heatmap(A,'Colormap',bone)
    
end
size(core(tt_x,2))

heatmap(unfold_H(core(tt_x,2))'*unfold_H(core(tt_x,2)),'Colormap',bone)
fprintf(' ━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp   ('       TT norm test')
fprintf(' ━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')

norm(tt_y-tt_x)

fprintf(' ━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp('      Matrix norm test')
fprintf(' ━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')

norm(full(tt_y)-full(tt_x))

fprintf(' ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp   ('   Element comparison: (1,1,1,1)')
fprintf(' ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')

tt_x(1,1,1,1)
tt_y(1,1,1,1)
disp('------------------------------------------------')

function H=unfold_H(core)
core=squeeze(core);
  if(ismatrix(core))
      H=core;
    return;
  end
%   A=size(core);
%   if(A(1)~= max(size(core)))
%       index_max=find(A==max(size(core)));
%       temp=A(1);
%       A(1)=max(size(core));
%       A(index_max)=temp;
%   end
% 
%   [x,y,z]=deal(A(1),A(2),A(3))

%4,2,2
%2,2,4 
  [x y z]=size(core);
  H=reshape(core,[x,y*z]);
end

function V=unfold_V(core)
core=squeeze(core);
  if(ismatrix(core))
      V=core;
    return;
  end
  
  [x y z]=size(core);
  V=reshape(core,[x*y z]);
end
