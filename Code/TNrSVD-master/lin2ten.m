function ten=lin2ten(lin,n,varargin)
% ten=lin2ten(lin,n,d) or ten=lin2ten(lin,n)
% ------------------------------------------
% Converts a linear index into a tensor multi-index, which counts from 1.
% The decomposition of the multi-index is specified by n and d.
%
% ten           =	vector, tensor multi-index,
%
% lin           =	scalar, linear index,
%
% n     		=	vector, n(i) contains the ith dimension,
%
% d     		=	scalar, only provided when n is a scalar argument. The
%                   size of the axis (=dimension) is then n^d.
%
% Reference
% ---------
%
% 11-2016, Kim Batselier, Ngai Wong

% if isscalar(n)
%     d=varargin{1};
%     n=n*ones(1,d);
% elseif isvector(n)
%     d=length(n);
% end

d=length(n);

ten=zeros(1,d);
for i=d-1:-1:1
    if rem(lin,prod(n(1:i)))==0
        ten(i+1)=fix(lin/prod(n(1:i)));
        lin=lin-(ten(i+1)-1)*prod(n(1:i));
    else
        ten(i+1)=fix(lin/prod(n(1:i)))+1;
        lin=rem(lin,prod(n(1:i)));
    end
    
end
ten(1)=lin;


end