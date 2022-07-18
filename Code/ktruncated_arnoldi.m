function [Q,H] = k_arndoli(A,q1,m)
%ARNOLDI    Arnoldi iteration
%   [Q,H] = ARNOLDI(A,q1,M) carries out M iterations of the
%   Arnoldi iteration with N-by-N matrix A and starting vector q1
%   (which need not have unit 2-norm).  For M < N it produces
%   an N-by-(M+1) matrix Q with orthonormal columns and an
%   (M+1)-by-M upper Hessenberg matrix H such that
%   A*Q(:,1:M) = Q(:,1:M)*H(1:M,1:M) + H(M+1,M)*Q(:,M+1)*E_M',
%   where E_M is the M'th column of the M-by-M identity matrix.

n = length(A);
if nargin < 3, m = n; end
q1 = q1/norm(q1);
Q = zeros(n,m); Q(:,1) = q1;
H = zeros(min(m+1,m),n);
for i=1:m
    z = A*Q(:,i);
    for j=1:i
        H(j,i) = Q(:,j)'*z;
        z = z - H(j,i)*Q(:,j);
    end
    if j < n
       H(j+1,j) = norm(z);
       if H(j+1,j) == 0, return, end
       Q(:,j+1) = z/H(j+1,j);
   end
end

