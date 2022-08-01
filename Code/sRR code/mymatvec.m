function y = mymatvec(A,x)
%small matvec counting function for optimization
global nbit
y = A*x;
nbit = nbit + 1;
end
