function y = mymatvec(A,x)

global nbit
y = A*x;
nbit = nbit + 1;

end
