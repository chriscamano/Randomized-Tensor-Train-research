function [S] = SRFT(s,n)
S= zeros(s,n);
%consider making s dynamic 
for  i=1:s
    %consider adding condition to make sure columns are linearly
    %independent.
    
    %consider drawing a random integer between 1 and n twice 
    %
    z=zeros(n,1);
    z(randi([1,n]))=exp(1i*2*pi*rand);
    S(i,:)=fft(z);
end

end

