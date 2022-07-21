function [S] = SRFT(s,n)
S= zeros(s,n);

for  i=1:s
    z=zeros(n,1);
    z(randi([1,n]))=exp(1i*2*pi*rand);
    S(i,:)=fft(z);
end

end

