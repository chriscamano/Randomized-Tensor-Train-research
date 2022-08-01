
function [S] = SRFT(s,n)
%Implementation of the Subsampled Random Fourier Transform 
%Designed and implemented by Chris Camaño and Roel Van Beumeen based on the
%algorithm description provided by Yuji Nakatsukasa and Joel Tropp in Feburary 2022 preprint:
% "Fast & Accurate Randomized Algorithms For Linear Systems and Eigenvalue Problems. 

% s: target low rank sketching dimension
% n: size of original matrix

S= zeros(s,n);
for  i=1:s
    %consider adding condition to make sure columns are linearly  independent.
    %consider drawing a random integer between 1 and n twice 
  
    z=zeros(n,1);
    z(randi([1,n]))=exp(1i*2*pi*rand);                  % combine the Steinhaus distribution matrix and randomly selected diagonal projection matrix at the same time
    S(i,:)=fft(z);
end

end

