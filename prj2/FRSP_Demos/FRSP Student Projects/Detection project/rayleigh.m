function x=rayleigh(N,m)
%
% rayleigh
%
% M-file to generate a random sequence with a Rayleigh pdf
% and specified mean.
%
%  N = desired sequence length
%  m = desired mean of sequence
%
% M. A. Richards
% February 1997
%


% Generate uniforms, then transform.

x=rand(N,1);
s=4*(m^2)/pi;
x=sqrt(-s*log(x));

% Various tests:
% Compute sample mean and standard deviation

%sample_mean=mean(x)
%sample_standard_deviation=std(x)

%fprintf('Theoretical standard deviation = %g\n',sqrt((4-pi)/pi)*m)


% Show 100-bin histogram and also a theoretical Rayleigh
% curve with the desired mean 'm'.

%figure(2)
%[counts bins]=hist(x,100);
%pdf=((pi/(2*m^2))*bins).*exp(-(pi/(4*m^2))*(bins.^2));
%bar(bins,counts/sum(counts))
%hold
%plot(bins,pdf/sum(pdf),'r')
%xlabel('x')
%ylabel('density')
%hold off
