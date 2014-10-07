%
% ROC
%
%  Comparative receiver operating characteristic
%  curves for log-normal and Rayleigh interference
%
%  Written by  M. A. Richards
%  February 1997
%
clear, hold off
format compact
j = sqrt(-1);


% Form two complex interference signals, one with Rayleigh
% amplitude, the other with log-normal amplitude but
% the same mean and variance in both cases.

N = input('Enter desired noise sequence length: ');
ray = rayleigh(N,1);
mean_ray = mean(ray)
std_ray = std(ray)
nsig = sqrt((4-pi)/pi);  % std dev of a Rayleigh with mean = 1
logn = log_norm2(N,1,2*nsig);
mean_logn = mean(logn)
std_logn = std(logn)

nray = ray.*exp(j*2*pi*rand(N,1));
nlogn = logn.*exp(j*2*pi*rand(N,1));

% Now create another version of each sequence with a target
% added to each sample at a 20 dB SNR.  Then do a linear
% detection on each of the 4 sequences.

npow = nsig^2;
A = sqrt(20*npow);
snr = 10*log10(A^2/nsig^2)

sray = nray + A;
slogn = nlogn + A;

mnray = abs(nray);
msray = abs(sray);
mnlogn = abs(nlogn);
mslogn = abs(slogn);
mean_nray = mean(mnray)
std_nray = std(mnray)
mean_nlogn = mean(mnlogn)
std_nlogn = std(mnlogn)



% Now let's set up a series of 100 threshold settings and
% evaluate # of threshold crossings of absolute value of
% each of the 4 signals at each setting.

tmax=max(mslogn);
for k=1:100;
  t=k*tmax/100;
  pfa_ray(k)=sum(mnray>t)/N;
  pd_ray(k) = sum(msray>t)/N;
  pfa_logn(k)=sum(mnlogn>t)/N;
  pd_logn(k)=sum(mslogn>t)/N;
end

semilogx(pfa_logn,pd_logn,'y');
hold
semilogx(pfa_ray,pd_ray,'r');
hold off
xlabel('Pfa')
ylabel('Pd')
grid

