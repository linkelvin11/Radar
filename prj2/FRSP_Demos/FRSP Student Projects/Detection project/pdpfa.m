%
% PdPfa
%
%  Pd/Pfa via analysis and simulation for
%  nonfluctuating target and Rayleigh interference
%
%  Written by  M. A. Richards
%  February 1997
%  Updated April 2000
%  Update again March 2004
%  And again Nov. 2006
% 
clear, hold off
close all
format compact
j = sqrt(-1);

% Form complex interference signal with Rayleigh amplitude
% and unit mean

N = input('Enter desired noise sequence length: ');
snrdB = input('Enter desired SNR (dB): ');

beta = 2/sqrt(pi)*1;
cnoise = beta/sqrt(2)*randn(N,1)+j*beta/sqrt(2)*randn(N,1);
ray = abs(cnoise);
mean_ray = mean(ray)  % just checking to make sure mean ~= 1

% Now create another sequence with a target
% added to each sample at the specified SNR.

npow = beta^2;
S = sqrt((10^(snrdB/10))*npow);

snoise = cnoise + S;

sray = abs(snoise);

% The matched filter output is simply the input multiplied
% by conj(S) (but S is real here); the test statistic is the magnitude of that
% quantity (see Fig. 6-5 in text).  However, since we usually don't
% know S in advance, I will assume S = 1 for scaling the data and
% establishing the threshold.

% Now let's set up a series of threshold settings to achieve
% desired Pfa's using analytic formula from text

Pfa = [1e-7 2e-7 5e-7 1e-6 2e-6 5e-6 1e-5 4e-5 1e-4 4e-4 1e-3  1e-2]';
T = beta*sqrt(-log(Pfa));

% Now let's measure the Pfa with these thresholds and see how
% the results stack up against the design goal

for k=1:length(T)
  Pfa_sim(k) = sum( ray>T(k) )/N;
end
fprintf('\n\n Threshold    Theoretical Pfa    Simulated Pfa\n')
for k=1:length(T)
  fprintf('%9.6g          %5.3g         %5.3g\n',T(k),Pfa(k),Pfa_sim(k))
end

figure
loglog(Pfa_sim,Pfa);xlabel('Theoretical Pfa');ylabel('Simulated Pfa');
grid


% Now let's predict Pd for the same thresholds analytically,
% and again compare against simulation

% Do two versions, one using Marcum Q to predict Pd, the other
% using Albersheim with N=1 (and rearranged to compute Pd instead
% of SNR)

% Some caution needed here in using the marcum function. For the
% simulation, where I was mimicking the actual processing of data, I
% couldn't assume I knew S, so I just used S = 1 in the threshold
% calculation.  To use the Marcum function, however, I do have to plug in
% the actual SNR, which means using the actual value of S, and I need to
% use the threshold I would have computed if I had known S.  That threshold
% would be larger than my actual threshold by a factor of S, so where the
% Marcum function uses sqrt(2)*T(k)/S/beta for the second argument, I just
% want to use sqrt(2)*T(k)/beta, because my threshold is already smaller by
% the factor S.

for k=1:length(T)
   A=log(0.62/Pfa(k));
   Z=snrdB/(6.2+(4.54/sqrt(1.44)));
   B=(10^Z-A)/(1.7+0.12*A);
   Pd_Alb(k)=1/(1+exp(-B));
   
   Pd_Q(k)=marcum(sqrt(2)*S/beta,sqrt(2)*T(k)/beta,0.001);
   
   Pd_sim(k) = sum( sray>T(k) )/N;
end
fprintf('\n\n Threshold    Pd (Marcum Q)    Pd (Albersheim)     Simulated Pd\n')
for k=1:length(T)
  fprintf('%9.6g        %6.4g              %6.4g              %6.4g\n',T(k),Pd_Q(k),Pd_Alb(k),Pd_sim(k))
end

% Finally, let's plot Pd vs. Pfa in  "receiver operating characteristic"
% (ROC) curve.  Plot both predicted and simulated Pd.

figure
semilogx(Pfa,Pd_Alb,Pfa,Pd_Q,Pfa,Pd_sim); grid; xlabel('Pfa'); ylabel('Pd');
legend('Albersheim','Marcum Q','Simulated',-1)
title(['SNR=',int2str(snrdB),' dB'])

figure
hist([ray,sray],100)

