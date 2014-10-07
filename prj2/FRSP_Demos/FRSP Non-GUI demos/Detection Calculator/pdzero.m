function x = pdzero(SNRdB,N,pfa,swerlingcase,pd)

%
% Function used with MATLAB's 'fzero' to find the value
% of SNR in dB that gives a specified Pd, based on
% Meyer and Mayer (M&M) equations.
%
% Mark A. Richards, May 2010
%

x = Pd(N, pfa, SNRdB, swerlingcase) - pd;


