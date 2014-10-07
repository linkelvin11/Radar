function x = pfazero(T,N,Pfa);
%
% Function used with MATLAB's 'fzero' to find the value
% of threshold T that gives a specified Pfa, based on
% Meyer and Mayer (M&M) equations.
%
% M&M (2.17) gives Pfa = 1 - inverse_gamma_fcn(T/sqrt(N),N-1).
% However, the MATLAB definition of the inverse gamma
% function differs from the Pearson's form used in M&M.
% The equivalent statement in MATLAB is
% Pfa = 1 - gammainc(T,N).
%
% Thus in MATLAB 1 - gammainc(T,N) - Pfa = 0.
% MATLAB's fzero will call this expression with T varying and
% N & Pfa fixed until it finds the value of T that
% makes this equation hold (within the tolerance set by fzero).
%
% Mark A. Richards, May 2002
% Modified from code developed by Douglas Dougerty, NSWC
%

x = 1 - gammainc(T,N) - Pfa;