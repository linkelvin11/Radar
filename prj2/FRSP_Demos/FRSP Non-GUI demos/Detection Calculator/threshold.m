function T=threshold(Pfa,N)

% Takes matrix of Pfa values and a matrix of N 
% values and finds correspdoning matrix of threshold values.
%
% If Pfa and N are scalar, it finds one threshold T.
%
% If Pfa is a vector of values and N is a vector of
% the same length filled with the same value N, this
% gives you vector of thresholds for different Pfas 
% and fixed N.  Similarly, you can make the Pfa vector
% a constant and the N vector changing, and get thresholds
% for fixed Pfa and varying N.
%
% Finally, by using a matrix, you can get thresholds for a
% combination of different Pfas and Ns.
%
% A test case, verified one entry at a time using Dougherty's code:
%
% If Pfa = [1e-4 1e-6 1e-8] and N = [1 10 20]
%          [2e-4 2e-6 2e-8]         [5 15 25]
%
% (i.e., Pfa and N are the 2x3 matrics shown), then you should get
%
%    T = [ 9.2103   32.71   55.956]
%        [16.898    39.981  62.75 ]
%
% with the display format set to 'short G'.
%
% User functions called: pfazero.m
%
% Mark A. Richards, June 2002.
% Modified from code developed by Douglas Dougerty, NSWC

if size(Pfa) ~= size(N)
    error('threshold: Error: Pfa and N not the same dimensions');
    return
end

% use Fehlner's initial approximations for T, the detection threshold.
% The units are in volts, relative to the noise floor.

n = log(2)*ones(size(Pfa)) ./ Pfa;   % Fehlner's "false alarm number"

T=zeros(size(Pfa));
smallN = N<=12;
bigN = ~smallN;

T(smallN) = N(smallN) .* ( 1 + 2.2 * log10(n(smallN)) ./ ...
    (N(smallN).^( 2/3 + 0.015 * log10(n(smallN)) )) );
T(bigN) = N(bigN) .* ( 1 + 2.2 * log10(n(bigN)) ./ (N(bigN).^( 2/3 + 0.015 * log10(n(bigN)) )) );


% Use the MATLAB fzero function to refine the initial estimate above and
% find the value of T that gives the desired Pfa.  Unfortunately, fzero
% is a scalar function, so we have to iterate through the combinations
% of Pfa, N, and initial T to compute the trhesolds one at a time.

[R,C]=size(Pfa);
for r = 1:R
    for c = 1:C
        T(r,c)=fzero(@pfazero,T(r,c),optimset('disp','off'), N(r,c), Pfa(r,c) );
    end
end
