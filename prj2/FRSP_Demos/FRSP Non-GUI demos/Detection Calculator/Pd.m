function [P,T]= Pd( N, Pfa, SNR_dB, Swerling_case )

% Function Pd(  N, Pfa, SNR_dB, Swerling_case )
%
%  
% P = Pd(N,Pfa,SNR_dB,Swerling_case)
%   Calculates Pd for non-coherent pulse integration of Swerling target
%   models in Gaussian I/Q noise with a square-law detector.  Based on the
%   equations in the appendix of "Radar Target Detection- Handbook of Theory
%   and Practice" by Daniel P. Meyer and Herbert A. Mayer (M&M) (Academic
%   Press, 1973) to calculate the  threshold and probability of detection
%   for the four standard Swerling cases and the nonfluctuating case.  Also
%   included is the Albersheim approximation to the nonfluctuating case.
%
%   Inputs: N = number of pulses to noncoherently integrate (matrix)
%           Pfa = probability of false alarm (matrix, same dimensions as N)
%           SNR_dB = target signal to noise ratio, in dB (matrix, same dimensions as N)
%           Swerling_case = 0 or 5 (no fluctuation), 1, 2, 3, 4 (Swerling cases),
%                           or 6 (Albersheim's equation) (scalar)
%
% [P,T] = Pd(N,Pfa,SNR_dB,Swerling_case)
%   Also returns the matrix of thresholds, in volts normalized to the
%   noise level, for each combinaton of N, Pfa, and SNR_dB.  T is a
%   matrix of the same dimensions as N, Pfa, and SNR_dB.
%
% A test case:
%
% If N = [1 10 20] and Pfa = [1e-4 1e-6 1e-8] and SNR_dB = [0  3  6]
%        [5 15 25]           [2e-4 2e-6 2e-8]              [6  3  0]
%
% (i.e., Pfa, N and SNR_dB are the 2x3 matrices shown), and if the
% display format is set to 'short g' with the 'format' command,
% then you should get
%
%    P = [0.0036813     0.32791      0.99785 ]  for Swerling Case 0 or 5
%        [0.89305       0.70377      0.077277]  (nonfluctuating)
%
%    P = [0.01          0.32594      0.63331 ]  for Swerling Case 1
%        [0.54213       0.43493      0.22943 ]
%
%    P = [0.01          0.34918      0.98859 ]  for Swerling Case 2
%        [0.74559       0.63916      0.10643 ]
%
%    P = [0.0063839     0.34661      0.76136 ]  for Swerling Case 3
%        [0.63859       0.49941      0.21419 ]
%
%    P = [0.006564      0.34248      0.99698 ]  for Swerling Case 4
%        [0.8044        0.66632      0.092689]
%
%    P = [0.056577      0.34833      0.99992 ]  for Case 6 (Albersheim's Equation)
%        [0.9225        0.74649      0.1175  ]  (nonfluctuating, **linear** detector)
%
%  User functions called: threshold.m
%
%  Required MATLAB toolobxes: Signal Processing (for 'marcumq')
%
%   Mark A. Richards, June 2002
%  
%
%   Modified extensively from Douglas Dougherty's "meyerfun.m", version 1.2
%       (by Douglas Dougherty,  NSWC DD Code T45,  4/14/99)
%
%   Updated to use marcumq from SP Toolbox, June 2010

if nargin ~= 4 ,  error('Pd: Error: incorrect number of arguments');  end ;

if ( (size(N)~=size(Pfa)) | (size(N)~=size(SNR_dB)) | (size(Pfa)~=size(SNR_dB)) )
    error('Pd: Error: Dimensions of N, Pfa, and/or SNR_dB do not agree')
end


% calculate the detection thresholds as linear S/N ratio normalized to noise:
T = threshold(Pfa,N);

% some preliminary calculations

iterations = 100;             % approximation for infinity
SNR = 10.^( SNR_dB / 10 );    % convert db S/N to linear S/N ratio
[R,C] = size(N);              % number of rows and columns in input matrices
P = zeros(R,C);
k = 0: iterations ;          % approximation for k = 0 to infinity

%
% calculate the Pd based on the Swerling case chosen:
%
%
% Swerling Case 0 (or 5): Non-fluctuating target model
%
% Based on M&M eqn. (3-37)
%
if ( Swerling_case == 0  | Swerling_case == 5 )
    
    NSNR = N .* SNR ;
	A = sqrt(2 * NSNR) ;
	B = sqrt(2 * T) ;
    % can we vectorize this? besseli can be vectorized, what about Marcum?

	for r=1:R
        for c=1:C
            Trc=T(r,c); Nrc=N(r,c); NSNRrc=NSNR(r,c);
            g = 0;
            for k=2:Nrc
                g = g + ( (Trc/NSNRrc)^((k-1)/2) ) * besseli(k-1,2*sqrt(NSNRrc*Trc)) ;
            end
            P(r,c) = marcumq(A(r,c),B(r,c)) + exp(-Trc-NSNRrc) * g ;
        end
    end

%
% Swerling Case 1: Large (exponential fluctuation pdf),
%                  slow (scan-to-scan) RCS fluctuations
%
% Based on M&M (3-56)
%
elseif ( Swerling_case == 1 )
  
NSNR = N.*SNR;
A = 1+1./NSNR ;
AN = A.^(N-1) ;
P = 1 - gammainc(T,N-1) + AN.*gammainc(T./A,N-1).*exp(-T./(1+NSNR)) ;


%
% Swerling Case 2: Large (exponential fluctuation pdf),
%                  rapid (pulse-to-pulse) RCS fluctuations (or frequency agility)
%
% Based on M&M (3-60)
%
elseif ( Swerling_case == 2 )
    
	P = 1 - gammainc(T./(1+SNR),N);
  
  
%
% Swerling Case 3: Small (chi-square fluctuation pdf),
%                  slow (scan-to-scan) RCS fluctuations
%
% Based on M&M eqn. (3-69)
%
% Dougherty's code based on M&M eqn. (A-85) is preserved as an option
%
elseif ( Swerling_case == 3 )
    
%
% Approximation based on M&M (3-69).  This version is exact (therefore better than
% the Dougherty version) for N=1 and N=2, and seems to match on the test case for
% the other values of N.  And is vectorized and thus faster.  A version of the
% Dougherty code is commented out immediately below and may be used if preferred.
%
	NSNR = N.*SNR;
	APN = ( 1 + 2./NSNR ).^(N-2) ;
	AN2 = ( 1 + NSNR/2) ;
	P = APN .* ( 1 + (T./AN2) - 2*(N-2)./NSNR ) .* exp(-T./AN2) ;

%
% This commented-out block of code is a direct adaptation of Dougherty's
% code for Swerling 3, and is based on M&M eqn. (A-85)
%
%   logT = log(T) ;
%   G = 1 ./ ( 1 + (N .* SNR / 2) ) ;
%   G1 = 1 - G ;
%   
%   for r=1:R
%       for c=1:C
%           Trc=T(r,c); Nrc=N(r,c); logTrc=logT(r,c); Grc=G(r,c);  G1rc=G1(r,c);
%           fact = [0, cumsum(log(1:(Nrc + 1 + iterations)))] ;  % log factorials,
%                                                                %     offset by 1 for indexing
% 
%   % you can't use a negative index into fact, so do it manually for N = 1.
%           if Nrc == 1
%               f = 1 ;
%           else
%               f = fact(Nrc - 2 + 1) ;
%           end
%           if Nrc <= 2
%               m = 0 ;
%           else
%               m = 0:(Nrc - 2) ;
%           end
%           
%           a = exp(-Trc) * Grc * exp( (Nrc-1) * logTrc - f ) ;
%           b = exp(-Trc) * sum( exp( m * logTrc - fact(m + 1) ) ) ;
%           g = exp(-Grc * Trc) / ( G1rc^(Nrc - 2) ) ;
%           d = 1 - (Nrc - 2) * Grc / G1rc + Grc * Trc ;
%           e = 1 - exp(-Trc * G1rc) * sum( exp( m * (logTrc + log(G1rc)) - fact(m + 1) ) ) ;
%           P(r,c) = a + b + g * d * e ;
%       end
%   end

%
% Swerling Case 4: Small (chi-square fluctuation pdf), 
%                  rapid (pulse-to-pulse) RCS fluctuations
%
% Based on M&M Eqns. A-107 and A-111 and Dougherty's code
%
elseif ( Swerling_case == 4 )
    
	G = 1 ./ (1 + (SNR / 2) ) ;
	GT = G .* T ;
	logGT = log(GT) ;
	A = log( ( 1 - G ) ./ G ) ;

		for r=1:R
			for c=1:C
				Trc=T(r,c); Grc=G(r,c); GTrc=GT(r,c); Nrc=N(r,c); logGTrc=logGT(r,c); Arc=A(r,c);
				fact = [0, cumsum(log(1:(Nrc + 1 + iterations)))] ;  % log factorials,
                                                                   %     offset by 1 for indexing
				k = 0 : Nrc ;
				if Trc > Nrc *( 2-Grc )       % M&M eqn. (A-107)
					m = 0:(2 * Nrc - 1) ;
					b = exp(-GTrc) * cumsum( exp( m * logGTrc - fact(m+1) ) ) ;  % 2N-1-k +1(index) = 2N-k
					P(r,c) = Grc^Nrc * sum( exp( fact(Nrc+1) - fact(k+1) - fact(Nrc-k+1) + ...
                        (Nrc-k) * Arc ) .* b(2*Nrc-k) ) ;
                    
				else      % Trc <= Nrc *( 2-Grc ), so use M&M eqn. (A-111)
					m = (Nrc + iterations):-1:Nrc ;
					b(m) = exp(-GTrc) * cumsum( exp( m * logGTrc - fact(m+1) ) ) ;  % m=N -> b(N)
					P(r,c) = 1 - Grc^Nrc * sum( exp( fact(Nrc+1) - fact(k+1) - fact(Nrc-k+1) +  ...
                        (Nrc-k) * Arc ) .* b(2*Nrc-k) ) ;
				end
			end
		end


%
% Case 6: Albersheim's Equation for nonfluctuating, **linear** detector
%
elseif (Swerling_case == 6 )

    A = log(0.62./Pfa);
    Z = ( SNR_dB + 5.*log10(N) ) ./ ( 6.2 + ( 4.54 ./ sqrt(N+0.44) ) );
    B = ( 10.^Z - A ) ./ (1.7 + (0.12).*A );
    P = 1 ./ ( 1 + exp(-B) );
    
    
else
  error('Pd: Error: Invalid Swerling case.') ;
  
end

%
% Limit probabilities to the interval [0,1].  Controls
% minor roundoff errors,mainly in the Marcum function used
% for the Swerling 0/5 case
%

P = max(0,min(P,1));