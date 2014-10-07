function   x = git_chirp( T, W, p )
%CHIRP      generate a sampled chirp signal
%    X = git_chirp( T, W, <P> )
%      X:  N=pTW samples of a "chirp" signal
%            exp(j(W/T)pi*t^2)   -T/2 <= t < +T/2
%      T:  time duration from -T/2 to +T/2
%      W:  swept bandwidth from -W/2 to +W/2
%    optional:
%      P:  samples at P times the Nyquist rate (W)
%            i.e., sampling interval is 1/(PW)
%            default is P = 1
%
if nargin < 3
  p = 1;
end    
J = sqrt(-1);
%--------------
delta_t = 1/(p*W);
N = round( p*T*W );    %--same as T/delta_t, but rounded
nn = [0:N-1]';
% x = exp( J*pi*W/T * (delta_t*nn - T/2).^2 );  % old version
x = exp( J*pi*W/T * (delta_t*nn - (N-1)/2/p/W).^2 );  % symmetric version


% even older version
%%%%% alf = 1/(2*p*p*T*W);
%%%%% git_chirp = exp( J*2*pi*alf*((nn-N/2).*(nn-N/2)) );
