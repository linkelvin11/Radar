function [ ncc ] = mclutter( CN )
%MCLUTTER Creates clutter at a given CNR
%   Detailed explanation goes here
hold off
format compact
J = sqrt(-1);
% close all

% Get root file name for saving results

% file=input('Enter root file name for data and listing files: ','s');


% form radar chirp pulse

T = 10e-6;     % pulse length, seconds
W = 10e6;      % chirp bandwidth, Hz
fs = 12e6;     % chirp sampling rate, Hz; oversample by a little

% fprintf('\nPulse length = %g microseconds\n',T/1e-6)
% fprintf('Chirp bandwidth = %g Mhz\n',W/1e6)
% fprintf('Sampling rate = %g Msamples/sec\n',fs/1e6)
s = git_chirp(T,W,fs/W);
% plot((1e6/fs)*(0:length(s)-1),[real(s) imag(s)])
% title('Real and Imaginary Parts of Chirp Pulse')
% xlabel('time (usec)')
% ylabel('amplitude')
% grid

Np = 20;              % 20 pulses
jkl = 0:(Np-1);       % pulse index array
PRF = 10.0e3;         % PRF in Hz
PRI = (1/PRF);        % PRI in sec
T_0 = PRI*jkl;        % relative start times of pulses, in sec
g = ones(1,Np);       % gains of pulses
T_out = [12 40]*1e-6; % start and end times of range window in sec
T_ref = 0;            % system reference time in usec
fc = 10e9;            % RF frequency in Hz; 10 GHz is X-band

% fprintf('\nWe are simulating %g pulses at an RF of %g GHz',Np,fc/1e9)
% fprintf('\nand a PRF of %g kHz, giving a PRI of %g usec.',PRF/1e3,PRI/1e-6)
% fprintf('\nThe range window limits are %g to %g usec.\n', ...
%     T_out(1)/1e-6,T_out(2)/1e-6)

% Compute unambiguous Doppler interval in m/sec
% Compute unambiguous range interval in meters

vua = 3e8*PRF/(2*fc);
rmin = 3e8*T_out(1)/2;
rmax = 3e8*T_out(2)/2;
rua = 3e8/2/PRF;

% fprintf('\nThe unambiguous velocity interval is %g m/s.',vua)
% fprintf('\nThe range window starts at %g km.',rmin/1e3)
% fprintf('\nThe range window ends at %g km.',rmax/1e3)
% fprintf('\nThe unambiguous range interval is %g km.\n\n',rua/1e3)

% Define number of targets, then range, SNR, and
% radial velocity of each.  The SNR will be the actual SNR of the target in
% the final data; it will not be altered by relative range.

Ntargets = 4;
del_R = (3e8/2)*( 1/fs )/1e3;                   % in km
ranges = [2 3.8 4.4 4.4]*1e3;               % in km
SNR =    [-3 5 10 7];               % dB
vels = [-0.4 -0.2 0.2 0.4]*vua;          % in m/sec

% From SNR, we compute relative RCS using the idea that SNR is proportional
% to RCS/R^4.  Students will be asked to deduce relative RCS.
rel_RCS = (10.^(SNR/10)).*(ranges.^4);
rel_RCS = db(rel_RCS/max(rel_RCS),'power')


% fprintf('\nThere are %g targets with the following parameters:',Ntargets)
% for i = 1:Ntargets
%   fprintf('\n  range=%5.2g km, SNR=%7.3g dB, rel_RCS=%7.3g dB, vel=%9.4g m/s', ...
%            ranges(i)/1e3,SNR(i),rel_RCS(i),vels(i) )
% end

% Now form the range bin - pulse number data map

% disp(' ')
% disp(' ')
% disp('... forming signal component')
y = radar(s,fs,T_0,g,T_out,T_ref,fc,ranges,SNR,vels);

% add thermal noise with unit power

% disp('... adding noise')
%randn('seed',77348911);
[My,Ny] = size(y);
nzz = (1/sqrt(2))*(randn(My,Ny) + J*randn(My,Ny));
y = y + nzz;

% create log-normal (ground) "clutter" with specified C/N and
% log-normal standard deviation for amplitude, uniform phase
% Clutter is uncorrelated in range, fully correlated in pulse #

disp('... creating clutter')
% CN = 20;         % clutter-to-noise ratio in first bin (dB)
SDxdB = 3;       % in dB (this is NOT the sigma of the complete clutter)
ncc=10 .^((SDxdB*randn(My,Ny))/10);
ncc = ncc.*exp( J*2*pi*rand(My,Ny) );

% Force the power spectrum shape to be Gaussian

disp('... correlating and adding clutter')
G  = exp(-(0:4)'.^2/1.0);
G = [G;zeros(Ny-2*length(G)+1,1);G(length(G):-1:2)];

for i=1:My
  ncc(i,:)=ifft(G'.*fft(ncc(i,:)));
end
 
% rescale clutter to have desired C/N ratio
pcc = var(ncc(:));
ncc = sqrt((10^(CN/10))/pcc).*ncc;
% 10*log10(var(ncc(:))/var(nzz(:)))  % check actual C/N

% Now weight the clutter power in range for assume R^2 (beam-limited) loss
cweight = T_out(1)*((T_out(1) + (0:My-1)'*(1/fs)).^(-1));
cweight = cweight*ones(1,Np);
ncc = ncc.*cweight;

end

