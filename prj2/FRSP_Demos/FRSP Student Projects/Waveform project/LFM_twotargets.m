% LFM two_targets
%
% M-file for ECE6272 computer project #2 on waveforms
%
% Mark A. Richards, Sept. 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % set up pulse length and sampling rate
% T = input('Enter pulse length (sec): ');  % 100 us pulse length
% 
% % use supplied chirp function to create complex chirp of
% % specified bandwidth and length
% B = input('Enter chirp swept bandwidth (Hz): ');
% OS = input('Enter chirp oversampling factor: ');

T = 100e-6; B = 1e6; OS = 2; c = 3e8;

disp(['Chirp BT product =',num2str(B*T)]);
Fs = OS*B;  % fast time sampling rate
Ts = 1/Fs;  % fast time sampling interval
xc = git_chirp(T,B,OS);
N = length(xc);

T1 = 101;  T2 =141;

P2 = 180;  % phase shift of 2nd target in degrees
P2 = (pi/180)*P2;

% create simple pulse of same duration and sampling rate
xs=ones(N,1);

% create two-target echo data for both types of pulse
ys = zeros(401,1);
ys(T1:T1+N-1) = xs;
ys(T2:T2+N-1) = ys(T2:T2+N-1) + exp(j*P2)*xs;

yc = zeros(401,1);
yc(T1:T1+N-1) = xc;
yc(T2:T2+N-1) = yc(T2:T2+N-1) + exp(j*P2)*xc;

% do matched filter output of each waveform by explicit convolution.  In
% the LFM case, include Hamming window for range sideobe suppression
hs = xs;
zs = conv(ys,hs);

hc = conj(xc(end:-1:1)).*hamming(N);
zc = conv(yc,hc);

range = (0:400+N-1)*(c*Ts/2) + 20e3 - (199)*(c*Ts/2);
figure(1);
plot(range/1e3,abs([zs zc])); xlabel('Range (km)'); ylabel('amplitude');
grid

