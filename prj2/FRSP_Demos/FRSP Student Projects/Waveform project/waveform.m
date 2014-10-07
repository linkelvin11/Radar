% waveform
%
% M-file for ECE6272 computer project #2 on waveforms
%
% Mark A. Richards, Feb. 2000
% updated Feb. 2002
% updated Sep. 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up pulse length and sampling rate
T = input('Enter pulse length (sec): ');  % 100 us pulse length

% use supplied chirp function to create complex chirp of
% specified bandwidth and length
B = input('Enter chirp swept bandwidth (Hz): ');
OS = input('Enter chirp oversampling factor: ');
disp(['Chirp BT product =',num2str(B*T)]);
Fs = OS*B;  % fast time sampling rate
Ts = 1/Fs;  % fast time sampling interval
xc = git_chirp(T,B,OS);
N = length(xc);

% create simple pulse of same duration and sampling rate
xs=ones(N,1);

% do matched filter output (autocorrelation function) of each
% waveform and compare on same scale
xcc = xcorr(xc);
xss = xcorr(xs);
time = (-N+1:+N-1)*Ts;
figure(1);
plot(time,abs([xss xcc])); xlabel('time (sec)'); ylabel('amplitude');

% normalize and repeat plot on log scale
xcc = xcc/max(abs(xcc));
xss = xss/max(abs(xss));
figure(2);
plot(time,20*log10(abs([xss xcc])));
axis([min(time) max(time) -50 0]);
xlabel('time (sec)'); ylabel('amplitude (dB)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now we want to look at range resolution.  Set up two scatterers
% 10 us apart; the received signals are then duplicates of the
% transmitted signal, overlapped with appropriate delay.
offset = round(10e-6/Ts);
rc = [xc;zeros(offset,1)] + [zeros(offset,1);xc];
rs = [xs;zeros(offset,1)] + [zeros(offset,1);xs];

% now do the match filtering by correlating with the respective
% reference pulse

rcc = xcorr(rc,xc);
rss = xcorr(rs,xs);
time = (-N-offset+1:+N+offset-1)*Ts;
figure(3);
plot(time,abs([rss rcc])); xlabel('time (sec)'); ylabel('amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now we want to do sidelobe suppression on the chirp, so we switch to 
% frequency domain matched filtering.

% compute the spectrum of the chirp signal and fftshift it to put
% the origin in the middle of the array (this shift is mainly
% for plotting convenience, it is not necessary).  We use a 16384
% point fft to get good spectrum definition.
Nfft = 4*4096;
Nf2 = Nfft/2;
XC = fftshift(fft(xc,Nfft));
freq = ((0:Nfft-1)-Nf2)/Nfft;

% compute the length of the spectral segment that contains
% most of the spectrum energy, namely the number of samples
% that cover a range of B Hz.  Keep it integer.  Generate the
% Hamming window and get it centered in an array the same size
% as the signal spectrum.
nspec = round(Nfft/OS);
h = hamming(nspec);
hpad=[zeros(Nf2-1-floor(nspec/2),1);h;zeros(Nf2+1-ceil(nspec/2),1)];
figure(4)
plot(freq,[abs(XC)/max(abs(XC)) hpad]);
xlabel('normalized frequency (cycles)'); ylabel('amplitude');

% apply the window to the matched filter frequency response
HC = hpad.*conj(XC);

% Do the matched filter with and without weighting,
% inverse transform, and plot the result.
% Compare to the unwindowed time-domain matched filter.
% Plot on dB scale to facilitate looking at sidelobes.
% Also pad out and align time-domain unweighted matched filter
% for comparison plotting.  Plot only for +/- T seconds.
YC = HC.*XC;
yc = fftshift(ifft(fftshift(YC)));
ZC = conj(XC).*XC;
zc = fftshift(ifft(fftshift(ZC)));
time = Ts*(-Nf2:Nf2-1);
indexT = (abs(time)<=T);
figure(5);
plot(time(indexT),20*log10(abs([yc(indexT) zc(indexT)])));
big = max(20*log10(abs(zc))); axis([-T +T big-60 big]);
xlabel('time (sec)'); ylabel('amplitude (dB)')

% check loss in processing gain
LPG_measured = 20*log10(max(abs(yc))) - 20*log10(max(abs(zc)))
LPG_formula = 10*log10((sum(h)/length(h))^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Repeat the frequency domain sidelobe suppression on the chirp with a
% 10% increase in span of the window

% compute the spectrum of the chirp signal and fftshift it to put
% the origin in the middle of the array (this shift is mainly
% for plotting convenience, it is not necessary).  We use a 16384
% point fft to get good spectrum definition.
Nfft = 4*4096;
Nf2 = Nfft/2;
XC = fftshift(fft(xc,Nfft));
freq = ((0:Nfft-1)-Nf2)/Nfft;

% compute the length of the spectral segment that contains
% most of the spectrum energy, namely the number of samples
% that cover a range of B Hz.  Keep it integer.  Generate the
% Hamming window and get it centered in an array the same size
% as the signal spectrum.
nspec = round(1.1*Nfft/OS);
h = hamming(nspec);
hpad=[zeros(Nf2-1-floor(nspec/2),1);h;zeros(Nf2+1-ceil(nspec/2),1)];
figure(6)
plot(freq,[abs(XC)/max(abs(XC)) hpad]);
xlabel('normalized frequency (cycles)'); ylabel('amplitude');

% apply the window to the matched filter frequency response
HC = hpad.*conj(XC);

% Do the matched filter with and without weighting,
% inverse transform, and plot the result.
% Compare to the unwindowed time-domain matched filter.
% Plot on dB scale to facilitate looking at sidelobes.
% Also pad out and align time-domain unweighted matched filter
% for comparison plotting.  Plot only for +/- T seconds.
YCP = HC.*XC;
ycp = fftshift(ifft(fftshift(YCP)));
ZC = conj(XC).*XC;
zc = fftshift(ifft(fftshift(ZC)));
time = Ts*(-Nf2:Nf2-1);
indexT = (abs(time)<=T);
figure(7);
plot(time(indexT),20*log10(abs([ycp(indexT) zc(indexT)])));
big = max(20*log10(abs(zc))); axis([-T +T big-60 big]);
xlabel('time (sec)'); ylabel('amplitude (dB)')

% check loss in processing gain
LPG_measured = 20*log10(max(abs(ycp))) - 20*log10(max(abs(zc)))
LPG_formula = 10*log10((sum(h)/length(h))^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now try time domain sidelobe suppression on the chirp.  We go back to
% time domain matched filtering, since that is more convenient now

% do matched filter output (autocorrelation function) of waveform with and
% without hamming weighting
xcc = xcorr(xc);
xcw = xcorr(xc,hamming(N).*xc);
time = (-N+1:+N-1)*Ts;
figure(8);
plot(time,abs([xcw, xcc]));
xlabel('time (sec)'); ylabel('amplitude');

figure(9);
plot(time,20*log10(abs([xcw, xcc])));
big = max(20*log10(abs(xcc))); axis([-T +T big-60 big]);
xlabel('time (sec)'); ylabel('amplitude (dB)')

% check loss in processing gain
LPG_measured = 20*log10(max(abs(xcc))) - 20*log10(max(abs(xcw)))
LPG_formula = 10*log10((sum(h)/length(h))^2)

% Do a plot comparing frequency- and time-domain weighted responses
figure(10);
yc = yc(indexT); yc = yc(2:end-1);
plot(time,20*log10(abs([yc xcw])));
big = max(20*log10(abs(zc))); axis([-T +T big-60 big]);
xlabel('time (sec)'); ylabel('amplitude (dB)')

% Do another plot to follow sidelobe level trends of the two
[pxcw txcw] = peaks(20*log10(abs(xcw)));
[pyc tyc] = peaks(20*log10(abs(yc)));
figure(11)
plot(time(pyc),20*log10(yc(pyc))); hold on
plot(time(pxcw),20*log10(xcw(pxcw)),'g'); hold off
big = max(20*log10(abs(zc))); axis([-T +T big-60 big]);
xlabel('time (sec)'); ylabel('amplitude (dB)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%range-Doppler coupling
%
% first create a Doppler shift term corresponding to an expected time
% shift due to range-Doppler coupling of 10 us
doppler=exp(j*2*pi/10/OS*(0:N-1)');

% now create doppler-shifted response to our original up chirp and process
% through the matched filter for the original waveform
xc1 = doppler.*xc;
xcc1 = xcorr(xc1,xc);

% to create a down chirp, all we need to do is conjugate the original
% up chirp; so repeat above using conj(xc)
xc2=doppler.*conj(xc);
xcc2 = xcorr(xc2,conj(xc));

% plot the two responses on a common axis
time = (-N+1:+N-1)*Ts;
figure(12);
plot(time,abs([xcc1 xcc2])); xlabel('time (sec)'); ylabel('amplitude');


