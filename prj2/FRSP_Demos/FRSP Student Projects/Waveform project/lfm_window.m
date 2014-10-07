% M file for experimenting with LFM range sidelobe control by windowing in both
% the time and frequency domains.
%
% Mark Richards
% March 2006

clear all, close all

% create a complex chirp with BT=100 and 10x oversampling for good
% definition in the plots.  It doesn't really matter so long as their
% product equals 100, but let's choose a pulse length of 1 microsecond and
% a swept bandwidth of 100 MHz to make it sound realistic.
oversample = 10;
tau = 1e-6;
BW = 100e6;
BT = BW*tau;

% The pulse length will be 1000 (BT product times the oversampling factor),
% so I'll use 4K FFTs for good definition in freq domain plots
Kfft = 4096;

% OK, create the basic LFM chirp and its spectrum
x = git_chirp(tau,BW,oversample);
X = fft(x,Kfft);
Xplot = abs(fftshift(X));

% set up time and freq abscissa plot vectors
t = [1:BT*oversample]'/BW/oversample-tau/2;  % in seconds; center it at t=0
f = [-Kfft/2:Kfft/2-1]'/Kfft*BW*oversample;  % in Hz; center it at f=0

figure(1)
subplot(1,3,1)
plot(t,real(x));
xlabel('Time (sec)'); ylabel('Re\{\itx\rm[\itn\rm]\}')
axis([-tau/2,tau/2,-1,1])
title('Real Part of LFM Waveform')

subplot(1,3,2)
plot(t,imag(x));
xlabel('Time (sec)'); ylabel('Re\{\itx\rm[\itn\rm]\}')
axis([-tau/2,tau/2,-1,1])
title('Imaginary Part of LFM Waveform')

subplot(1,3,3)
plot(f,Xplot);
xlabel('Frequency (Hz)')
ylabel('|\itX\rm(\itf\rm)|')
title('Magnitude of LFM Spectrum')

% OK, that was the easiest part.  Let's do 3 pulse compression cases now:
% no window, window applied in the time domain, and window
% applied in the frequency domain.

% First the standard matched filter with no sidelobe control.

hmf = conj(x(end:-1:1));  % matched filter impulse response

ymf = abs(conv(x,hmf));  % output of ideal matched filter;
                         % this will be twice as long as the individual signals
t2 = [1:length(ymf)]/BW/oversample-tau;

figure(2)
subplot(3,4,1)
plot(t2,ymf)
xlabel('Time (sec)'); ylabel('|\ity_{mf}\rm[\itn\rm]|')
grid
title('Magnitude of Matched Filter Ouput')

subplot(3,4,5)
plot(t2,db(ymf,'voltage'))
axis([-tau,tau,20*log10(BT*oversample)-60,20*log10(BT*oversample)])
xlabel('Time (sec)'); ylabel('|\ity_{mf}\rm[\itn\rm]|')
grid
title('Magnitude of Matched Filter Ouput (dB)')

subplot(3,4,9)
plot(t2,db(ymf/max(ymf),'voltage'))
axis([-tau,tau,-60,0])
xlabel('Time (sec)'); ylabel('|\ity_{mf}\rm[\itn\rm]|')
grid
title('Normalized Magnitude of Matched Filter Ouput (dB)')


% Now let's use a Hamming window, applied in the frequency domain, with
% its cutoffs at +/- BW/2. We need to use some caution here to line things
% up so the windowed spectrum has appropriate symmetry.

spread = 1;  % window will cover a bandwidth of spread*BW;
             % this allows me to play around with the window cutoff a little.
hmf = conj(x(end:-1:1));  % matched filter impulse response
Hmf = fft(hmf,Kfft);
kmax = floor(spread*Kfft/2/oversample);  % last DFT index corresponding to
                                         % window bandwidth or less
w = hamming(2*kmax+1);
wfd = [w(kmax+1:2*kmax+1);zeros(Kfft-length(w),1);w(1:kmax)];  % build frequency-domain window
Hfd = Hmf.*wfd;  % apply the window
hfd = ifft(Hfd,Kfft);
hfd = hfd(1:length(x));  % enforce a finite length on hfd

figure(3)  % check plot of window and spectrum alignment
plot(f,abs(fftshift([Hmf, max(abs(Hmf))*wfd, Hfd])));
xlabel('Frequency (Hz)')
ylabel('Spectral Magnitude')
title('Magnitude of Window, Before & After Spectra, Freq. Domain Windowing')

figure(4)  % check plot of mis-matched filter impulse response
subplot(2,2,1)
plot(real(hfd)); title('Real Part of Impulse Response')
subplot(2,2,2)
plot(imag(hfd)); title('Imaginary Part of Impulse Response')
subplot(2,2,3)
plot(abs(hfd)); title('Magnitude of Impulse Response')
subplot(2,2,4)
plot(unwrap(angle(hfd))); title('Unwrapped Phase of Impulse Response')
suptitle('Frequency Domain Windowing')

yfd = abs(conv(x,hfd));  % output of mis-matched filter;
                         % this will be twice as long as the individual signals

figure(2)
subplot(3,4,2)
plot(t2,yfd)
xlabel('Time (sec)'); ylabel('|\ity_{fd}\rm[\itn\rm]|')
grid
title('Magnitude of Matched Filter Ouput')

subplot(3,4,6)
plot(t2,db(yfd,'voltage'))
axis([-tau,tau,20*log10(BT*oversample)-60,20*log10(BT*oversample)])
xlabel('Time (sec)'); ylabel('|\ity_{fd}\rm[\itn\rm]|')
grid
title('Magnitude of Matched Filter Ouput (dB)')

subplot(3,4,10)
plot(t2,db(yfd/max(yfd),'voltage'))
axis([-tau,tau,-60,0])
xlabel('Time (sec)'); ylabel('|\ity_{fd}\rm[\itn\rm]|')
grid
title('Normalized Magnitude of Matched Filter Ouput (dB)')


% Now let's use a Hamming window, applied in the **time** domain.

hmf = conj(x(end:-1:1));  % matched filter impulse response
w = hamming(length(hmf));
hft = w.*hmf;
Hft = fft(hft, Kfft);

figure(5)  % informational plot of matched and mis-matched filter spectra
plot(f,abs(fftshift([Hmf, Hft])));
xlabel('Frequency (Hz)')
ylabel('Spectral Magnitude')
title('Magnitude of Before & After Spectra, Time-Domain Windowing')

figure(6)  % check plot of mis-matched filter impulse response
subplot(2,2,1)
plot(real(hft)); title('Real Part of Impulse Response')
subplot(2,2,2)
plot(imag(hft)); title('Imaginary Part of Impulse Response')
subplot(2,2,3)
plot(abs(hft)); title('Magnitude of Impulse Response')
subplot(2,2,4)
plot(unwrap(angle(hft))); title('Unwrapped Phase of Impulse Response')
suptitle('Time Domain Windowing')

yft = abs(conv(x,hft));  % output of mismatched filter;
                         % this will be twice as long as the individual signals

figure(2)
subplot(3,4,3)
plot(t2,yft)
xlabel('Time (sec)'); ylabel('|\ity_{ft}\rm[\itn\rm]|')
grid
title('Magnitude of Matched Filter Ouput')

subplot(3,4,7)
plot(t2,db(yft,'voltage'))
axis([-tau,tau,20*log10(BT*oversample)-60,20*log10(BT*oversample)])
xlabel('Time (sec)'); ylabel('|\ity_{ft}\rm[\itn\rm]|')
grid
title('Magnitude of Matched Filter Ouput (dB)')

subplot(3,4,11)
plot(t2,db(yft/max(yft),'voltage'))
axis([-tau,tau,-60,0])
xlabel('Time (sec)'); ylabel('|\ity_{ft}\rm[\itn\rm]|')
grid
title('Normalized Magnitude of Matched Filter Ouput (dB)')

% One more variation: amplitude modulate the LFM waveform itself with the
% square root of the window in time, then use matched filter for that case.
% I will increase the amplitude of the AM'[ed LFM to have the same energy
% as all of the other pulses.

xsq = x.*sqrt(hamming(length(x)));
xsq = xsq*norm(x)/norm(xsq);
hsq = conj(xsq(end:-1:1));  % matched filter impulse response
Hsq = fft(hsq, Kfft);

figure(7)  % informational plot of filter spectrum
plot(f,abs(fftshift(Hsq)));
xlabel('Frequency (Hz)')
ylabel('Spectral Magnitude')
title('Magnitude of Spectra, Amplitude Modulation')

figure(8)  % check plot of filter impulse response
subplot(2,2,1)
plot(real(hsq)); title('Real Part of Impulse Response')
subplot(2,2,2)
plot(imag(hsq)); title('Imaginary Part of Impulse Response')
subplot(2,2,3)
plot(abs(hsq)); title('Magnitude of Impulse Response')
subplot(2,2,4)
plot(unwrap(angle(hsq))); title('Unwrapped Phase of Impulse Response')
suptitle('Time Domain Amplitude Modulation')

ysq = abs(conv(xsq,hsq));  % output of ideal matched filter;
                           % this will be twice as long as the individual signals

figure(2)
subplot(3,4,4)
plot(t2,ysq)
xlabel('Time (sec)'); ylabel('|\ity_{sq}\rm[\itn\rm]|')
grid
title('Magnitude of Matched Filter Ouput')

subplot(3,4,8)
plot(t2,db(ysq,'voltage'))
axis([-tau,tau,20*log10(BT*oversample)-60,20*log10(BT*oversample)])
xlabel('Time (sec)'); ylabel('|\ity_{sq}\rm[\itn\rm]|')
grid
title('Magnitude of Matched Filter Ouput (dB)')

subplot(3,4,12)
plot(t2,db(ysq/max(ysq),'voltage'))
axis([-tau,tau,-60,0])
xlabel('Time (sec)'); ylabel('|\ity_{sq}\rm[\itn\rm]|')
grid
title('Normalized Magnitude of Matched Filter Ouput (dB)')

suptitle('Matched, No AM           Mismatched, Freq          Mismatched, Time              AM')