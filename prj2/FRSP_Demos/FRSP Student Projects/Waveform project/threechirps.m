% threechirps
%
% M-file for ECE6272 computer project #2 on waveforms
%
% Mark A. Richards, Feb. 2000


% set up pulse length and oversampling rate
T = 100e-6;  % pulse length (sec)
OS = 1.2;  %  chirp oversampling factor

% use supplied chirp function to create three complex chirps of
% specified bandwidth and length and oversample ratio

% note that each signal will be of a different length because the
% bandwidths are not the same, therefore neither are the sampling rates
xc10 = git_chirp(T,10/T,OS);      % BT=10
xc100 = git_chirp(T,100/T,OS);    % BT=100
xc1000 = git_chirp(T,1000/T,OS);  % BT=1000

% power in each signal also is proportional to bandwidth, because
% sampling rate is proportional to bandwidth but time duration
% is constant.  For  "niceness" of the plot, normalize by square 
% root to length so all spectra will have about the same amplitude.
xc10 = xc10/sqrt(length(xc10));
xc100 = xc100/sqrt(length(xc100));
xc1000 = xc1000/sqrt(length(xc1000));

% compute the spectra of all three, shifting the origin to the
% center of the plot, and plotting against normalized frequency

N=length(xc1000);
XC10 = abs(fftshift(fft(xc10,N)));
XC100 = abs(fftshift(fft(xc100,N)));
XC1000 = abs(fftshift(fft(xc1000,N)));

% define frequency variable and do the plot
freq = ((0:N-1)/N)-0.5;
plot(freq,[XC10 XC100 XC1000]);
xlabel('normalized frequency (cycles)'); ylabel('spectrum amplitude')