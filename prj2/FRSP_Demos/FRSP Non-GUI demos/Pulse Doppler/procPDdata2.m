%
% procPDdata
%
%  process pulse Doppler radar project data
%
%  Written by J. H. McClellan
%  Modified by M. A. Richards
%
clear, hold off
format compact
J = sqrt(-1);
close all

% choose color or black and white plots

%bw=input('Color plots? (y/n): ','s');
bw = 'y';

% Get root file name for reading results

%file=input('Enter root file name for data file: ','s');
file = 'test';

eval(['load ',file,'.mat'])

fprintf('\nPulse length = %g microseconds\n',T/1e-6)
fprintf('Chirp bandwidth = %g Mhz\n',W/1e6)
fprintf('Sampling rate = %g Msamples/sec\n',fs/1e6)
figure(1)
plot((1e6/fs)*(0:length(s)-1),[real(s) imag(s)])
title('Real and Imaginary Parts of Chirp Pulse')
xlabel('time (usec)')
ylabel('amplitude')
grid

PRI = 1/PRF;
fprintf('\nWe are simulating %g pulses at an RF of %g GHz',Np,fc/1e9)
fprintf('\nand a PRF of %g kHz, giving a PRI of %g usec.',PRF/1e3,PRI/1e-6)
fprintf('\nThe range window limits are %g to %g usec.\n', ...
    T_out(1)/1e-6,T_out(2)/1e-6)


% Compute unambiguous Doppler interval in m/sec
% Compute unambiguous range interval in meters

vua = 3e8*PRF/(2*fc);
rmin = 3e8*T_out(1)/2;
rmax = 3e8*T_out(2)/2;
rua = rmax-rmin;

fprintf('\nThe unambiguous velocity interval is %g m/s.',vua)
fprintf('\nThe range window starts at %g km.',rmin/1e3)
fprintf('\nThe range window ends at %g km.',rmax/1e3)
fprintf('\nThe unambiguous range interval is %g km.\n\n',rua/1e3)


% Convert range samples to absolute range units.

[My,Ny]=size(y);
range=(3e8/2)*((0:My-1)*(1/fs) + T_out(1))/1e3;

% Force oversize FFT, and compute doppler scale factor
Lfft = 32;
doppler = (((0:Lfft-1)/Lfft)-0.5)*vua;


% Plot overlay of individual range traces

disp(' ')
disp(' ')
disp('...plotting overlay of range traces')
figure(2)
if (bw ~= 'y')
  plot(range,db(y,'voltage'),'k')
else
  plot(range,db(y,'voltage'))
end
title('OVERLAY OF RANGE TRACES')
xlabel('distance (km)')
ylabel('amplitude (dB)')
grid


% Noncoherently integrate the range traces and display

disp('...plotting integrated range trace')
figure(3)
if (bw ~= 'y')
  plot(range,db(sum((abs(y).^2)')','power'),'k')
else
  plot(range,db(sum((abs(y).^2)')','power'))
end
title('NONCOHERENTLY INTEGRATED RANGE TRACE')
xlabel('range bin')
ylabel('power')
grid

% Doppler process and square-law detect the whole
% unprocessed array and display mesh.
% Use Hamming window throughout.

disp('...computing raw range-Doppler map')
Y=fft(conj(y').*(hamming(Ny)*ones(1,My)),Lfft);

% note we take mag-squared of Y here also
Y=git_rotate(Y.*conj(Y),Lfft/2);

YdB=db(abs(Y)/max(max(abs(Y))),'power');
figure(4)
if (bw ~= 'y')
   colormap([0 0 0]);
else
   colormap('default');
end
mesh(doppler,range,YdB')
title('RANGE-DOPPLER PLOT OF UNPROCESSED DATA')
ylabel('range (km)')
xlabel('velocity (m/s)')

levels=([-1 -5 -10 -15 -20 -25 -30]);
figure(5)
if (bw ~= 'y')
  contour(doppler,range,YdB',levels,'k')
else
  contour(doppler,range,YdB',levels)
end
title('RANGE-DOPPLER CONTOUR PLOT OF UNPROCESSED DATA')
ylabel('range (km)')
xlabel('velocity (m/s)')
grid



% % Apply three-pulse canceller in each range bin to raw data
% 
% disp('...performing 3-pulse clutter cancellation')
% yc = y * circulant([1;-2;1;zeros(Ny-3,1)],Ny-1);
% 
% % Doppler process and square-law detect the whole
% % clutter-cancelled array and display mesh
% 
% disp('...computing clutter-cancelled range-Doppler map')
% [Myc,Nyc]=size(yc);
% YC=fft(conj(yc').*(hamming(Nyc)*ones(1,Myc)),Lfft);
% % note we take mag-squared of YC here
% YC=git_rotate(YC.*conj(YC),Lfft/2);
% 
% YCdB=db(abs(YC)/max(max(abs(YC))),'power');
% figure(6)
% if (bw ~= 'y')
%    colormap([0 0 0]);
% else
%    colormap('default');
% end
% mesh(doppler,range,YCdB')
% title('RANGE-DOPPLER PLOT OF CLUTTER-CANCELLED DATA')
% ylabel('range (km)')
% xlabel('velocity (m/s)')
% 
% levels=([-1 -5 -10 -15 -20 -25]);
% figure(7)
% if (bw ~= 'y')
%   contour(doppler,range,YCdB',levels,'k')
% else
%   contour(doppler,range,YCdB',levels)
% end
% title('RANGE-DOPPLER CONTOUR PLOT OF CLUTTER-CANCELLED DATA')
% ylabel('range (km)')
% xlabel('velocity (m/s)')
% grid

yc = y;
% Perform matched filtering in range on
% clutter-cancelled data

Ls = length(s);
disp('...performing matched filtering')
h = conj( s(Ls:-1:1) );
ycf = fftfilt( h, yc );

[Mycf,Nycf]=size(ycf);
YCF=fft(conj(ycf').*(hamming(Nycf)*ones(1,Mycf)),Lfft);
% note that we take mag-squared of YCF here
YCF=git_rotate(YCF.*conj(YCF),Lfft/2);

% compute new range and time scales here to take filter delay of Ls-1
% samples into account.

ranged = range - (Ls-1)*(1/fs)*3e8/2/1e3;
timed = (1e3*ranged)*2/3e8;

% At this point I'm going to compensate for R^4 target power loss.
% Note that this will overcompensate for R^2 clutter power loss.
% Also note that YCF is already in power units.
cweight =(timed/T_out(1)).^4;
cweight = ones(Lfft,1)*cweight;
YCF = YCF.*cweight;

% Do the mesh and contour plots

YCFdB=db(abs(YCF)/max(max(abs(YCF))),'power');
figure(8)
if (bw ~= 'y')
   colormap([0 0 0]);
else
   colormap('default');
end
mesh(doppler,ranged,YCFdB')
title('RANGE-DOPPLER PLOT OF FULLY-PROCESSED DATA')
ylabel('range (km)')
xlabel('velocity (m/s)')

levels=([-1 -5 -10 -15 -20 -25 -30]);
figure(9)
if (bw ~= 'y')
  contour(doppler,ranged,YCFdB',levels,'k')
else
  contour(doppler,ranged,YCFdB',levels)
end
title('RANGE-DOPPLER CONTOUR PLOT OF FULLY-PROCESSED DATA')
ylabel('range (km)')
xlabel('velocity (m/s)')
grid

if (bw ~= 'y')
  colormap('default')
end

% ---- cut point
%  Find range bins for peaks.
%  Based on non-adaptive, manually set threshold
%  relative to highest peak

[x k] = max(YCFdB);
figure(10)
if (bw ~= 'y')
  plot(ranged,x,'k')
else
  plot(ranged,x)
end
xlabel('range (km)')
ylabel('power (dB)')
title('MAXIMUM DOPPLER RESPONSE VS. RANGE')
grid

y=max(x,-16); % this basically thresholds out all the noise and sidelobes
Ly=length(y);

% This loop identifies which range bins have local peaks;
% for each peak, it uses the index data from the 'max'
% function to get the peak locations in Doppler, and
% the surrounding two Doppler samples for each one.

YCFamp = sqrt(YCF);  % linear scale range-Doppler matrix
spikes = [];
for i=2:Ly-1
  if (y(i) > y(i+1)) && (y(i) > y(i-1))
    spikes=[spikes;[i,k(i),YCFamp(k(i)-1,i),YCFamp(k(i),i),YCFamp(k(i)+1,i)]];
  end
end
[Mspikes,Nspikes]=size(spikes);

% For each local peak found by the above, use 'peakinterp'
% to refine the amplitude and location estimate

target = [];
for i=1:Mspikes
  [amp del_k] = peakinterp(spikes(i,3:5));
  target = [target;[amp, ranged(spikes(i,1)), spikes(i,2)+del_k]];
end

% Put Doppler peaks into m/s units

target(:,3) = (((target(:,3)-1)/Lfft)-0.5)*vua;

% Adjust target amplitudes for MTI filter,
% and then renormalize.  Since the target 

target(:,1) = target(:,1)./(4*(sin(target(:,3)*(pi/vua))).^2);
target(:,1) = target(:,1)/max(target(:,1));
target(:,1) = db(target(:,1),'voltage');


% List out detected target amplitudes, ranges, Dopplers

fprintf('\n\nESTIMATED PARAMETERS OF DETECTED TARGETS:\n')
fprintf('\nNumber     Rel Amp (dB)     Range (km)     Vel (m/s)')

for i = 1:Mspikes
  fprintf('\n  %2.0g       %7.3g          %6.3g       %9.4g', ...
    i,target(i,:))
end
disp(' ')
disp(' ')

shg
