function pdatanpc()
%
% procdata
%
%  process pulse Doppler radar project data
%
%  Written by J. H. McClellan
%  Modified by M. A. Richards
%
% Updated by M. A. Richards, Oct. 2006
%

% clear, 
hold off
format compact
J = sqrt(-1);
close all

% Get root file name for reading results

% file=input('Enter root file name for data file: ','s');
file = 'tmp';

eval(['load ',file,'.mat'])

fprintf('\nPulse length = %g microseconds\n',T/1e-6)
fprintf('Chirp bandwidth = %g Mhz\n',W/1e6)
fprintf('Sampling rate = %g Msamples/sec\n',fs/1e6)
figure
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
rwin = rmax-rmin;
rua = 3e8/2/PRF;

fprintf('\nThe unambiguous velocity interval is %g m/s.',vua)
fprintf('\nThe range window starts at %g km.',rmin/1e3)
fprintf('\nThe range window ends at %g km.',rmax/1e3)
fprintf('\nThe sampled range window is %g km long.',rwin/1e3)
fprintf('\nThe unambiguous range interval is %g km.\n\n',rua/1e3)


% Convert range samples to absolute range units.
[My,Ny]=size(y);
range=(3e8/2)*((0:My-1)*(1/fs) + T_out(1))/1e3;

pulse = (1:Ny);

% Force oversize FFT, and compute doppler scale factor
Lfft = 2^(nextpow2(Ny)+3);
doppler = (((0:Lfft-1)/Lfft)-0.5)*vua;

fprintf('\nThe Doppler increment is %g Hz.',PRF/Lfft)
fprintf('\nThe velocity increment is %g m/s.',3e8*PRF/Lfft/2/fc)


% Start with a few plots to examine the data


% plot power of raw data in dB
ydB=db(abs(y)/max(max(abs(y))),'voltage');
figure
mesh(pulse,range,ydB)
title('FAST-TIME/SLOW-TIME PLOT OF RAW DATA')
ylabel('range (km)')
xlabel('pulse number')

% Plot overlay of individual range traces
disp(' ')
disp(' ')
disp('...plotting overlay of range traces')
figure
plot(range,db(y,'voltage'))
title('OVERLAY OF RANGE TRACES')
xlabel('distance (km)')
ylabel('amplitude (dB)')
grid

% Noncoherently integrate the range traces and display
disp('...plotting integrated range trace')
figure
plot(range,db(sum((abs(y).^2)')','power'))
title('NONCOHERENTLY INTEGRATED RANGE TRACE')
xlabel('range bin')
ylabel('power')
grid

% Doppler process and square-law detect the whole
% unprocessed array and display mesh.
% Use Hamming window throughout.

disp('...computing raw range-Doppler map')
Y=fft(conj(y').*(hamming(Ny)*ones(1,My)),Lfft);
Y=git_rotate(Y.*conj(Y),Lfft/2);  % note we take mag-squared of Y here also
YdB=db(abs(Y),'power');
figure
mesh(doppler,range,YdB')
title('RANGE-DOPPLER PLOT OF UNPROCESSED DATA')
ylabel('range (km)')
xlabel('velocity (m/s)')

levels=(max(YdB(:))+[-1 -5 -10 -15 -20 -25 -30]);
figure
contour(doppler,range,YdB',levels)
title('RANGE-DOPPLER CONTOUR PLOT OF UNPROCESSED DATA')
ylabel('range (km)')
xlabel('velocity (m/s)')
grid

% Now start processing the data ...

% Pulse compression first.  Use time-domain Hamming weighting of the
% impulse response for range sidelobe control
Ls = length(s);
disp('...performing matched filtering')
h = conj( s(Ls:-1:1) );
h = h.*hamming(length(h));  % time-domain Hamming window for range sidelobe control
yp = zeros(My+length(h)-1,Ny);
for i=1:Ny
    yp(:,i) = conv(h,y(:,i));
end
[Myp,Nyp]=size(yp);
% yp = fftfilt( h, y ); % using fftfilt instead of conv because it filters
%                       % multiple columns with one call

% compute new range and time scales here to take account of increased range
% length due to convolution and offset due to filter delay of Ls-1 samples.
rangep = (3e8/2)*(((0:Myp-1)-(Ls-1))*(1/fs) + T_out(1))/1e3;
timep = (1e3*rangep)*2/3e8;


ypdB=db(abs(yp),'voltage');
figure
mesh(pulse,rangep,ypdB)
title('FAST-TIME/SLOW-TIME PLOT OF PULSE-COMPRESSED DATA')
ylabel('range (km)')
xlabel('pulse number')


levels=(max(ypdB(:))+[-1 -5 -10 -15 -20]);
figure
contour(pulse,rangep,ypdB,levels)
title('FAST-TIME/SLOW-TIME CONTOUR PLOT OF PULSE-COMPRESSED DATA')
ylabel('range (km)')
xlabel('pulse number')
grid

% Range-Doppler plots of pulse-compressed data
YP=fft(conj(yp').*(hamming(Nyp)*ones(1,Myp)),Lfft);
YP=git_rotate(YP.*conj(YP),Lfft/2);
YPdB=db(abs(YP),'power');
figure
mesh(doppler,rangep,YPdB')
title('RANGE-DOPPLER PLOT OF PULSE-COMPRESSED DATA')
ylabel('range (km)')
xlabel('velocity (m/s)')


levels=(max(YPdB(:))+[-1 -5 -10 -15 -20 -25 -30]);
figure
contour(doppler,rangep,YPdB',levels)
title('RANGE-DOPPLER CONTOUR PLOT OF PULSE-COMPRESSED DATA')
ylabel('range (km)')
xlabel('velocity (m/s)')
grid

% Apply three-pulse canceller in each range bin to raw data

disp('...performing 3-pulse clutter cancellation')
% h = [1 -2 1]';
h = 1;
ypm = zeros(Myp,Nyp+length(h)-1);
for i=1:Myp
    ypm(i,:) = conv(h,yp(i,:));
end
[Mypm,Nypm]=size(ypm);

% Doppler process and square-law detect the whole
% clutter-cancelled array and display mesh

disp('...computing clutter-cancelled range-Doppler map')
YPM=fft(conj(ypm').*(hamming(Nypm)*ones(1,Mypm)),Lfft);
YPM=git_rotate(YPM.*conj(YPM),Lfft/2);

YPMdB=db(abs(YPM),'power');
% figure
% mesh(doppler,rangep,YPMdB')
% title('RANGE-DOPPLER PLOT OF CLUTTER-CANCELLED DATA')
% ylabel('range (km)')
% xlabel('velocity (m/s)')

levels=(max(YPMdB(:))+[-1 -5 -10 -15 -20 -25]);
% figure
% contour(doppler,rangep,YPMdB',levels)
% title('RANGE-DOPPLER CONTOUR PLOT OF CLUTTER-CANCELLED DATA')
% ylabel('range (km)')
% xlabel('velocity (m/s)')
% grid


% OK, now let's do a 1D range-only CFAR to set a threshold map.
% Define CFAR window 23 in range, including test cell and 2 guard cells, so
% 18 actual averaging cells
cfar = ones(23,1)/18;
cfar(10:14)=0;
% Now convolve it in range dimension with range-Doppler map YPM
% to get average noise map; discard half-transients to maintain data array size
N = zeros(size(YPM));
for i=1:Lfft
    temp = conv(YPM(i,:),cfar);
    N(i,:) = temp(12:end-11);
end

% plot the resulting noise map
NdB=db(abs(N),'power');
figure
mesh(doppler,rangep,NdB')
title('RANGE-DOPPLER PLOT OF CFAR NOISE ESTIMATE')
ylabel('range (km)')
xlabel('velocity (m/s)')

% multiply the noise by a factor of 16 (20 dB) and then threshold the data
threshold = 16*N;
detected = NaN*ones(size(YPM));
detected(YPM>threshold) = YPM(YPM>threshold);

% plot the resulting detection map on a dB scale
figure
mesh(doppler,rangep,db(detected','power'))
title('RANGE-DOPPLER DETECTION MAP')
ylabel('range (km)')
xlabel('velocity (m/s)')

% To search for range bins with targets, first noncoherently integrate
% across the frequency bins
YPMrange = sum(YPM);
% noise floor estimate based on median to avoid elevation of estimate by
% target responses
Nrange = median(YPMrange);
Trange = 8*Nrange;  % threshold 8x (9 DB) above noise estimate

% This loop identifies which range bins have local peaks above the
% threshold. It also sets up a vector for plotting convenience to circle
% the peaks that are found.
spikesr = [];
marker_r = NaN*ones(1,length(YPMrange));
for i=2:Mypm-1
    if ((YPMrange(i) > YPMrange(i+1)) & (YPMrange(i) > YPMrange(i-1)) & (YPMrange(i) > Trange))
        spikesr=[spikesr;i,YPMrange(i)];
        marker_r(i) = YPMrange(i);
    end
end
[Mspikesr,Nspikesr]=size(spikesr);

figure
plot(rangep,db([YPMrange;Nrange*ones(1,length(YPMrange)); ...
    Trange*ones(1,length(YPMrange))],'power'));
hold on
plot(rangep,db(marker_r,'power'),'-ro')
hold off
xlabel('range (km)')
ylabel('power (dB)')
title('RANGE PEAKS')
grid

targets = [];
% Now find the Doppler peak(s) for each range bin having a target(s).  Keep
% adjoining Doppler values as well to support subsequent interpolation.
for i = 1:Mspikesr
    rb = spikesr(i,1)  % current range bin
    % search in Doppler in this range bin.  Have to do this modulo the
    % DFT size to allow for Doppler peaks near the spectrum edges.
    spikesd = [];
    marker_d = NaN*ones(1,Lfft);
    for k=1:Lfft
        km1 = k-1;
        if km1 < 1
            km1 = km1+Lfft;
        end
        kp1 = k+1;
        if kp1 > Lfft
            kp1 = kp1-Lfft;
        end
        if ((YPM(k,rb) > YPM(kp1,rb)) & (YPM(k,rb) > YPM(km1,rb)) & (YPM(k,rb) > Nrange/2))
            spikesd=[spikesd;k,YPM(km1,rb),YPM(k,rb),YPM(kp1,rb)];
            marker_d(k) = YPM(k,rb);
        end
    end  % end of loop over Doppler bins
    [Mspikesd,Nspikesd]=size(spikesd);

    figure
    plot(doppler,db([YPM(:,rb),Nrange*ones(Lfft,1)],'power'));
    hold on
    plot(doppler,db(marker_d,'power'),'-ro')
    hold off
    xlabel('velocity (m/s)')
    ylabel('power (dB)')
    title(['DOPPLER PEAKS FOR RANGE BIN ',int2str(rb)])
    grid

    % For each local peak in Doppler, use 'peakinterp'
    % to refine the amplitude and location estimate using a quadratic
    % interpolation

    for i=1:Mspikesd
        % note that peakinterp works on magnitude, not magnitude-squared
        [amp del_k] = peakinterp(sqrt(spikesd(i,2:4)));
        targets = [targets;[amp, rangep(rb), spikesd(i,1)+del_k]];
    end
end  % end of loop over range bins containing peaks

[Mtargets, Ntargets] = size(targets);

fprintf('\n\nINTERIM PARAMETERS OF DETECTED TARGETS:\n')
fprintf('\nNumber       Power         Range (km)       DFT Index')

for i = 1:Mtargets
    fprintf('\n  %2.0g       %7.3g          %6.3g       %9.4g', ...
        i,targets(i,:))
end
disp(' ')
disp(' ')

% OK, now 'targets' contains a list, we hope, of all targets with the
% range bin, interpolated peak magnitude, and interpolated Doppler bin of
% each.  Now we need to get in the desired units and adjust the amplitudes
% for a couple of factors.

% Put Doppler peaks into m/s units

targets(:,3) = (((targets(:,3)-1)/Lfft)-0.5)*vua;

% Adjust target amplitudes for MTI filter,
% and then convert to relative RCS

targets(:,1) = targets(:,1)./(4*(sin(targets(:,3)*(pi/vua))).^2);


fprintf('\n\nINTERIM PARAMETERS OF DETECTED TARGETS AFTER MTI FILTER CORRECTION:\n')
fprintf('\nNumber       Power         Range (km)       Vel (m/s)')

for i = 1:Mtargets
    fprintf('\n  %2.0g       %7.3g          %6.3g       %9.4g', ...
        i,targets(i,:))
end
disp(' ')
disp(' ')


targets(:,1) = targets(:,1).*targets(:,2).^2;
targets(:,1) = targets(:,1)/max(targets(:,1));
targets(:,1) = db(targets(:,1),'voltage');

% List out detected target amplitudes, ranges, Dopplers

fprintf('\n\nESTIMATED PARAMETERS OF DETECTED TARGETS:\n')
fprintf('\nNumber     Rel RCS (dB)     Range (km)     Vel (m/s)')

for i = 1:Mtargets
    fprintf('\n  %2.0g       %7.3g          %6.3g       %9.4g', ...
        i,targets(i,:))
end
disp(' ')
disp(' ')
