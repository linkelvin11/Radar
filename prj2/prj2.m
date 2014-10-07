%% Kelvin Lin Radar Project 2
% This matlab file serves as a test bench, and was only used for ease of
% data collection.
clc; clear all; close all;
%% Graph Changes In Clutter vs CNR
close all;
% mclutter generates clutter (cut from makedata)

figure
n = 3;
N = n*n;
d = 5;
for ii = d:d:d*N;
    y = mclutter(ii);
    y = abs(sort(y));
    len = length(y);
    x1 = y(floor(len/4));
    x2 = y(floor(len/2));
    x3 = y(floor(3*len/4));
    subplot(n,n,ii/d)
    hist(y,len); hold on;
    M = 15;
    plot(line([x3 x3],[0 M],'LineWidth',1,'Color',[0 0 0]));
    plot(line([x2 x2],[0 M],'LineWidth',1,'Color',[0 0 0]));
    plot(line([x1 x1],[0 M],'LineWidth',1,'Color',[0 0 0]));
    title(strcat('CNR = ',num2str(ii)));
    ylim([0 M]);
    xlim([min(min(y)) max(max(y))]);
end
suptitle('Clutter Amplitude vs CNR');

%% Test Effectiveness of MTI - OFF
close all;

% There are 4 targets with the following parameters:
%   range=    2 km, SNR=     -3 dB, rel_RCS=  -26.7 dB, vel=      -60 m/s
%   range=  3.8 km, SNR=      5 dB, rel_RCS=  -7.55 dB, vel=      -30 m/s
%   range=  4.4 km, SNR=     10 dB, rel_RCS=      0 dB, vel=       30 m/s
%   range=  4.4 km, SNR=      7 dB, rel_RCS=     -3 dB, vel=       60 m/s

mdatac(10);
pdatanpc;
close all;
% At CNR = 10 is when a non-MTI system seems to break down. Even at CNR =
% 9, pdatanpc was still able to resolve only 4 targets for a sizeable
% majority of the time. However, at CNR = 10, there were many false alarms.

% ESTIMATED PARAMETERS OF DETECTED TARGETS: @CNR = 10
% 
% Number     Rel RCS (dB)     Range (km)     Vel (m/s)
%    1             0            1.24         0.01448
%    2         -46.7             1.3         -0.2123
%    3           -83            1.36           1.827
%    4         -21.5            1.66         0.06558
%    5          -130               2          -59.97
%    6          -103             3.8             -30
%    7           -95             4.4           29.87
%    8          -106             4.4            60.1 
 
%% Test Effectiveness of MTI - ON
close all;

mdatac(50);
pdata;

% The pulse canceller is extremely effective for high values of CNR. Even
% at CNR = 45, the algorithm is still able to resolve all four targets.

% At CNR = 49, the processing will occasionally drop the target at 3.8 km
% -30m/s.

% ESTIMATED PARAMETERS OF DETECTED TARGETS: @CNR = 49
% 
% Number     Rel RCS (dB)     Range (km)     Vel (m/s)
%    1         -25.1               2          -62.31
%    2             0             4.4           30.45
%    3         -2.76             4.4           60.39 

% At CNR = 50, we see that two targets are dropped

% ESTIMATED PARAMETERS OF DETECTED TARGETS:
% 
% Number     Rel RCS (dB)     Range (km)     Vel (m/s)
%    1             0             4.4           30.88
%    2         -2.68             4.4           60.86 

% Take note that while at CNR49 the slow moving targets were the first to
% suffer, it was the target at the lower range which suffered first.
% However, once CNR was increased, the target that suffered next was
% actually the target at 2km. This suggests that targets at lower range
% suffer more than targets at low velocity.

% looking at the clutter generating code, it becomes obvious why. The
% clutter is correlated with the range due to the radar equation

%% Confirm Previously stated hypothesis
close all;

cn = 55
mdatacr(cn);        %generates a grid of targets at SNR = 10 and the given CNR
pdatacr;

%data gathered stored in txt file

% we see that as the CNR increases, the lower targets are the first to
% suffer. This is due to the fact that the stopband of the
% 3-pulse-canceller actually widens with increasing CNR. Additionally, the
% clutter is observed to be widening in doppler as the CNR increases. The
% furthest targets are the last ones to be affected due to the fact that
% the clutter is correlated with the range. As range increases, clutter
% power decreases. This is to be expected, as the radar equation includes a
% 1/R^4 term.

%% Move targets slower
close all;
%input to mdataslow(v) generates a lowest absolute velocity of v*15m/s
mdataslow(1);
pdatacr;

% slower targets are still affected first, and closer targets
% suffer before further ones. 

%% change bandwidth
close all;
tic
mdatabw(60,1e8);
pdatabw;
toc

% an increase in chirp bandwidth reduces ambiguity in each range bin. This
% results in taller peaks, allowing you to detect targets at higher CNR
% values.
% However, increasing the chrip bandwidth also requires an increase in the
% sampling rate. Computational time increases significantly with increasing
% bandwidth.
% in the BW = 1e6 case, we can see that this bandwidth is in fact "good
% enough" -- it is able to detect targets at CNR values of up to 46. Any
% lower bandwidth will simply result in an inability to resolve targets at
% any CNR.