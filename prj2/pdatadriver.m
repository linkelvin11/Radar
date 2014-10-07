%% Description of custom functions

% makePDdata3(CN) generates data with CNR = CN

% makePDdata3(CN,LV) generates data with CNR = CN, and lowest velocity LV

% procPDdata2(CN,LV) processes data generated with a CNR of CN and lowest 
%   velocity with velocity = LV. this function calls makePDdata3(CN,LR) MTI
%   processing is left intact

% procPDdata3(CN,LV) processes data generated with CNR = CN (calls
%   makePDdata3(CN,LV). MTI processing is removed.

% procPDdata4 contains no changes other than using step instead of
% quadratic interpolation

% procPDdata5 uses linear interpolation


%% generate histograms for varying levels of clutter

figure
n = 3;
N = n*n;
d = 10;
for ii = d:d:d*N;
    y = makePDdata2(ii);
    y = sort(y);
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
suptitle('Histogram Changes vs CNR');

%The black lines denote the 25%, 50% and 75% points
%note that the 25%, 50% and 75% points move with the CNR

%% Test effectiveness of MTI
CNR = 55;
LV = 0.1;

procPDdata2(CNR,LV);
procPDdata3(CNR,LV);
close all;
% check in command window for detected targets
% CNR = 01: All targets resolved
% CNR = 02: 1 fa w/o TPC, target @ 3.1 md w/ TPC
% CNR = 05: 2 fa w/o TPC, target @ 3.1 md w/ TPC
% CNR = 10: 2 fa w/o TPC, target @ 3.1 md w/ TPC
% CNR = 15: 2 fa w/o TPC, target @ 3.1 md w/ TPC
% CNR = 20: 32 fa w/o TPC, target @ 3.1 md w/ TPC

% after this point, we will discard results from non-TPC test
% the following results refer to the TPC active case

% CNR = 25: target @ 3.1 md w/ TPC
% CNR = 30: target @ 3.1 md w/ TPC
% CNR = 35: 1 fa @ 3.52km 93.67m/s, target @ 3.1 md w/ TPC
% CNR = 40: target @ 3.1 md w/ TPC
% CNR = 45: 1 fa @ 3.52km 94.76m/s, target @ 3.1 md w/ TPC
% CNR = 50: 1 fa @ 4.05km 20m/s, target @ 3.1 md w/ TPC
% CNR = 55: 35 fa. target @ 3.1 md w/ TPC

% as CNR increases, the clutter spreads in doppler, causing the algorithm
% to report many false alarms. This is the data from the CNR = 55 case:
% ESTIMATED PARAMETERS OF DETECTED TARGETS:
% Number     Rel Amp (dB)     Range (km)     Vel (m/s)
%    1         -27.3             2.3          -112.6
%    2         -5.99            2.46          -18.31
%    3         -5.09            2.61          -17.31
%    4         -16.7             2.7          -55.46
%    5         -7.58            2.75           19.91
%    6        -0.316            2.95           14.87
%    7         -1.41            2.99           -13.9
%    8         -5.24            3.01           18.85
%    9         -5.17            3.09           19.83
%   1e+01         -6.47            3.11          -21.71
%   1e+01             0            3.14          -13.94
%   1e+01         -5.26            3.17          -17.74
%   1e+01         -6.11            3.21           19.36
%   1e+01         -6.35            3.25           18.96
%   2e+01         -5.39            3.31           20.06
%   2e+01         -7.23            3.36           19.77
%   2e+01         -3.99            3.43          -18.56
%   2e+01         -5.01            3.46          -19.44
%   2e+01         -17.1             3.5           93.58
%   2e+01         -5.64            3.54          -18.63
%   2e+01         -1.14            3.57           13.93
%   2e+01         -6.76             3.6          -20.22
%   2e+01        -0.978            3.65            18.2
%   2e+01         -1.75            3.68          -18.46
%   3e+01         -1.57             3.7          -19.59
%   3e+01         -1.64            3.78          -18.82
%   3e+01         -6.75            3.82           21.69
%   3e+01        -0.674            3.89           -17.5
%   3e+01        -0.996            3.94              17
%   3e+01         -6.13            3.96          -19.04
%   3e+01         -4.75            3.99           20.95
%   3e+01         -7.51            4.03           21.99
%   3e+01          -2.3            4.12          -20.09
%   3e+01         -3.02            4.16          -18.93
%   4e+01         -5.64             4.2          -20.73 
% as seen above, the false alarm targets have significantly lower
% velocities as compared to the genuine targets. ~14-22m/s
% additionally, in the CNR = 35,45,50 cases, the false alarm was simply a
% single target detected twice. This is probably due to the spread in
% doppler of the clutter (forced the target to peak twice)

%% find a suitable test for moving targets
close all;

CNR = 11;
LV = 0.1;

if LV > 0.25
    LV = 0.25;
end

procPDdata2(CNR, LV);
%%
procPDdata3(CNR, LV);



% Interesting Result: moving the target closer to the clutter is actually
% worse for the case with MTI enabled. The pulse cancellation seems to
% eliminate the slower target along with the clutter. (slow targets are
% indistinguishable from clutter
% for the target to be detected, velocity must be ~>45 m/s


%% changing chirp bandwidth
%default bandwidth 10e6
clear
makePDdata4(8e6);
procPDdata

% as the bandwidth increases, so does the computational time required to
% resolve the targets.
% however, this increase in computational time seems due to the actual
% plotting of data, not in the processing of said data.

% increasing chirp bandwidth does not help us detect targets any better,
% they simply have a smaller unambiguous range

% lowering the bandwidth down to half of the current W generates a false
% alarm around the target at 3.5 km.
% Number     Rel Amp (dB)     Range (km)     Vel (m/s)
%    1         -10.3             2.3          -112.8
%    2        -0.243             2.7          -56.44
%    3             0             3.5           93.74
%    4         -15.4            3.55           93.61 