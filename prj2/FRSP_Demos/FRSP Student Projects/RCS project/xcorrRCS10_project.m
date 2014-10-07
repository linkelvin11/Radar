% xcorrRCS10_project

% same as RCS_project but for computing correlation around zero degrees
%
% Mark Richards, September 2006
%
% Updated for averaging over multiple target, September 2010

clear all

close all

Nt=10; % # of targets to average over

% Set number and locations of scatterers.  They are uniformly
% distributed within a box 5 m by 10 m centered at the origin.

% N=input('Enter number of scatterers to use: ');
N=50;

% each column of 'x' or 'y' are the scatterer corrdinates for a single
% target
y=5*rand(N,Nt)-2.5; x=10*rand(N,Nt)-5;

% Amplitudes are fixed = 1
z=1;

% Input number of angles and nominal range and frequency
% thetamax=input('Enter +/- range of angles (degrees): ');
% Dtheta = input('Enter angle increment (degrees): ');
% M = round(thetamax/Dtheta);
% R=input('Enter nominal range (m): ');
% f=input('Enter frequency (Hz): ');

thetamax = 3; Dtheta = 0.02; R = 10e3; f = 10e9; M = round(thetamax/Dtheta);

% Loop over aspect angle to do complex voltage estimates
q=zeros(2*M+1,1); echo=zeros(2*M+1,Nt);

% Loop over radar-target aspect angles
for k=1:2*M+1
    q(k)=(pi/180)*((k-M-1)*Dtheta); % current aspect angle
    % Loop over targets and individual point scatterers
    % each column of 'echo' is a different target
    for r = 1:Nt
        for p=1:N
            phasor=z*exp(j*4*pi*f*norm([x(p,r)-R*cos(q(k)),y(p,r)-R*sin(q(k))])/3e8);
            echo(k,r)=echo(k,r)+phasor;
        end
    end
end

% autocorr of the complex data for each target
s = zeros(4*M+1,Nt);
lag = Dtheta*(-2*M:+2*M);
for r = 1:Nt
    s(:,r) = xcorr(echo(:,r));
    % plot autocorr of each target with markers
    figure(1)
    plot(lag,abs(s(:,r))/max(abs(s(:,r)))); xlabel('aspect angle (degrees)');
    ylabel('normalized magnitude of autocorrelation');
    title('Autocorrelation, 1 Target')
    grid;
    % plot theoretical markers
    thetacorr = 3e8/2/5/f*180/pi;
    hold on
    plot([-thetacorr,-thetacorr],[0,1],'-r');
    plot([thetacorr,thetacorr],[0,1],'-r');
    hold off
    pause(0.5)
end

% Now do the average complex autocorrelation and plot and mark its
% magnitude
s_avg = mean(s,2);
figure(2)
plot(lag,abs(s_avg)/max(abs(s_avg))); xlabel('aspect angle (degrees)');
ylabel('normalized magnitude of autocorrelation');
title('Autocorrelation, 10 Targets')
grid;
% plot theoretical markers
thetacorr = 3e8/2/5/f*180/pi;
hold on
plot([-thetacorr,-thetacorr],[0,1],'-r');
plot([thetacorr,thetacorr],[0,1],'-r');
hold off


pause



% Now do the autocorrelation with frequency

% fmax = input('Enter +/- range of freequencies: ')
% df = input('Enter frequency step: ')
fmax = 30e6;
df = 1e6;
Mf = round(fmax/df);

% Loop over radar frequencies
fx = zeros(2*Mf+1,1);
echof=zeros(2*Mf+1,Nt);

for k=1:2*Mf+1;
    fx(k) = (k-Mf-1)*df;
    % Loop over targets and individual point scatterers
    % each column of 'echof' is a different target
    for r = 1:Nt
        for p=1:N
            phasor=z*exp(j*4*pi*fx(k)*norm([x(p,r)-R,y(p,r)])/3e8);
            echof(k,r)=echof(k,r)+phasor;
        end
    end
end

% autocorr of the complex data for each target
sf = zeros(4*Mf+1,Nt);
lagf = df*(-2*Mf:+2*Mf);
for r = 1:Nt
    sf(:,r) = xcorr(echof(:,r));
    % plot autocorr of each target with markers
    figure(3)
    plot(lagf/1e6,abs(sf(:,r))/max(abs(sf(:,r)))); xlabel('RF Frequency (MHz)');
    ylabel('normalized magnitude of autocorrelation');
    title('Autocorrelation, 1 Target')
    grid;
    % plot theoretical markers
    fcorr = 3e8/2/10;
    hold on
    plot([-fcorr,-fcorr]/1e6,[0,1],'-r');
    plot([fcorr,fcorr]/1e6,[0,1],'-r');
    hold off
    pause(0.5)
end

% Now do the average complex autocorrelation and plot and mark its
% magnitude
sf_avg = mean(sf,2);
figure(4)
plot(lagf/1e6,abs(sf_avg)/max(abs(sf_avg))); xlabel('RF Frequency (MHz)');
ylabel('normalized magnitude of autocorrelation');
title('Autocorrelation, 10 Targets')
grid;
% plot theoretical markers
fcorr = 3e8/2/10;
hold on
plot([-fcorr,-fcorr]/1e6,[0,1],'-r');
plot([fcorr,fcorr]/1e6,[0,1],'-r');
hold off