% RCS10_project

% Does a simple-minded relative RCS polar plot for a "complex
% target" consisting of a random collection of point scatterers.
% No multiple bounce or any such nonsense.
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

% plot the scatterer distribution for the first target as an example
figure(1)
plot(x(:,1),y(:,1),'*');
xlabel('x'), ylabel('y')
title('Distribution of Scatterers')
axis('equal'); axis([-5,5,-2.5,2.5]);

% Amplitudes are fixed = 1
z=1;

% Input number of angles and nominal range and frequency
% M=input('Enter # of angles: ');
% R=input('Enter nominal range (m): ');
% f=input('Enter frequency (Hz): ');
M=round(360/0.5); R=10e3; f = 10e9;
% Loop over aspect angle to do complex voltage estimates
q=zeros(M); echo=zeros(M,Nt);

% Loop over radar-target aspect angles
for k=1:M
    q(k)=(2*pi/M)*(k-1); % current aspect angle
    % Loop over targets and individual point scatterers
    % each column of 'echo' is a different target
    for r = 1:Nt
        for p=1:N
            phasor=z*exp(j*4*pi*f*norm([x(p,r)-R*cos(q(k)),y(p,r)-R*sin(q(k))])/3e8);
            echo(k,r)=echo(k,r)+phasor;
        end
    end
end

% The magnitude-squared of the total complex echo
% is the RCS to within a constant
RCS=abs(echo).^2;
RCSdB=10*log10(RCS);

% Plot the dB data for the first target as an example
figure(2)
plot((180/pi)*q,RCSdB(:,1)); xlabel('aspect angle (degrees)');
ylabel('relative RCS (dB)');
title('RCS vs. Aspect Angle, Many-Scatterer Case')
axis([0,360,10*log10(N*z)-30,10*log10(N*z)+10]);

% Compute a 100-bin histogram of the data for each target;
% plot the theoretical exponential distribution as well
%
density = zeros(100,Nt);
centers = [];
for r = 1:Nt
    if (isempty(centers))
        [count,centers]=hist(RCS(:,r),100);
    else
        [count,centers]=hist(RCS(:,r),centers);
    end
    bin_size = centers(2)-centers(1);
    area = sum(count)*bin_size;
    density(:,r) = count/area;
    % plot the histogram of the data for each target, and an overlaid
    % exponential for each target, based only on that target's data
    figure(3);
    bar(centers,density(:,r),1);
    hold on
    msig=mean(RCS(:,r));  % compute mean of current linear-scale RCS data
    exp_pdf = exp(-centers/msig)/msig;  % theoretical pdf on same x-axis values
    plot(centers,exp_pdf,'r','LineWidth',2)
    xlabel('radar cross section (m^2)')
    ylabel('relative probability')
    title('PDF of RCS vs. Aspect, Many-Scatterer Case')
    hold off
    pause(0.5)
end

% Now average the separate histograms and plot with an overlaid exponential
% based on the mean of all the data for all of the targets
figure(4);
bar(centers,mean(density,2),1);
hold on
msig=mean(RCS(:));  % compute mean of all linear-scale RCS data
exp_pdf = exp(-centers/msig)/msig;  % theoretical pdf on same x-axis values
plot(centers,exp_pdf,'r','LineWidth',2)
xlabel('radar cross section (m^2)')
ylabel('relative probability')
title('PDF of RCS vs. Aspect, Many-Scatterer Case, All Targets')
hold off

% As a sanity check, do a single histogram of all of the RCS data for all
% targets, using the sam bin boundaries as above.  This should match figure
% 4 exactly.
figure(5)
[count,centers]=hist(RCS(:),centers);
bin_size = centers(2)-centers(1);
area = sum(count)*bin_size;
density_all = count/area;
bar(centers,density_all,1);
hold on
plot(centers,exp_pdf,'r','LineWidth',2)
xlabel('radar cross section (m^2)')
ylabel('relative probability')
title('PDF of RCS vs. Aspect, Many-Scatterer Case, All Targets')
hold off


% Now do the variation with frequency

% Input number of angles and nominal range and frequency
% Mf=input('Enter # of frequencies: ');
% f1=input('Enter start frequency (Hz): ');
Mf = 300;
f1 = 10e9;
df = 3e8/2/10;
% f2=input(['Enter freq step (Hz) (minimum=',num2str(df)/1e6,' MHz): ']);
f2 = 15e6;
fv = f1+(0:Mf-1)*f2;

% Loop over radar frequencies
echof=zeros(Mf,Nt);
for k=1:Mf
    % Loop over targets and individual point scatterers
    % each column of 'echof' is a different target
    for r = 1:Nt
        for p=1:N
            phasor=z*exp(j*4*pi*fv(k)*norm([x(p,r)-R,y(p,r)])/3e8);
            echof(k,r)=echof(k,r)+phasor;
        end
    end
end

% The magnitude-squared of the total complex echo
% is the RCS to within a constant
RCSf=abs(echof).^2;
RCSfdB=10*log10(RCSf);

% Plot the dB data for the first target as an example
figure(6)
plot(fv,RCSfdB(:,1)); xlabel('RF frequency (Hz)');
ylabel('relative RCS (dB)');
title('RCS vs. RF Frequency, Many-Scatterer Case')
axis([0,360,10*log10(N*z)-30,10*log10(N*z)+10]);

% Compute a 100-bin histogram of the data for each target;
% plot the theoretical exponential distribution as well
%
densityf = zeros(100,Nt);
centersf = [];
for r = 1:Nt
    if (isempty(centersf))
        [countf,centersf]=hist(RCSf(:,r),100);
    else
        [countf,centersf]=hist(RCSf(:,r),centersf);
    end
    bin_sizef = centersf(2)-centersf(1);
    areaf = sum(countf)*bin_sizef;
    densityf(:,r) = countf/areaf;
    % plot the histogram of the data for the first target, and an overlaid
    % exponential for that target, based only on that target's data
        figure(7);
        bar(centersf,densityf(:,r),1);
        hold on
        msigf=mean(RCSf(:,r));  % compute mean of current linear-scale RCS data
        exp_pdf = exp(-centersf/msigf)/msigf;  % theoretical pdf on same x-axis values
        plot(centersf,exp_pdf,'r','LineWidth',2)
        xlabel('radar cross section (m^2)')
        ylabel('relative probability')
        title('PDF of RCS vs. RF, Many-Scatterer Case, 1 Target')
        hold off
        pause(0.5)
end

% Now average the separate histograms and plot with an overlaid exponential
% based on the mean of all the data for all of the targets
figure(8);
bar(centersf,mean(densityf,2),1);
hold on
msigf=mean(RCSf(:));  % compute mean of all linear-scale RCS data
exp_pdf = exp(-centersf/msigf)/msigf;  % theoretical pdf on same x-axis values
plot(centersf,exp_pdf,'r','LineWidth',2)
xlabel('radar cross section (m^2)')
ylabel('relative probability')
title('PDF of RCS vs. RF, Many-Scatterer Case, All Targets')
hold off

% As a sanity check, do a single histogram of all of the RCS data for all
% targets, using the sam bin boundaries as above.  This should match figure
% 8 exactly.
figure(9)
[countf,centersf]=hist(RCSf(:),centersf);
bin_sizef = centersf(2)-centersf(1);
areaf = sum(countf)*bin_sizef;
density_allf = countf/areaf;
bar(centersf,density_allf,1);
hold on
plot(centersf,exp_pdf,'r','LineWidth',2)
xlabel('radar cross section (m^2)')
ylabel('relative probability')
title('PDF of RCS vs. RF, Many-Scatterer Case, All Targets')
hold off

