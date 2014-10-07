% RCS_project

% Sample solution to the RCS project for ECE6272
% This code only does the pdf portions of the project.  Pdfs are obtained
% by collecting data for each of some number (typically 10 or more for good
% results) of sample random targets, then using all the data for all the
% targets to determine the pdf.
%
% Mark Richards, April 2009

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Input Section
%
% For variation with angle:
% Input number of angles and nominal range and frequency
% M=input('Enter # of angles: ');
% R=input('Enter nominal range (m): ');
% f=input('Enter frequency (Hz): ');
M=1000; R=10e3; f = 10e9;
Ntgt =input('Enter number of sample targets to average over: ');
N=50; % # of scatterers per target
%
% For variation with frequency:
% Input number of angles and nominal range and frequency
% Mf=input('Enter # of frequencies: ');
% f1=input('Enter start frequency (Hz): ');
Mf = 300;
f1 = 10e9;
df = 3e8/2/10;
% f2=input(['Enter freq step (Hz) (minimum=',num2str(df)/1e6,' MHz): ']);
f2 = 15e6;
fv = f1+(0:Mf-1)*f2;
%
%
% End User Input Section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(' '); disp('Many small scatterers case:')
x = zeros(N+1,Ntgt);
y = zeros(N+1,Ntgt);
echo=zeros(M,Ntgt);
RCS=zeros(M,Ntgt);
RCSdB=zeros(M,Ntgt);
phase=zeros(M,Ntgt);
count=zeros(1,Ntgt);
for n = 1:Ntgt
    % Set locations of scatterers.  They are uniformly
    % distributed within a box 5 m by 10 m centered at the origin.
    y(1:N,n)=5*rand(N,1)-2.5; x(1:N,n)=10*rand(N,1)-5;
    if (n==1)
        figure(1)
        plot(x(1:N,n),y(1:N,n),'*');
        xlabel('x'), ylabel('y')
        title('Sample Distribution of Scatterers, Many-Scatterer Case')
        axis('equal'); axis([-5,5,-2.5,2.5]);
    end
    % Amplitudes are fixed =1
    z=1;

    % Loop over aspect angle to do complex voltage estimates
    q=zeros(M,1);
    for k=1:M
        q(k)=(2*pi/M)*(k-1); % current aspect angle
        % Loop over individual point scatterers
        for p=1:N
            phasor=z*exp(j*4*pi*f*norm([x(p,n)-R*cos(q(k)),y(p,n)-R*sin(q(k))])/3e8);
            echo(k,n)=echo(k,n)+phasor;
        end
    end

    % The magnitude-squared of the total complex echo
    % is the RCS to within a constant
    RCS(:,n)=abs(echo(:,n)).^2;
    RCSdB(:,n)=10*log10(RCS(:,n));
    phase(:,n)=(180/pi)*angle(echo(:,n));
    if (n==1)
        figure(2)
        plot((180/pi)*q,RCSdB(:,n)); xlabel('aspect angle (degrees)');
        ylabel('relative RCS (dB)');
        title('Sample RCS vs. Aspect Angle, Many-Scatterer Case')
        axis([0,360,10*log10(N*z)-30,10*log10(N*z)+10]);

        figure(3)
        plot((180/pi)*q,phase(:,n));
        xlabel('aspect angle (degrees)'); ylabel('phase (deg)');
        title('Sample Phase vs. Aspect Angle, Many-Scatterer Case')
        axis([0,360,-180,180]);

        h=waitbar(0,'Fraction of targets completed ...');
    end
    waitbar(n/Ntgt,h);
end
close(h)

% Compute mean of angle RCS data from all trials, display on linear and dB scale
% Note that you do not average the dB data to get the average dB; average
% the linear data, then convert
RCS_avg = mean(RCS(:));
RCSdB_avg = 10*log10(RCS_avg);
disp(['RCS average over angle, many-scatterer-case, = ',num2str(RCS_avg)])
disp(['RCS average over angle in dB, many-scatterer-case, = ',num2str(RCSdB_avg)])

% Compute and plot histogram, overlay the theoretical exponential pdf
figure(4);
[count,centers]=hist(RCS(:),100);
bin_size = centers(2)-centers(1);
area = sum(count)*bin_size;
density = count/area;
bar(centers,density,1); shg
hold on

exp_pdf = exp(-centers/RCS_avg)/RCS_avg;  % theoretical pdf on same x-axis values
plot(centers,exp_pdf,'r','LineWidth',2)
xlabel('radar cross section (m^2)')
ylabel('relative probability')
title('PDF of RCS vs. Angle, Many-Scatterer Case')
hold off

pause

% Now do the variation with frequency
echof=zeros(Mf,Ntgt);
RCSf=zeros(Mf,Ntgt);
RCSfdB=zeros(Mf,Ntgt);
phasef=zeros(Mf,Ntgt);
countf=zeros(1,Ntgt);
densityf=zeros(100,Ntgt);
for n = 1:Ntgt

    % Use same scatterer distributions as before. Amplitudes are still
    % fixed = 1.
    z=1;

    % Loop over frequency to do complex voltage estimates
    for k=1:Mf
        % Loop over individual point scatterers
        for p=1:N
            phasor=z*exp(j*4*pi*fv(k)*norm([x(p,n)-R,y(p,n)])/3e8);
            echof(k,n)=echof(k,n)+phasor;
        end
    end

    RCSf(:,n)=abs(echof(:,n)).^2;
    RCSfdB(:,n)=10*log10(RCSf(:,n));
    phasef(:,n)=(180/pi)*angle(echof(:,n));
    if (n==1)
        figure(6)
        plot(fv/1e9,RCSfdB(:,n));
        xlabel('frequency (GHz)'); ylabel('relative RCS (dB)');
        title('Sample RCS vs. Frequency, Many-Scatterer Case')

        figure(7)
        plot(fv/1e9,phasef(:,n));
        xlabel('frequency (GHz)'); ylabel('phase (deg)');
        title('Sample Phase vs. Frequency, Many-Scatterer Case')
        axis([min(fv)/1e9,max(fv)/1e9,-180,180]);

        h=waitbar(0,'Fraction of targets completed ...');
    end
    waitbar(n/Ntgt,h);
end
close(h)

% Compute mean of frequency RCS data from all trials, display on linear and dB scale
% Note that you do not average the dB data to get the average dB; average
% the linear data, then convert
RCSf_avg = mean(RCSf(:));
RCSfdB_avg = 10*log10(RCSf_avg);
disp(['RCS average over frequency, many-scatterer-case, = ',num2str(RCSf_avg)])
disp(['RCS average over frequency in dB, many-scatterer-case, = ',num2str(RCSfdB_avg)])


figure(8);
[countf,centers]=hist(RCSf(:),100);
bin_size = centers(2)-centers(1);
area = sum(countf)*bin_size;
densityf = countf/area;
bar(centers,densityf,1); shg
hold on

exp_pdf = exp(-centers/RCSf_avg)/RCSf_avg;  % theoretical pdf on same x-axis values
plot(centers,exp_pdf,'r','LineWidth',2)
xlabel('radar cross section (m^2)')
ylabel('relative probability')
title('PDF of RCS vs. Frequency, Many-Scatterer Case')
hold off

pause

% Now repeat the whole exercise for the case of a
% dominant scatterer plus many small ones; and RCS
% of dominant = 10*(RCS of all the small ones together)

% Again I use the same sets of small scatterers and
% just add a different dominant scatterer to each set.  Also use same range,
% frequency, number of aspect angles.

disp(' '); disp('One dominant + many small scatterers cases:')

% We still have composite echo of the small scatterers
% as a function of aspect angle in the array 'echo';
% we just need to add the echo from the dominant.  However.
% I'll keep that in a separate array from the Rayleigh scatterers
% because I want to repeat this for a different dominant
% scatterer a little later

echo_dom = zeros(size(echo));
echof_dom = zeros(size(echof));

% Small scatterer amplitudes are fixed z = 1;
% big scatterer amplitude is therefore a2*N*z
% where a2 = ratio of big to sum of small
a2 = 10;
w = sqrt(a2*N*z);

% Loop over targets
for n = 1:Ntgt

    % Generate the position of the dominant scatterer
    y(N+1,n)=5*rand(1,1)-2.5; x(N+1,n)=10*rand(1,1)-5;
    if n==1
        figure(9)
        plot(x(1:N,n),y(1:N,n),'b*'); hold;
        plot(x(N+1,n),y(N+1,n),'ro'); hold off;
        xlabel('x'), ylabel('y')
        axis('equal'); axis([-5,5,-2.5,2.5]);
        title('Sample Distribution, One-Dominant + Many-Scatterer Case')
    end

    % Compute dominan's contributions over aspect angle
    for k=1:M
        % compute the dominant scatterer ontribution at each angle
        echo_dom(k,n)=w*exp(j*4*pi*f*norm([x(N+1,n)-R*cos(q(k)),y(N+1,n)-R*sin(q(k))])/3e8);
    end

    % Compute dominant's contributions over frequency
    for k=1:Mf
        echof_dom(k,n)=w*exp(j*4*pi*fv(k)*norm([x(N+1,n)-R,y(N+1,n)])/3e8);
    end

    % Compute RCS vs. angle and freq. with the dominant included
    RCS_dom(:,n)=abs(echo(:,n)+echo_dom(:,n)).^2;
    RCSdB_dom(:,n)=10*log10(RCS_dom(:,n));
    RCSf_dom(:,n)=abs(echof(:,n)+echof_dom(:,n)).^2;
    RCSfdB_dom(:,n)=10*log10(RCSf_dom(:,n));
    if (n==1)
        figure(10)
        plot((180/pi)*q,RCSdB_dom);
        xlabel('aspect angle (degrees)'); ylabel('relative RCS (dB)');
        title('RCS vs. Angle, Large Dominant+Many Case')
        axis([0,360,-10+20*log10(w),10+20*log10(w)]);

        figure(11)
        plot(fv/1e9,RCSfdB_dom);
        xlabel('frequency (GHz)'); ylabel('relative RCS (dB)');
        title('RCS vs. Frequency, Large Dominant+Many Case')
        axis([min(fv)/1e9,max(fv)/1e9,-10+20*log10(w),10+20*log10(w)]);
        
        h = waitbar(0,'Fraction of targets completed ...');
    end
    waitbar(n/Ntgt,h);
end
close(h)

% Compute a 100-bin histogram of the angle data;
% plot the theoretical chi-square distribution as well
figure(12);
[count,centers]=hist(RCS_dom(:),100);
bin_size = centers(2)-centers(1);
area = sum(count)*bin_size;
density = count/area;
bar(centers,density,1);
hold on

msig=mean(RCS_dom(:));  % compute mean of current linear-scale RCS data
rice_pdf = ((1+a2)/msig)*exp(-a2 - (1+a2)*centers/msig) ...
    .*besseli(0,2*sqrt(a2)*sqrt((1+a2)*(centers/msig)));  % theoretical pdf
plot(centers,rice_pdf,'r','LineWidth',2)
xlabel('radar cross section (m^2)')
ylabel('relative probability')
title('PDF of RCS vs. Aspect, Large Dominant+Many Case')
hold off

pause

% repeat for frequency data
figure(13);
[count,centers]=hist(RCSf_dom(:),100);
bin_size = centers(2)-centers(1);
area = sum(count)*bin_size;
density = count/area;
bar(centers,density,1);
hold on

msig=mean(RCSf_dom(:));  % compute mean of current linear-scale RCS data
rice_pdf = ((1+a2)/msig)*exp(-a2 - (1+a2)*centers/msig) ...
    .*besseli(0,2*sqrt(a2)*sqrt((1+a2)*(centers/msig)));  % theoretical pdf
plot(centers,rice_pdf,'r','LineWidth',2)
xlabel('radar cross section (m^2)')
ylabel('relative probability')
title('PDF of RCS vs. Frequency, Large Dominant+Many Case')
hold off

pause


% Once more, this time for the chi-square case Swerling uses
a2 = 1+sqrt(2);  % the difference is this particular value of 'a2'
w = sqrt(a2*N*z);

% Loop over targets
for n = 1:Ntgt

    % Generate the position of the dominant scatterer
    y(N+1,n)=5*rand(1,1)-2.5; x(N+1,n)=10*rand(1,1)-5;
    if n==1
        figure(14)
        plot(x(1:N,n),y(1:N,n),'b*'); hold;
        plot(x(N+1,n),y(N+1,n),'ro'); hold off;
        xlabel('x'), ylabel('y')
        axis('equal'); axis([-5,5,-2.5,2.5]);
        title('Sample Distribution, One-Dominant + Many-Scatterer Case')
    end

    % Compute dominan's contributions over aspect angle
    for k=1:M
        % compute the dominant scatterer ontribution at each angle
        echo_dom(k,n)=w*exp(j*4*pi*f*norm([x(N+1,n)-R*cos(q(k)),y(N+1,n)-R*sin(q(k))])/3e8);
    end

    % Compute dominant's contributions over frequency
    for k=1:Mf
        echof_dom(k,n)=w*exp(j*4*pi*fv(k)*norm([x(N+1,n)-R,y(N+1,n)])/3e8);
    end

    % Compute RCS vs. angle and freq. with the dominant included
    RCS_dom(:,n)=abs(echo(:,n)+echo_dom(:,n)).^2;
    RCSdB_dom(:,n)=10*log10(RCS_dom(:,n));
    RCSf_dom(:,n)=abs(echof(:,n)+echof_dom(:,n)).^2;
    RCSfdB_dom(:,n)=10*log10(RCSf_dom(:,n));
    if (n==1)
        figure(15)
        plot((180/pi)*q,RCSdB_dom(:,n));
        xlabel('aspect angle (degrees)'); ylabel('relative RCS (dB)');
        title('RCS vs. Angle, Large Dominant+Many Case')
        axis([0,360,-10+20*log10(w),10+20*log10(w)]);

        figure(16)
        plot(fv/1e9,RCSfdB_dom(:,n));
        xlabel('frequency (GHz)'); ylabel('relative RCS (dB)');
        title('RCS vs. Frequency, Large Dominant+Many Case')
        axis([min(fv)/1e9,max(fv)/1e9,-10+20*log10(w),10+20*log10(w)]);
        
        h = waitbar(0,'Fraction of targets completed ...');
    end
    waitbar(n/Ntgt,h);
end
close(h)

% Compute a 100-bin histogram of the angle data;
% plot the theoretical chi-square distribution as well
figure(17);
[count,centers]=hist(RCS_dom(:),100);
bin_size = centers(2)-centers(1);
area = sum(count)*bin_size;
density = count/area;
bar(centers,density,1);
hold on

msig=mean(RCS_dom(:));  % compute mean of current linear-scale RCS data
rice_pdf = ((1+a2)/msig)*exp(-a2 - (1+a2)*centers/msig) ...
    .*besseli(0,2*sqrt(a2)*sqrt((1+a2)*(centers/msig)));  % theoretical pdf
plot(centers,rice_pdf,'r','LineWidth',2)
xlabel('radar cross section (m^2)')
ylabel('relative probability')
title('PDF of RCS vs. Aspect, Swerling Small Dominant+Many Case')
hold off

pause

% repeat for frequency data
figure(18);
[count,centers]=hist(RCSf_dom(:),100);
bin_size = centers(2)-centers(1);
area = sum(count)*bin_size;
density = count/area;
bar(centers,density,1);
hold on

msig=mean(RCSf_dom(:));  % compute mean of current linear-scale RCS data
rice_pdf = ((1+a2)/msig)*exp(-a2 - (1+a2)*centers/msig) ...
    .*besseli(0,2*sqrt(a2)*sqrt((1+a2)*(centers/msig)));  % theoretical pdf
plot(centers,rice_pdf,'r','LineWidth',2)
xlabel('radar cross section (m^2)')
ylabel('relative probability')
title('PDF of RCS vs. Frequency, Swerling Small Dominant+Many Case')
hold off