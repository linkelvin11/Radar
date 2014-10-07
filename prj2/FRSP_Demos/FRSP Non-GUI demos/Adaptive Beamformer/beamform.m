% beamform.m

% for experimenting with simple beamforming calculations
%
% M. A. Richards, October 2004

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Input Section

lambda = 0.03; % wavelength
d = lambda/2;  % element spacing
dl = d/lambda;
N = 16;  % # of array elements
aoa_max = asin(1/2/dl);  % maximum "real space" AOA (radians)
window_on = false;  % true or false
threshold_beamspace = true;  % true or false
thresh = 0.05;

Nangle = 1000;  % # of angles for evaluating beam pattern
Nm1 = Nangle-1;

t_aoa = pi/180*(0); % target AOA (radians)
j_aoa1 = pi/180*(18); % jammer #1 AOA (radians)
j_aoa2 = pi/180*(-33); % jammer #2 AOA (radians)

SNR = 0;  % signal to noise ratio (dB)
JSR1 = +50;  % jammer #1 to noise ratio (dB)
JSR2 = +30;  % jammer #2 to noise ratio (dB)

%  End user input section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_n = 1;  % noise power
p_t = p_n*(10^(SNR/10));  % target power
p_j1 = p_t*(10^(JSR1/10)) ; % jammer power
p_j2 = p_t*(10^(JSR2/10)) ; % jammer power

%%%%%%%%%%%%%%%%%%%
% Pre-beamforming (raw data) SIR
SIR_before = p_t/(p_n + p_j1 + p_j2),'power';
disp(['SIR before beamforming =',num2str(db(SIR_before,'power')),' dB']);

% compute signal vectors

target = sqrt(p_t)*exp(j*2*pi*(0:N-1)'*d*sin(t_aoa)/lambda);
if (window_on)
    t = target.*taylorwin(length(target),4,-30);
else
    t = target;
end

j1 = sqrt(p_j1)*exp(j*2*pi*(0:N-1)'*d*sin(j_aoa1)/lambda);
j2 = sqrt(p_j2)*exp(j*2*pi*(0:N-1)'*d*sin(j_aoa2)/lambda);

% compute covariance matrix with and without jammers

R = p_n*eye(N);
Rj = p_n*eye(N) + p_j1*j1*(j1') + p_j2*j2*(j2');
disp(['Covariance matrix rank with jammers = ',num2str(rank(R))]);

% compute beamformer weight vector with and without jammers

w = R\conj(t);
wj = Rj\conj(t);

% compute and display beampattern

theta = -aoa_max + 2*aoa_max/Nm1*(0:Nm1);

W = zeros(Nangle,1);
Wj = zeros(Nangle,1);
for p = 1:Nangle
    W(p) = w.'*exp(-j*2*pi*(0:N-1)'*dl*sin(theta(p)));
    Wj(p) = wj.'*exp(-j*2*pi*(0:N-1)'*dl*sin(theta(p)));
end
Wp = db(abs(W),'voltage');
scale = 10*log10(N^2);
% scale = 0;
Wp = Wp - scale;
Wjp = db(abs(Wj),'voltage') - scale;
figure(1)
plot(180/pi*theta,Wp)
axis([-180*aoa_max/pi +180*aoa_max/pi -60  0])
% grid
xlabel('Angle of Arrival (degrees)');
ylabel('Normalized Array Response (dB)')
title('Unadapted Array Pattern')
vline(180/pi*[j_aoa1, j_aoa2])

figure(2)
plot(180/pi*theta,Wjp)
axis([-180*aoa_max/pi +180*aoa_max/pi -60  0])
% grid
xlabel('Angle of Arrival (degrees)');
ylabel('Normalized Array Response (dB)')
title('Fully Adaptive Array Pattern')
vline(180/pi*[j_aoa1, j_aoa2])

% Compute and display SIR after beamforming

SIR_after = real( norm(wj'*target)^2/(wj'*inv(Rj)*wj) );
disp(['SIR after beamforming = ',num2str(SIR_after),' = ',num2str(10*log10(SIR_after)),' dB']);
SIR_gain = SIR_after/SIR_before;
disp(['SIR gain =',num2str(SIR_gain),' = ',num2str(db(SIR_gain,'power')),' dB'])

% Now try unity gain ("distortionless") constraint solution from Guerci

kappa = t'*transpose(inv(Rj))*conj(t);
wj1 = wj/kappa;

Wj1 = zeros(Nangle,1);
for p = 1:Nangle
    Wj1(p) = wj1.'*exp(-j*2*pi*(0:N-1)'*dl*sin(theta(p)));
end
Wjp1 = db(abs(Wj1),'voltage');
figure(3)
plot(180/pi*theta,Wjp1)
axis([-180*aoa_max/pi +180*aoa_max/pi -60  0])
% grid
xlabel('Angle of Arrival (degrees)');
ylabel('Normalized Array Response (dB)')
title('Distortionless Beamformer Array Pattern')
vline(180/pi*[j_aoa1, j_aoa2])

SIRdis = real( norm(wj1'*target)^2/(wj1'*inv(Rj)*wj1) );
disp(['SIR distortionless = ',num2str(SIRdis),' = ',num2str(10*log10(SIRdis)),' dB']);


% Now try a post-DFT beamformer
% First, create a DFT matrix of the desired size

K = 16;
D = ones(K,N);
n = 0:N-1;
for k = 1:K-1
    D(k+1,:) = exp(-j*2*pi*k*n/K);
end

% Now transform the covariance matrix and form a new weight vector
Rjd = conj(D)*Rj*D.';
td = D*t;
wjd = Rjd\conj(td);
if threshold_beamspace
    null = find(abs(wjd) < thresh*max(abs(wjd)))  % discard DOFs with little weight
    wjd(null) = 0;
end
Wjd = zeros(Nangle,1);
for p = 1:Nangle
    Wjd(p) = wjd.'*(D*exp(-j*2*pi*(0:N-1)'*dl*sin(theta(p))));
end

scale = 10*log10(N^2);
Wjpd = db(abs(Wjd),'voltage') - scale;
figure(4)
plot(180/pi*theta,Wjpd)
axis([-180*aoa_max/pi +180*aoa_max/pi -60  0])
% grid
xlabel('Angle of Arrival (degrees)');
ylabel('Normalized Array Response (dB)')
title('Post-DFT (Beamspace) Array Pattern')
vline(180/pi*[j_aoa1, j_aoa2])

SIRdft = real( norm(wjd'*td)^2/(wjd'*inv(Rjd)*wjd) )/K^2;
disp(['SIRdft = ',num2str(SIRdft),' = ',num2str(10*log10(SIRdft)),' dB']);


