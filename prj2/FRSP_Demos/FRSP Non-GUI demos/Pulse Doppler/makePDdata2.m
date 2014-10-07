%
% makePDdata
%
%  make pulse Doppler radar project data
%
%  Written by J. H. McClellan
%  Modified by M. A. Richards
%
close all;
clear, hold off
format compact
J = sqrt(-1);


% User Input Section  ########################################
% ############################################################

% Get root file name for saving results
%file=input('Enter root file name for data and listing files: ','s');
file = 'test';
T = 10e-6;     % pulse length, seconds
W = 10e6;      % chirp bandwidth, Hz
fs = 12e6;     % chirp sampling rate, Hz; oversample by a little


Np = 20;              % # of pulses
jkl = 0:(Np-1);       % pulse index array
PRF = 25.0e3;         % PRF in Hz
PRI = (1/PRF);        % PRI in sec
T_0 = PRI*jkl;        % relative start times of pulses, in sec
g = ones(1,Np);       % gains of pulses
T_out = [12 38]*1e-6; % start and end times of range window in sec
T_ref = 0;            % system reference time in usec
fc = 10e9;            % RF frequency in Hz; 10 GHz is X-band

% Compute unambiguous Doppler interval in m/sec
% Compute unambiguous range interval in meters

vua = 3e8*PRF/(2*fc);
rmin = 3e8*T_out(1)/2;
rmax = 3e8*T_out(2)/2;
rua = rmax-rmin;

% Define number of targets, then range, amplitude, and
% radial velocity of each

Ntargets = 4;
del_R = (3e8/2)*( 1/fs )/1e3;                   % in km
ranges = [2.3 2.7 3.1 3.5]*1e3;               % in km
SNR =    [10 20 15 20];               % dB
vels = [-0.3  -0.15  +0.1  +0.25]*vua;          % in m/sec

% End User Input Section  ####################################
% ############################################################



% form radar chirp pulse

fprintf('\nPulse length = %g microseconds\n',T/1e-6)
fprintf('Chirp bandwidth = %g Mhz\n',W/1e6)
fprintf('Sampling rate = %g Msamples/sec\n',fs/1e6)
s = git_chirp(T,W,fs/W);
figure(1)
plot((1e6/fs)*(0:length(s)-1),[real(s) imag(s)])
title('Real and Imaginary Parts of Chirp Pulse')
xlabel('time (usec)')
ylabel('amplitude')
grid

fprintf('\nWe are simulating %g pulses at an RF of %g GHz',Np,fc/1e9)
fprintf('\nand a PRF of %g kHz, giving a PRI of %g usec.',PRF/1e3,PRI/1e-6)
fprintf('\nThe range window limits are %g to %g usec.\n', ...
    T_out(1)/1e-6,T_out(2)/1e-6)

fprintf('\nThe unambiguous velocity interval is %g m/s.',vua)
fprintf('\nThe range window starts at %g km.',rmin/1e3)
fprintf('\nThe range window ends at %g km.',rmax/1e3)
fprintf('\nThe unambiguous range interval is %g km.\n\n',rua/1e3)


fprintf('\nThere are %g targets with the following parameters:',Ntargets)
for i = 1:Ntargets
  fprintf('\n  range=%5.2g km, SNR=%7.3g dB, vel=%9.4g m/s', ...
           ranges(i)/1e3,SNR(i),vels(i) )
end

% Now form the range bin - pulse number data map

disp(' ')
disp(' ')
disp('... forming signal component')
y = radar(s,fs,T_0,g,T_out,T_ref,fc,ranges,SNR,vels);

% add thermal noise with unit power

disp('... adding noise')
%randn('seed',77348911);
[My,Ny] = size(y);
nzz = (1/sqrt(2))*(randn(My,Ny) + J*randn(My,Ny));
y = y + nzz;

% create log-normal (ground) "clutter" with specified C/N and
% log-normal standard deviation for amplitude, uniform phase
% Clutter is uncorrelated in range, fully correlated in pulse #

disp('... creating clutter')
CN = 20;         % clutter-to-noise ratio in first bin
SDxdB = 3;       % in dB (this is NOT the sigma of the complete clutter)
ncc=10 .^((SDxdB*randn(My,Ny))/10);
ncc = ncc.*exp( J*2*pi*rand(My,Ny) );

% Force the power spectrum shape to be Gaussian

disp('... correlating and adding clutter')
G  = exp(-(0:4)'.^2/1.0);
G = [G;zeros(Ny-2*length(G)+1,1);G(length(G):-1:2)];

for i=1:My
  ncc(i,:)=ifft(G'.*fft(ncc(i,:)));
end
 
% rescale clutter to have desired C/N ratio
pcc = mean(mean(abs(ncc).^2));
ncc = sqrt((10^(CN/10))/pcc)*ncc;

% Now weight the clutter power in range for assume R^2 (beam-limited) loss
cweight = T_out(1)*((T_out(1) + (0:My-1)'*(1/fs)).^(-1));
cweight = cweight*ones(1,Np);
ncc = ncc.*cweight;

y = y + ncc;

[My,Ny]=size(y);
d=(3e8/2)*((0:My-1)*(1/fs) + T_out(1))/1e3;
figure(2)
plot(d,db(y,'voltage'))
xlabel('distance (km)')
ylabel('amplitude (dB)')
grid

figure
hist(db(y,'voltage'));

% Save the data matrix in specified file.
% Save the student version in the mystery file.  The student version doesn
% not include the target parameters, jsut the radar parameters and the
% data.
% Also save all parameter value displays in corresponding file

data_file=[file,'.mat'];
mystery_file=[file,'_mys.mat'];
listing_file=[file,'.lis'];

eval(['save ',data_file,' J T W fs s Np PRF PRI T_out fc vua', ...
    ' rmin rmax rua Ntargets ranges vels SNR y']);

eval(['save ',mystery_file,' J T W fs s Np PRF T_out fc y']);

fid=fopen(listing_file,'w');

fprintf(fid,['\rDESCRIPTION OF DATA IN FILE ',file,'.mat AND ',file,'_mys.mat\r\r']);
fprintf(fid,'\rPulse length = %g microseconds\r',T/1e-6);
fprintf(fid,'Chirp bandwidth = %g Mhz\r',W/1e6);
fprintf(fid,'Sampling rate = %g Msamples/sec\r',fs/1e6);
fprintf(fid,'\rWe are simulating %g pulses at an RF of %g GHz',Np,fc/1e9);
fprintf(fid,'\rand a PRF of %g kHz, giving a PRI of %g usec.',PRF/1e3,PRI/1e-6);
fprintf(fid,'\rThe range window limits are %g to %g usec.\r', ...
    T_out(1)/1e-6,T_out(2)/1e-6);
fprintf(fid,'\rThe unambiguous velocity interval is %g m/s.',vua);
fprintf(fid,'\rThe range window starts at %g km.',rmin/1e3);
fprintf(fid,'\rThe range window ends at %g km.',rmax/1e3);
fprintf(fid,'\rThe unambiguous range interval is %g km.\r\r',rua/1e3);
fprintf(fid,'\rThere are %g targets with the following parameters:', ...
  Ntargets);
for i = 1:Ntargets
  fprintf(fid,'\r  range=%5.2g km, SNR=%7.3g dB, vel=%9.4g m/s', ...
           ranges(i)/1e3,SNR(i),vels(i) );
end

fclose(fid);

fprintf(['\n\nData is in file ',data_file])
fprintf(['\nStudent data is in file ',mystery_file])
fprintf(['\nListing is in file ',listing_file,'\n\n'])

shg
