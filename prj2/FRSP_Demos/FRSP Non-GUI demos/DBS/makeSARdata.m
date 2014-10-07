%
% makeSARdata
%
%  make SAR project data
%
%  Written by M. A. Richards
%
clear all, close all, hold off
format compact
clc

c = 3e8;                               % velocity of light in m/sec


%######################################################
% User input section  #################################

% need a file name to store data, e.g. 'myfile' or 'temp' or 'sardata'.
% Data will then be in 'myfile.dat', etc.
file=input('Enter root file name for data and listing files: ','s');

DCR = 10;           % cross-range resolution, m
DR = 10;            % range resolution, m
Rcrp = 30000;        % range to swath center
v = 150;             % platform velocity, m/s
fc = 10e9;           % RF frequency in Hz
lambda = c/fc;        % wavelength, m
tau = 5e-6;           % pulse length, seconds (bandwidth set by resolution)
Daz = 0.2;            % antenna azimuth size, m
thetaaz = lambda/Daz; % azimuth beamwidth, radians
BWdopp = 2*v/lambda*thetaaz;     % slow time Doppler bandwidth, Hz
Ls = 3000;            % swath depth, m
oversample_st = 3;    % slow time oversample factor; higher makes
                      % prettier pictures but larger data sets
oversample_ft = 3;    % fast time oversample factor, similar to slow time    

% Define target locations, one row of (x,R) coordinates per target
% Ranges are relative to the CRP range (Rcrp) above.

% coords = [0,0];         % a single point target at the CRP

coords = ...            % a grid of 9 point targets
    [-1000,-1000;
    -1000,0;
    -1000,+1000;
    0,-1000;
    0,0;
    0,+1000;
    +1000,-1000;
    +1000,0;
    +1000,+1000];


% End user input section #####################################
% ############################################################

% check consistency of cross-range resolution and aperture size
if (DCR < Daz/2)
    disp('Requested cross-range resolution below minimum');
    disp(['Requested cross-range resolution = ',num2str(DCR)]);
    disp(['Minimum cross-range resolution = (Daz/2) = ', num2str(Daz/2)]);
    disp('')
    return
end


% Compute required aperture time
Ta = lambda*Rcrp/2/v/DCR;
Rc = (v*Ta)^2/8/Rcrp;        % nominal range curvature, m

fprintf('\nRange resolution = %g meters',DR);
fprintf('\nCross-range resolution = %g meters',DCR);
fprintf('\nRange to central reference point = %g km',Rcrp/1e3);
fprintf('\nSwath length = %g m',Ls);
fprintf('\nAperture time = %g seconds',Ta);
fprintf('\nRange curvature at CRP = %g meters\n\n',Rc);

% Compute PRF, see if both range ambiguity and Doppler bandwidth
% constraints can be met

PRFgoal = oversample_st*BWdopp;     % PRF in Hz
PRF = min(PRFgoal,c/2/Ls);
if (PRF < PRFgoal)
    if (PRF < BWdopp)
        disp('No viable PRF');
        disp(['Minimum PRF (Doppler BW constraint) = ',num2str(BWdopp)]);
        disp(['Maximum PRF = (Swath ambiguity constraint) = ', ...
                num2str(2*Ls/c)]);
        disp('')
        return
    end
end
oversample_st = PRF/BWdopp;
PRI = (1/PRF);        % PRI in sec

% compute number of pulses, aperture positions

Npulses = round(Ta/PRI);    % # of pulses
jkl = 0:(Npulses-1);        % pulse index array
T_0 = PRI*jkl;              % relative start times of pulses, in sec
T_0 = T_0 - max(T_0)/2;     % recenters pulse times symmetrically +/- around zero
u = v*T_0;                  % aperture positions

% form radar chirp pulse

W = c/2/DR;       % chirp bandwidth, Hz
fs = oversample_ft*W;       % chirp sampling rate, Hz
s = git_chirp(tau,W,fs/W);
Ns = length(s);

fprintf('\nPulse length = %g microseconds\n',tau/1e-6)
fprintf('Chirp bandwidth = %g Mhz\n',W/1e6)
fprintf('Time-bandwidth product = %g\n',W*tau)
fprintf('Fast time sampling rate = %g Msamples/sec\n',fs/1e6)

figure(1)
plot((1e6/fs)*(0:length(s)-1),[real(s) imag(s)])
title('Real and Imaginary Parts of Chirp Pulse')
xlabel('time (usec)')
ylabel('amplitude')
grid

fprintf('\nWe are simulating %g pulses at an RF of %g GHz',Npulses,fc/1e9)
fprintf('\nand a PRF of %g kHz, giving a PRI of %g usec.\n\n',PRF/1e3,PRI/1e-6)

Ntargets = size(coords,1);

fprintf('\nThere are %g targets with the following locations relative to the CRP:',Ntargets)
for nt = 1:Ntargets
    fprintf('\n  x=%+8.1f m, R=%+8.1f m',coords(nt,1),coords(nt,2));
end
fprintf('\n\n')

% Precompute range vs. pulse number.  Radar coordinates given by x=u,r=0;
% target coordinates relative to the CRP are in array 'coords'.

range = zeros(Ntargets,Npulses);
for nt = 1:Ntargets
    for m = 1:Npulses
        range(nt,m) = sqrt((coords(nt,1)-u(m))^2 + (Rcrp+coords(nt,2))^2);
    end
end
figure(2)
plot(u,range)
xlabel('aperture position (m)')
ylabel('range (m)')
title('Range vs. Aperture Position for each Scatterer')
grid

Rmin = min(min(range));
Rmax = max(max(range));

% check to see if min and max range fall in the defined swath
if ( (Rmin < Rcrp-Ls/2) | (Rmax > Rcrp+Ls/2) )
    disp('Target ranges outside of swath!')
    disp(['Swath limits are ',num2str(Rcrp-Ls/2),' to ',num2str(Rcrp+Ls/2)])
    disp(['Min and max range are ',num2str(Rmin),' and ',num2str(Rmax)])
    return
end

% Define the range window
T_out = [2*(Rcrp-Ls/2)/c, 2*(Rcrp+Ls/2)/c+tau];

fprintf('\nThe range window limits are %g to %g usec.\n', ...
    T_out(1)/1e-6,T_out(2)/1e-6)
fprintf('\nThe range window starts at %g km.',(Rcrp-Ls/2)/1e3)
fprintf('\nThe range window ends at %g km.',(Rcrp+Ls/2)/1e3)
fprintf('\nThe unambiguous range interval is %g km.\n\n',c/2/PRF/1e3)


% Now build the fast time/slow time data matrix
disp('Now forming signal matrix')

delta_t = 1/fs;                        % sampling interval (sec)
t_y = [ T_out(1):delta_t:T_out(2) ]';  % output sampling times (sec)
T_p = Ns*delta_t;                      % length of input pulse (sec)

% ensure that all vectors are column vectors

s=s(:);  T_0=T_0(:);

% determine the quadratic phase modulation parameters for
% later interpolation of pulse samples

t_s = delta_t*[0:(Ns-1)]';
s_ph = unwrap(angle(s));
warning off MATLAB:polyfit:RepeatedPointsOrRescale  % kills an annoying MATLAB warning about
% the polyfit function
q = polyfit(t_s,s_ph,2);

% check result using correlation coefficient
sfit = polyval(q,t_s);
if (s_ph'*sfit)/norm(s_ph)/norm(sfit) < 0.99
    disp('pulse phase is not linear or quadratic')
    disp('')
    return
end

%
%---  Form (initially empty) output matrix ---
%
Nr = length(t_y);       % output samples in a matrix
y = zeros(Nr,Npulses);

for m = 1:Npulses      % loop over aperture positions
    if (mod(m,100)==1)
        disp(['  ... processing pulse #',int2str(m),' of ',int2str(Npulses)])
    end
    for n = 1:Ntargets   % loop over targets
        r = range(n,m);
        
        % Compute start and end time of reflected pulse at receiver,
        % ensure that it falls entirely within the range (time) window
        tmin = 2*r/c;tmax = tmin + T_p;
        if ( (tmin < T_out(1)) | (tmax > T_out(2)) )
            fprintf('\nEcho from target #%g at range %g km',n,r/1e3)
            fprintf('\nis COMPLETELY OUT OF the range window')
            fprintf('\non pulse #%g.\n',m)
            return
        end
        
        % Figure out which sample locations in the output grid contain
        % reflected pulse
        t_vals = t_y - tmin;
        n_out = find( t_vals >= 0  &  t_vals < T_p );
        
        % Place pulse into output matrix.
        % Stop-and-hop assumed, so only phase modulation is due to geometry.
		% All pulses have unit amplitude.
		% There is no adjustment for R^4, either.
      
        y(n_out,m) = y(n_out,m) + ...
            exp( -j*2*pi*fc*tmin ) ...
            .* exp( j*polyval(q,t_vals(n_out)) );
        
    end  % end of loop over targets
end    % end of loop over pulses

figure(3)
imagesc(1:Npulses,t_y,real(y))
xlabel('pulse number')
ylabel('fast time (sec)')
title('Real Part of Data Matrix')


% Save the data matrix in specified file.
% Save the student version in the mystery file.
% Also save all parameter value displays in corresponding file

data_file=[file,'.mat'];
mystery_file=[file,'_mys.mat'];
listing_file=[file,'.lis'];

eval(['save ',data_file,' c tau W fs s Npulses PRF PRI T_out fc', ...
        ' Ta v lambda Rc u Rcrp Ls DR DCR Ntargets coords y t_y Nr']);

eval(['save ',mystery_file,' c tau W fs s Npulses PRF PRI T_out fc', ...
        ' Ta v lambda Rc u Rcrp Ls DR DCR y t_y Nr']);

fid=fopen(listing_file,'w');

fprintf(fid,['\rDESCRIPTION OF DATA IN FILE ',file,'.mat AND ',file,'_mys.mat\r\r']);
fprintf(fid,'\nRange resolution = %g meters',DR);
fprintf(fid,'\nCross-range resolution = %g meters',DCR);
fprintf(fid,'\nRange to central reference point = %g km',Rcrp/1e3);
fprintf(fid,'\nAperture time = %g seconds',Ta);
fprintf(fid,'\nRange curvature at CRP = %g meters\n\n',Rc);
fprintf(fid,'\rPulse length = %g microseconds\r',tau/1e-6);
fprintf(fid,'\rChirp bandwidth = %g Mhz\r',W/1e6);
fprintf(fid,'Time-bandwidth product = %g\n',W*tau);
fprintf(fid,'\rWe are simulating %g pulses at an RF of %g GHz',Npulses,fc/1e9);
fprintf(fid,'\r     and a PRF of %g kHz, giving a PRI of %g usec.',PRF/1e3,PRI/1e-6);
fprintf(fid,'\nRange resolution = %g meters',DR);
fprintf(fid,'\nCross-range resolution = %g meters',DCR);
fprintf(fid,'\nRange to central reference point = %g km',Rcrp/1e3);
fprintf(fid,'\nSwath length = %g m',Ls);
fprintf(fid,'\nRange curvature at CRP = %g meters\n\n',Rc);
fprintf(fid,'\rThe range window limits are %g to %g usec.\r', ...
    T_out(1)/1e-6,T_out(2)/1e-6);
fprintf(fid,'\rThe range window starts at %g km.',c*T_out(1)/2/1e3);
fprintf(fid,'\rThe range window ends at %g km.',c*T_out(2)/2/1e3);
fprintf(fid,'\nAperture time = %g seconds',Ta);
fprintf(fid,'\rThere are %g scatterers with the following coordinates:', ...
    Ntargets);
for nt = 1:Ntargets
    fprintf(fid,'\r  x=%+8.1f m, R=%+8.1f m',coords(nt,1),coords(nt,2) );
end

fclose(fid);

fprintf(['\n\nData is in file ',data_file])
fprintf(['\nStudent data is in file ',mystery_file])
fprintf(['\nListing is in file ',listing_file,'\n\n'])
