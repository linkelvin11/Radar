%
% procSARdata_DBS_6272v3
%
%  process SAR project data using DBS algorithm
%
%  Written by M. A. Richards
%  December 2010
%

clear all; close all; hold off
format compact

%#######################################################
% User input section ###################################

% algorithm control parameters
dechirp = false;        % use azimuth dechirp step or not
oversample_freq = 5;    % oversampling in Doppler; bigger makes better
                        % picture but needs more memory and time
fix_geometry = false;   % perform geometric corrections or not

% End user input section ###############################
% ######################################################

% get data and dislay its characteristics.  Un-comment one of the two
% following lines to load the desired data set

% load dbs_data_6272_hard;
load dbs_data_6272_easy;

fprintf(['\rDESCRIPTION OF DATA\r\r']);
fprintf('\nRange resolution = %g meters',DR);
fprintf('\nCross-range resolution = %g meters',DCR);
fprintf('\nRange to central reference point = %g km',Rcrp/1e3);
fprintf('\nAperture time = %g seconds',Ta);
fprintf('\nRange curvature at CRP = %g meters\n\n',Rc);
fprintf('\rPulse length = %g microseconds\r',tau/1e-6);
fprintf('\rChirp bandwidth = %g Mhz\r',W/1e6);
fprintf('\rTime-bandwidth product = %g',W*tau);
fprintf('\rFast time sampling rate = %g Msamples/sec\r',fs/1e6);
fprintf('\rWe are simulating %g pulses at an RF of %g GHz',Npulses,fc/1e9);
fprintf('\r     and a PRF of %g kHz, giving a PRI of %g usec.',PRF/1e3,PRI/1e-6);
fprintf('\nRange resolution = %g meters',DR);
fprintf('\nCross-range resolution = %g meters',DCR);
fprintf('\nRange to central reference point = %g km',Rcrp/1e3);
fprintf('\nSwath length = %g m',Ls);
fprintf('\nRange curvature at CRP = %g meters\n\n',Rc);
fprintf('\rThe range window limits are %g to %g usec.\r', ...
    T_out(1)/1e-6,T_out(2)/1e-6);
fprintf('\rThe range window starts at %g km.',c*T_out(1)/2/1e3);
fprintf('\rThe range window ends at %g km.',c*T_out(2)/2/1e3);
fprintf('\rThere are %g scatterers with the following coordinates:', ...
    Ntargets);
for nt = 1:Ntargets
    fprintf('\r  x=%+8.1f m, R=%+8.1f m',coords(nt,1),coords(nt,2) );
end
tic

% Process the data

R = c*t_y/2; % vector of range bin centers

% Pulse compression first.  Doing this in the frequency domain. No window.
Ns = length(s);
h = conj(s(Ns:-1:1));
Nconv = Nr+Ns-1;
H = fft(h,Nconv);
fprintf('\n\nBeginning fast-time processing\n')

for m = 1:Npulses
    if (mod(m,100)==1)
        disp(['  ... compressing pulse #',int2str(m),' of ',int2str(Npulses)])
    end
    ypc = ifft(fft(y(:,m),Nconv).*H);
    y(:,m) = ypc(Ns:end);
end

figure(1)
ydisp = db(max(abs(y),eps),'voltage');  % adding machine epsilon provents a "log of zero" warning later
ydisp = max(ydisp,max(ydisp(:))-60);
imagesc(u,(R-Rcrp)/1e3,ydisp);
xlabel('aperture position (m)');
ylabel('range relative to CRP (km)');
title('After Range Compression, Before Cross-Range Compression')
grid;

% Now azimuth "compression", which is just a DFT in DBS. Initially, scale
% all range bins using the range to the CRP.  Will fix this simplification
% later. Dechirp the data first if that option selected.

% if dechirp
%     for n = 1:Nr
%         y(n,:) = y(n,:).*exp(j*(4*pi/lambda)*u.^2/2/R(n));
%     end
% end

fprintf('\n\nBeginning slow-time processing\n')
Nf = oversample_freq*Npulses;
scale_cr = Rcrp*lambda/2/v*((0:Nf-1)/Nf-0.5)*PRF;
ysar = zeros(Nr,Nf);
for n = 1:Nr
    if (mod(n,100)==1)
        disp(['  ... compressing range bin #',int2str(n),' of ',int2str(Nr)])
    end
    ysar(n,:) = fftshift(fft(y(n,:),oversample_freq*Npulses));
end

figure(2)
ydisp = db(max(abs(ysar),eps),'voltage');
ydisp = max(ydisp,max(ydisp(:))-60);
imagesc(scale_cr/1e3,(R-Rcrp)/1e3,ydisp);
xlabel('cross-range (km)');
ylabel('range relative to CRP (km)');
title('Fully Compressed Image')
grid;

if (fix_geometry)
    
    % now resample to square things up. Will work just on the linear scale
    % magnitude image now; don't need complex data anymore.
    fprintf('\n\nBeginning cross-range resampling\n')
    
    ysarmag = max(abs(ysar),eps);
    
    % find cross-range coordinates at narrowest spot (min range). We will
    % interpolate all range bins to this narrowest grid so that we will not
    % have to interpolate outside of the cross-range extent for which we
    % have data at any range.
    Rmin = R(1);
    min_scale_cr = Rmin*lambda/2/v*((0:Nf-1)/Nf-0.5)*PRF;
    
    % now loop through the range bins, interpolating the cross-range
    % magnitude data at each range bin to the same cross-range sample
    % positions given in min_scale_cr
    
    ysarmag_rs1 = zeros(size(ysarmag)); % preallocate to increase speed, save memory
    
    for n = 1:Nr
        if (mod(n,100)==1)
            disp(['  ... resampling range bin #',int2str(n),' of ',int2str(Nr)])
        end
        Rn = R(n);
        n_scale_cr = Rn*lambda/2/v*((0:Nf-1)/Nf-0.5)*PRF;
        ysarmag_rs1(n,:) = interp1(n_scale_cr,ysarmag(n,:),min_scale_cr,'linear',NaN);
        % cubic interpolation in the above line is somewhat better quality, but
        % much slower and not worth it just for display purposes
    end
    
    figure(3)
    ydisp_rs1 = db(max(ysarmag_rs1,eps),'voltage');
    ydisp_rs1 = max(ydisp_rs1,max(ydisp_rs1(:))-60);
    imagesc(min_scale_cr/1e3,(R-Rcrp)/1e3,ydisp_rs1);
    xlabel('cross-range (km)');
    ylabel('range relative to CRP (km)');
    title('Resampled Image')
    grid;
    
    % It's a little more complicated in range, since we need a
    % cross-range-dependent shift in range, not just a linear scaling of the
    % axis like we had in cross-range.
    fprintf('\n\nBeginning range shifting\n')
    
    ysarmag_rs2 = zeros(size(ysarmag_rs1)); % preallocate to increase speed, save memory
    
    % Cycle through the cross-range positions
    for m = 1:Nf
        
        if (mod(m,100)==1)
            disp(['  ... shifting cross-range bin #',int2str(m),' of ',int2str(Nf)])
            pause(0.01)
            % this 'pause' allows output to comand window to update; wasn't
            % happening otherwise for some reason, even though it worked
            % fine in the three other similar segments above without the pause.
        end
        
        % for this cross-range position, do the nonlinear range shift
        Ru = sqrt(R.^2 + min_scale_cr(m)^2);
        
        for n = 1:length(R)
            nu = find(R>=Ru(n),1,'first'); %finds closest nu such that R(nu) >= Ru(n)
            if isempty(nu)
                ysarmag_rs2(n,m) = NaN;
            else
                % This line moves into the current range bin the nearest
                % range bin in the warped image, i.e. the shift is limited
                % to integer numbers of ragne bins.  Relatively fast, but
                % crude.
                ysarmag_rs2(n,m) = ysarmag_rs1(nu,m);
                
                % This line does a linear interpolation between the range
                % bin data above and one range bin earlier, so it improves
                % the result.  However, it runs about 28x slower! If you
                % want to try this, comment out the line above and
                % uncomment this one.
%                 ysarmag_rs2(n,m) = interp1(R,ysarmag_rs1(:,m),Ru(n),'linear',NaN);

            end
        end
    end
    
    figure(4)
    ydisp_rs2 = db(max(ysarmag_rs2,eps),'voltage');
    ydisp_rs2 = max(ydisp_rs2,max(ydisp_rs2(:))-60);
    imagesc(min_scale_cr/1e3,(R-Rcrp)/1e3,ydisp_rs2);
    xlabel('cross-range (km)');
    ylabel('range relative to CRP (km)');
    title('Resampled and Range-Shifted Image')
    grid;
    
end  % end of "fix_geometry" condition

toc
