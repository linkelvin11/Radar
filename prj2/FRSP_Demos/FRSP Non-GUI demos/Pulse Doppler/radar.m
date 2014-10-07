function y = radar( x, fs, T_0, g, T_out, T_ref, fc, r, snr, v )
% RADAR      simulate radar returns from a single pulse or burst
%            of identical pulses
%  usage:
%    R = radar( X, Fs, T_0, G, T_out, T_ref, Fc, R, SNR, V )
%      X:      baseband single pulse waveform (complex vector)
%      Fs:     sampling frequency of input pulse        [in Hz]
%      T_0:    start time(s) of input pulse(s)          [sec]
%                (number of pulses in burst assumed = length(g) )
%      G:      complex gain(s) of pulse(s)
%      T_out:  2-vector [T_min,T_max] defines output
%                window delay times w.r.t. start of pulse
%      T_ref:  system "reference" time, needed to simulate
%                 burst returns. THIS IS THE "t=0" TIME !!!
%      Fc:     center freq. of the radar.                [in Hz]
%      R:      vector of ranges to target(s)             [meters]
%                 (number of targets assumed = length(r) )
%    SNR:      vector of target SNRs (unit noise power assumed)
%      V:      vector of target velocities (optional)    [in m/sec]
%                 (positive velocities are towards the radar)
%
%  note(1): VELOCITY in meters/sec !!!
%           distances in m, times in sec, BW in Hz.

%  note(2): assumes each pulse is constant (complex) amplitude
%  note(3): will accomodate up to quadratic phase pulses
%  note(4): vector of ranges, R, allows DISTRIBUTED targets
%
%     (c) jMcClellan 7/28/90
%         Modified by M. A. Richards, August 1991

J = sqrt(-1);
c = 3e8;                               % velocity of light in m/sec
Mx = length(x);
delta_t = 1/fs;                        % sampling interval (sec)
t_y = [ T_out(1):delta_t:T_out(2) ]';  % output sampling times (sec)
T_p = Mx*delta_t;                      % length of input pulse (sec)

% Assume zero velocities (stationary targets) if no velocity
% vector provided

if nargin < 7
  v = zeros(r);
end

% ensure that all vectors are column vectors

x=x(:);  g=g(:);  T_0=T_0(:);  r=r(:); snr=snr(:);  v=v(:);

% determine the quadratic phase modulation parameters for
% later interpolation of pulse samples

t_x = delta_t*[0:(Mx-1)]';
x_ph = unwrap(angle(x));
q = polyfit(t_x,x_ph,2);

% check result using correlation coefficient

xfit = polyval(q,t_x);
if (x_ph'*xfit)/norm(x_ph)/norm(xfit) < 0.99
  disp('pulse phase is not quadratic')
  keyboard
end

%
%---  Form (initially empty) output matrix ---
%
Mr = length(t_y);  Nr = length(g);     % output samples in a matrix
y = zeros(Mr,Nr);

% Index 'i' loops over the number of targets
for  i = 1:length(r)
   ri = r(i);
   vi = v(i);
   f_doppler = 2*vi*fc/c;

% Index 'j' loops over the number of pulses
 for j = 1:length(g)
   r_at_T_0 = ri - vi*T_0(j);

% Compute start and end time of reflected pulse at receiver,
% ensure that it falls at least partially within the range (time) window
   tau = 2*r_at_T_0/(c+vi);   tmax = tau + T_p;
   if tau >= T_out(2) | tmax <= T_out(1)
      fprintf('\nEcho from target #%g at range %g km',i,ri)
      fprintf('\nis COMPLETELY OUT OF the range window')
      fprintf('\non pulse #%g.\n',j)
   else

% Figure out which sample locations in the output grid contain
% reflected pulse
      t_vals = t_y - tau;
      n_out = find( t_vals >= 0  &  t_vals < T_p );
      if tau < T_out(1)
        fprintf('\nEcho from target #%g at range %g km',i,ri)
        fprintf('\nSTARTS BEFORE the range window')  
        fprintf('\non pulse #%g.\n',j)
      end
      if tmax > T_out(2)
        fprintf('\nEcho from target #%g at range %g km',i,ri)
        fprintf('\nFINISHES AFTER the range window')  
        fprintf('\non pulse #%g.\n',j)
      end

% Place scaled, range-delayed, Doppler shifted pulse into output matrix
% Unit noise power and unit nominal pulse amplitude assumed to
% get amplitude from SNR.  Amplitude also adjusted for relative R^2 loss,
% which can be done equivalently as t^2.  this means the target has the
% specified SNR only if it's in the first range bin, otherwise it is less
% by the R^4 loss.

  amp = 10^(snr(i)/20);  amp = amp*((T_out(1)/tau)^2);
  y(n_out,j) = y(n_out,j) + ...
          ( amp * g(j) * exp( -J*2*pi*fc*tau ) ) ...
          .* [ exp( J*2*pi*f_doppler*t_y(n_out) ) ]  ...
          .* [ exp( J*polyval(q,t_vals(n_out)) ) ];

    end  % end of "if tau >= T_out(2) ...."
  end    % end of loop over j
end      % end of loop over i
