function AF=ambiguity(x,N,M,Tf,labels)
%
% ambiguity
%
% calling sequence: AF=ambiguity(x,N,M,Tf,labels)
%
% M-file to compute the ambiguity function of a real or
% complex vector
%
%  x = complex input vector, the signal to be analyzed
%  N = length of x
%  M = fft size.  M >= 2*N-1.
% Tf = fast time sampling interval in seconds
% labels = optional argument to control axis labels.
%          If labels=='normalize', axes are plotted relative to pulse
%          length.  If labels is any other value, absolute (seconds and Hz)
%          labels are used.
%
% Mark Richards
% March 2002
% Revised August 2004
%

% ensure that data is used in row vector form.  Be careful
% NOT to conjugate it at the same time it is transposed
% (MATLAB's transpose is a Hermitian)

[r c] = size(x);
if (c == 1)
  x = x.';
end

% ensure adequate FFT size is used

if M < 2*N-1
   disp('Error in ambiguity.m: FFT size too small. Must have M >= 2*N-1.')
   return
end

% determine whether to use absolute or relative labeling

abs_label = 0;
if (nargin < 5)
    abs_label = true;
else
    abs_label = ~strcmp(labels,'normalize');
end

% form matrix of overlapped signals.  Deliberately include one sample in
% each direction in time that will be zero; makes for better plots

AF = zeros(2*N+1,M);
for n=-N+2:N
  ll = max(1,n);
  ul = min(N,n+N);
  if (n <= 1)
      AF(n+N,1:n+N-1) = x(1:n+N-1).*conj(x(2-n:N));
  else
      AF(n+N,n:N) = x(n:N).*conj(x(1:N-n+1));
  end
end

% now transform each column separately, and shift to get
% frequency origin in the middle of the array

for n = 1:2*N+1
 AF(n,:) = fftshift(M*ifft(AF(n,:)));
end

% ambiguity function is magnitude (not squared) of result;

AF = abs(AF);
peak=max(max(AF));

% set up normalized and unnormalized time and frequency scales

f_abs = ((-M/2:M/2-1)/M)*(1/Tf);
t_abs = (-N:N)*Tf;

f_rel = ((-M/2:M/2-1)/M)*(1/Tf)*(N*Tf);
t_rel = (-N:N)/N;

%  do plots

figure(1);
if (abs_label)
    surf(f_abs,t_abs,AF);
    axis([min(f_abs),max(f_abs),min(t_abs),max(t_abs),0,max(max(AF))])
    xlabel('Doppler shift (Hz)')
	ylabel('delay (sec)')
else
    surf(f_rel,t_rel,AF);
    axis([min(f_rel),max(f_rel),-1,1,0,max(max(AF))])
    xlabel('normalized Doppler shift')
    ylabel('normalized delay')
end
colormap('hsv')
title('Ambiguity Surface')

figure(2);
axis('square')
c=[0.994*peak,0.708*peak,0.316*peak,0.1*peak];
if (abs_label)
	contour(f_abs,t_abs,AF,c);
    axis([min(f_abs),max(f_abs),min(t_abs),max(t_abs)])
	xlabel('Doppler shift (Hz)')
	ylabel('delay (sec)')
else
	contour(f_rel,t_rel,AF,c);
    axis([min(f_rel),max(f_rel),-1,1,])
    xlabel('normalized Doppler shift')
    ylabel('normalized delay')
end
title('-0.5, -3, -10, and -20 dB Contours')

figure(3);
axis('normal')
if (abs_label)
	plot(t_abs,AF(:,M/2+1))
    axis([min(t_abs),max(t_abs),0,max(AF(:,M/2+1))])
	xlabel('delay (sec)')
	ylabel('amplitude')
else
	plot(t_rel,AF(:,M/2+1))
    axis([-1,1,0,max(AF(:,M/2+1))])
	xlabel('normalized delay')
	ylabel('amplitude')
end
title('Zero-Doppler Cut')

figure(4);
if (abs_label)
	plot(f_abs,AF(N,:))
    axis([min(f_abs),max(f_abs),0,max(AF(N,:))])
	xlabel('Doppler shift (Hz)')
	ylabel('amplitude')
else
	plot(f_rel,AF(N,:))
    axis([min(f_rel),max(f_rel),0,max(AF(N,:))])
	xlabel('normalized Doppler shift')
	ylabel('amplitude')
end
title('Zero-Delay Cut')
