close all
clear all
kOS = 1.3;
x = git_chirp(1e6,100e-6,kOS);
X = fftshift(abs(fft(x,1024)));
f = (0:1023)/1024-0.5;

figure(1)
plot(f,X)

wx = x.*bartlett(length(x));
WX = fftshift(abs(fft(wx,1024)));
fco = 0.5/kOS;
indexu = find(f>=fco); indexl = find(f<=-fco);
indexu=indexu(1); indexl=indexl(end);
Wham = zeros(size(WX));
Wham(indexl:indexu) = bartlett(indexu-indexl+1);
Wham = Wham*max(WX);
figure(2)
plot(f,[WX,Wham])
