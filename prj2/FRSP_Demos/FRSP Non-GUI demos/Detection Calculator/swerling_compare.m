%
% M file for computing a figure to compare the 5 swerling cases + Albersheim
% for N=10 pulses, Pfa = 1e-8, and SNR from -2 to+15 dB
%
% Mark Richards
% July 2002

SNR_dB = linspace(-2,15);
Pfa = 1e-8*ones(size(SNR_dB));
N = 10*ones(size(SNR_dB));

% Step through the cases:

Pd0 = Pd(N,Pfa,SNR_dB,0); % nonfluctuating; also called Swerling 0 or 5 in some cases
Pd1 = Pd(N,Pfa,SNR_dB,1);
Pd2 = Pd(N,Pfa,SNR_dB,2);
Pd3 = Pd(N,Pfa,SNR_dB,3);
Pd4 = Pd(N,Pfa,SNR_dB,4);
Pd6 = Pd(N,Pfa,SNR_dB,6); % Albersheim's equation

% OK, now draw the results

plot(SNR_dB,[Pd0; Pd1; Pd2; Pd3; Pd4; Pd6])
axis([-2,15,0,1]);
xlabel('SNR (dB)'); ylabel('Pd'); grid;
legend('Nonfluctuating','Swerling 1','Swerling 2','Swerling 3', ...
    'Swerling 4','Albersheim','Location','Southeast');