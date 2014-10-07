%
% M file for computing a figure to illustrate the effect of
% noncoherent integration on Pd vs. SNR for Swerling 0,
% for Pfa = 1e-8, SNR from -2 to+15 dB, and N=1 to 100.
%
% Mark Richards
% July 2002

SNR_dB = linspace(0,15);  % set up a range of SNRs to consider
Pfa = 1e-8*ones(size(SNR_dB));  % Fixed value of Pfa

% Step through the cases:

Pd0 = Pd(1*ones(size(SNR_dB)),Pfa,SNR_dB,0);
Pd1 = Pd(2*ones(size(SNR_dB)),Pfa,SNR_dB,0);
Pd2 = Pd(5*ones(size(SNR_dB)),Pfa,SNR_dB,0);
Pd3 = Pd(10*ones(size(SNR_dB)),Pfa,SNR_dB,0);
Pd4 = Pd(20*ones(size(SNR_dB)),Pfa,SNR_dB,0);

% OK, now draw the results

plot(SNR_dB,[Pd0; Pd1; Pd2; Pd3; Pd4])
axis([0,15,0,1]);
xlabel('SNR (dB)'); ylabel('Pd'); grid;
legend('N=1','N=2','N=5','N=10','N=20','Location','SouthEast');