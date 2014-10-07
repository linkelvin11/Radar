%
% Example M file for computing a figure to illustrate the required SNR for
% a given Pd, Pfa, and N for each of the Swerling models.
%
% Mark Richards June 2010

clear all
close all

% fix Pfa
pfa = 1e-8;

% set up a grid of pd and N values
% pd = [(.1:.05:.75),(.8:.02:.88),(.9:.01:.99)];
pd = [(0.1:0.1:0.9),0.95,0.975,0.99];
% n = [1:9,10:2:18,20:5:50];
n = 1:50;

[PD,N] = meshgrid(pd,n);

% fix Pfa in this example
PFA = 1e-8*ones(size(PD));

% Do it for all the Swerling cases
S0 = SNRdB(PD,PFA,N,0); xlabel('Pd'), ylabel('N');
S1 = SNRdB(PD,PFA,N,1); xlabel('Pd'), ylabel('N');
S2 = SNRdB(PD,PFA,N,2); xlabel('Pd'), ylabel('N');
S3 = SNRdB(PD,PFA,N,3); xlabel('Pd'), ylabel('N');
S4 = SNRdB(PD,PFA,N,4); xlabel('Pd'), ylabel('N');
S5 = SNRdB(PD,PFA,N,5); xlabel('Pd'), ylabel('N');
SA = SNRdB(PD,PFA,N,6); xlabel('Pd'), ylabel('N');

% Plot difference in SNR vs. nonflucutating case
figure; mesh(pd,n,S0); title('Single-Sample SNR (dB), Marcum')
xlabel('Pd'), ylabel('N')
figure; mesh(pd,n,S0-S1); title('S1 diff'); xlabel('Pd'), ylabel('N');
figure; mesh(pd,n,S0-S2); title('S2 diff'); xlabel('Pd'), ylabel('N');
figure; mesh(pd,n,S0-S3); title('S3 diff'); xlabel('Pd'), ylabel('N');
figure; mesh(pd,n,S0-S4); title('S4 diff'); xlabel('Pd'), ylabel('N');
figure; mesh(pd,n,S0-S5); title('S5 diff'); xlabel('Pd'), ylabel('N');
figure; mesh(pd,n,S0-SA); title('Alb diff'); xlabel('Pd'), ylabel('N');

        