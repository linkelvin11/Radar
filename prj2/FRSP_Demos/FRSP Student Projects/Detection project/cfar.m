% cfar.m
%
% simple cell-averaging lead/lag CFAR
%
% Mark A. Richards
% April 2000
% Updated Nov. 2006

clear all, close all
load cfar_even.mat;

% Compute threshold for Pfa = 10^(-2) and mean power = 1.
% Because the mean power = 1, this is also the threshold multiplier
% for any mean power level.

Pfa = 10^(-2);
Tmult = -log(Pfa)

% Set up the ideal threshold

T_ideal = Tmult*ones(size(z));

% The CFAR has 50 each lead and lag cells, and 3 guard
% cells on each side of the cell under test.  

Nref = 50;
Nguard = 3;
Nz=length(z);

% Compute cell-averaging CFAR multiplier

Nc = 2*Nref;
alpha = Nc*((Pfa^(-1/Nc)) - 1)

% Do a brute force sliding window CFAR over the range of
% indices where the reference windows fully overlap the data.
% This is inefficient but easy to understand.
% Use NaNs ("not-a-numbers") for the undefined values in order
% to create a full length output sequence for plotting convenience.


first = Nref + Nguard + 1;
last = Nz - Nref - Nguard;
avg = NaN*ones(size(z));
for k=first:last
  avg(k) = mean([z(k-Nguard-Nref:k-Nguard-1);z(k+Nguard+1:k+Nguard+Nref)]);
end

T_cfar = alpha*avg;

% plot the results.  First, data and thresholds
figure
plot([z T_ideal T_cfar]);
xlabel('sample'); ylabel('power');grid

% now the threshold histogram
figure
hist(T_cfar,50); axis([0 20 0 1200]);
xlabel('threshold value'); ylabel('number of occurrences');
title('Histogram of CFAR Threshold Values')
% add vertical marker at the value of the ideal threshold
hold on; v = axis;
plot([Tmult,Tmult],[v(3),v(4)],'r')
hold off;

% Observed Pfa

Ncross = sum(z>T_cfar)
Pfa_obs = Ncross/(Nz-2*Nref-2*Nguard)

% Now repeat the whole business for the CFAR_uneven.mat data file

load cfar_uneven.mat;

% Set up the ideal thresholds

T_ideal1 = [Tmult*ones(Nz/2,1); NaN*ones(Nz-Nz/2,1)];
T_ideal10 = [NaN*ones(Nz/2,1); 10*Tmult*ones(Nz-Nz/2,1)];

% Do the CFAR.

avg = NaN*ones(size(z));
for k=first:last
  avg(k) = mean([z(k-Nguard-Nref:k-Nguard-1);z(k+Nguard+1:k+Nguard+Nref)]);
end

T_cfar = alpha*avg;

% plot the results.  First, data and thresholds
figure
plot([z T_ideal1 T_ideal10 T_cfar]);
xlabel('sample'); ylabel('power');grid

% re-plot, but just the transition region
figure
plot([z(2300:2700) T_ideal1(2300:2700) T_ideal10(2300:2700) T_cfar(2300:2700)]);
xlabel('sample'); ylabel('power');grid; axis([0 400 0 60])


% now the threshold histogram
figure
hist(T_cfar,100); axis([0 100 0 1200]);
xlabel('threshold value'); ylabel('number of occurrences');
title('Histogram of CFAR Threshold Values')
% add vertical marker at the values of the ideal thresholds
hold on; v=axis;
plot([Tmult,Tmult],[v(3),v(4)],'r')
plot(10*[Tmult,Tmult],[v(3),v(4)],'r')
hold off;



% Observed Pfa

Ncross = sum(z>T_cfar)
Pfa_obs = Ncross/(Nz-2*Nref-2*Nguard)


% Now repeat for the target case

load cfar_target.mat
Nz = length(z);

% Standard CA CFAR
first = Nref + Nguard + 1;
last = Nz - Nref - Nguard;
avg_ca = NaN*ones(size(z));
for k=first:last
  avg_ca(k) = mean([z(k-Nguard-Nref:k-Nguard-1);z(k+Nguard+1:k+Nguard+Nref)]);
end

% Repeat for SOCA CFAR
avg_soca = NaN*ones(size(z));
for k=first:last
   avg_soca(k) = min((sum(z(k-Nguard-Nref:k-Nguard-1))), ...
      sum(z(k+Nguard+1:k+Nguard+Nref)))/Nref;
end

% Repeat for GOCA CFAR
avg_goca = NaN*ones(size(z));
for k=first:last
   avg_goca(k) = max((sum(z(k-Nguard-Nref:k-Nguard-1))), ...
      sum(z(k+Nguard+1:k+Nguard+Nref)))/Nref;
end

alpha_soca = 5.175;
alpha_goca = 6.686;

T_ca = alpha*avg_ca;
T_soca = alpha_soca*avg_soca;
T_goca = alpha_goca*avg_goca;

T_ideal = Tmult*ones(size(z));

% plot the results.  First, data and thresholds
figure
plot(db([z T_ideal T_ca  T_soca  T_goca],'power'));
axis([450 650 0 25]);
xlabel('sample'); ylabel('power');grid
