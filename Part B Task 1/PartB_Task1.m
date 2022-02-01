clc;
clear;
close all;

load("Rx1.mat");
load("Rx2.mat");
load("Rx3.mat");
load("Rx4.mat");

% Cartesian coordinates of four RX
r1 = [0;0;0];
r2 = [60;-88;0];
r3 = [100;9;0];
r4 = [60;92;0];
rs = [r1, r2, r3, r4];
numOfRx = 4;

% Carrier Frequency
Fc = 2.4e9;
% Light velocity
c = 3e8;
% Symbol Duration
Tcs = 5e-9;
% Path loss exponent
alpha = 2;
% Number of Rx
N = 4;
% SNR in dB
SNR = 20;
snr = 10^(SNR/10);
% Noise power in dB
noiseVar = 5;
% Sampling period
Ts = 5e-9;

%% TOA
% Association Stage
figure();
plot(1:length(x1_Time), abs(x1_Time), 1:length(x1_Time), abs(x2_Time), ...
    1:length(x1_Time), abs(x3_Time), 1:length(x1_Time), abs(x4_Time));
title("Magnitude of signals received by four antennas");
xlabel("Time");
ylabel("Signal Magnitude");
legend("x1", "x2", "x3", "x4");

% % set one third of maximum magnitude as the threshold
% th1 = max(abs(x1_Time)) / 3;
% th2 = max(abs(x2_Time)) / 3;
% th3 = max(abs(x3_Time)) / 3;
% th4 = max(abs(x4_Time)) / 3;

% set the 60% of the useful signal magnitude as the threshold to detect the
% momenet when the transmitted signal is received.
noiseVar = 10 ^ (noiseVar/10);
signalPower = noiseVar * snr;
th = sqrt(signalPower) * 0.6;

t0 = 20 * Ts;
t1 = find(abs(x1_Time) >= th, 1) * Ts;
t2 = find(abs(x2_Time) >= th, 1) * Ts;
t3 = find(abs(x3_Time) >= th, 1) * Ts;
t4 = find(abs(x4_Time) >= th, 1) * Ts;
ts = [t1, t2, t3, t4];

rho1 = (t1-t0)*c;
rho2 = (t2-t0)*c;
rho3 = (t3-t0)*c;
rho4 = (t4-t0)*c;
rhos = [rho1, rho2, rho3, rho4];

% metric fusion stage
H_toa = [r2 r3 r4].';
b_toa = 0.5 * [(norm(r2))^2 - rho2^2 + rho1^2; ...
               (norm(r3))^2 - rho3^2 + rho1^2;
               (norm(r4))^2 - rho4^2 + rho1^2;];
rm_toa = pinv(H_toa) * b_toa;

% Plot positioning circles
figure();
hold on;
circlePlot(r1, rho1);   text(r1(1),r1(2),' r1');
circlePlot(r2, rho2);   text(r2(1),r2(2),' r2');
circlePlot(r3, rho3);   text(r3(1),r3(2),' r3');
circlePlot(r4, rho4);   text(r4(1),r4(2),' r4');
plot(rm_toa(1),rm_toa(2),'sk','MarkerFaceColor','r'); text(rm_toa(1),rm_toa(2),' Tx');
hold off;
daspect([1,1,1]);
xlabel("x");    ylabel("y");
title(sprintf("TOA results: r_m = [%.2f, %.2f, %.2f]^T", rm_toa(1), rm_toa(2), rm_toa(3)));


%% TDOA
% TDOA method does not suffer from MS clock synchronization errors
% Association Stage has been done in last TOA Section
% metric fusion stage
H_tdoa = rs(:, 2:1:numOfRx).';
rhoi1s = zeros(1, numOfRx);
b_tdoa = [];

% Estimate the unknown parameter rho1
syms rho1_tdoa;

for i = 2:1:numOfRx
    ti = ts(i);
    rhoi = rhos(i);
    rhoi1 = (ti-ts(1)) * c;
    rhoi1s(i) = rhoi1;
    b_tdoa = [b_tdoa; (norm(rs(:,i)))^2 - rhoi1^2 - 2*rhoi1 * rho1_tdoa];
end
b_tdoa = b_tdoa / 2;

eqn = (pinv(H_tdoa'*H_tdoa) * H_tdoa'*b_tdoa)'*pinv(H_tdoa'*H_tdoa)*H_tdoa'*b_tdoa - rho1_tdoa^2 == 0;

rho1_tdoa_roots = double(solve(eqn, rho1_tdoa));
% use the positive root as the result
rho1_tdoa_sol = rho1_tdoa_roots(rho1_tdoa_roots>0);

b_tdoa = subs(b_tdoa,rho1_tdoa,rho1_tdoa_sol);
b_tdoa = double(b_tdoa);

rm_tdoa = pinv(H_tdoa) * b_tdoa;

% Plot hyperbolic curves

figure();
circlePlot(r1, rho1_tdoa_sol);      text(r1(1),r1(2),'  r1');
title(sprintf("TDOA results : r_m = [%.2f, %.2f, %.2f]^T", rm_tdoa(1), rm_tdoa(2), rm_tdoa(3)));
daspect([1,1,1]);
hold on;

hyperbolicPlot(r1, r2, rhoi1s(2));  text(r2(1),r2(2),'  r2');
hyperbolicPlot(r1, r3, rhoi1s(3));  text(r3(1),r3(2),'  r3');
hyperbolicPlot(r1, r4, rhoi1s(4));  text(r4(1),r4(2),'  r4');
plot(rm_tdoa(1),rm_tdoa(2),'sk','MarkerFaceColor','r'); text(rm_toa(1),rm_toa(2),'  Tx');
xlabel("x");    ylabel("y");

hold off;

