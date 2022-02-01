clc;
clear;
close all;

load("Rx1.mat");
load("Rx2.mat");
load("Rx3.mat");
load("Rx4.mat");
Xs_RSS = {x1_RSS; x2_RSS; x3_RSS; x4_RSS};

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

% Association
% Power of transmitted signal
Ptx_dbm = 150;
Ptx_w = (1e-3) * 10^(Ptx_dbm/10);

% wavelength
lambda = c / Fc;

% Tx power gain = 1 due to isotropic antenna
Gtx = 1;

% Rx power gain = 1 due to isotropic antenna
Grx = 1;

% range estimate
rhos = zeros(N, 1);
for i = 1:1:N
    % Power of received signal
    Prxi_w = mean(abs(Xs_RSS{i,1}).^2);
    % Estimate the range
    rhos(i) = sqrt((Ptx_w/Prxi_w)*Gtx*Grx)*lambda/(4*pi);
end

% Metric Fusion Stage
H = rs(:, 2:1:end).';

b = zeros(N-1, 1);
for i = 2:1:N
    b(i-1,1) = 0.5 * ( (norm(rs(:,i)))^2 - rhos(i)^2 + rhos(1)^2 );
end

rm_rssi = pinv(H) * b;

% Plot positioning circles
figure();
hold on;
for i = 1:1:N
    txt = sprintf(" r%d", i);
    circlePlot(rs(:,i), rhos(i));   text(rs(1,i),rs(2,i),txt);
end

plot(rm_rssi(1),rm_rssi(2),'sk','MarkerFaceColor','r'); text(rm_rssi(1),rm_rssi(2),' Tx');
hold off;
daspect([1,1,1]);
xlabel("x");    ylabel("y");
title(sprintf("RSSI results: r_m = [%.2f, %.2f, %.2f]^T", rm_rssi(1), rm_rssi(2), rm_rssi(3)));


