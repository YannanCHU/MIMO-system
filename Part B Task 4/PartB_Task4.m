clc;
clear;
close all;

load("Xmatrix_LAA_1.mat");
load("Xmatrix_LAA_2.mat");
load("Xmatrix_LAA_3.mat");
load("Xmatrix_LAA_4.mat");
Xs_LAA = {x1_LAA; x2_LAA; x3_LAA; x4_LAA};

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
% wavelength
lambda = c / Fc;


% store useful signal eigen values without noise
lambdas = zeros(N, 1);
% store the ratio of different ranges and the reference range (i.e., rho1)
K = zeros(N-1, 1);

for i = 1:1:N
    % calculate the covariance matrix
    Rxx_ith = Xs_LAA{i} * Xs_LAA{i}' / size(Xs_LAA{i},2);
    
    % sort the eigen values of covariance matrix
    eigVals_ith = sort(eig(Rxx_ith), 'descend');
    
    % useful signal eigen values with noise
    gamma_ith = eigVals_ith(1);
    
    % noise eigen values
    noiseVar_ith = mean(eigVals_ith(2:end));
    
    % useful signal eigen values without noise by
    % subtracting from gamma the noise estimate
    lambdas(i,1) = gamma_ith - noiseVar_ith;
    
    if i >= 2
        K(i-1,1) = (lambdas(i,1) / lambdas(1,1)) ^ (1/(2*alpha));
    end
end

% Metric fusion stage - r is 3 by N-1 matrix
r = rs(:,2:end);
rx = r(1,:)';
ry = r(2,:)';
rz = r(3,:)';

H = [2 * (ones(N-1, 1) * r1' - r'), (ones(N-1,1) - K.^2)];
b = [(norm(r1))^2 * ones(N-1,1) - rx.^2 - ry.^2 - rz.^2];

rm_rho1 = pinv(H) * b;
rm_laa = rm_rho1(1:end-1);
rho1_sq = rm_rho1(end);
rho1 = sqrt(rho1_sq);

rhos = [rho1; zeros(N-1,1)];
for i = 2:1:N
    rhos(i) = K(i-1) * rho1;
end

% Plot positioning circles
figure();
hold on;
for i = 1:1:N
    txt = sprintf(" r%d", i);
    circlePlot(rs(:,i), rhos(i));   text(rs(1,i),rs(2,i),txt);
end

plot(rm_laa(1),rm_laa(2),'sk','MarkerFaceColor','r'); text(rm_laa(1),rm_laa(2),' Tx');
hold off;
daspect([1,1,1]);
xlabel("x");    ylabel("y");
title(sprintf("LAA results: r_m = [%.2f, %.2f, %.2f]^T", rm_laa(1), rm_laa(2), rm_laa(3)));
