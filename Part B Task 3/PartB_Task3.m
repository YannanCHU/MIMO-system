clc;
clear;
close all;

load("Xmatrix_1_DFarray.mat");
load("Xmatrix_2_DFarray.mat");
load("Xmatrix_3_DFarray.mat");
load("Xmatrix_4_DFarray.mat");
Xs_DOA = {x1_DOA; x2_DOA; x3_DOA; x4_DOA};

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

% Antenna Array
ruca = [0.1250  0.0625      -0.0625     -0.1250     -0.0625     0.0625; ...
    0       0.1083      0.1083      0           -0.1083     -0.1083; ...
    0       0           0           0           0           0];

% get the geometry distribution of array and convert the unit from meter to
% half wavelength
array =  ruca' / (lambda / 2);

% MUSIC algorithm to estimate DOAs
azimuth_angles = zeros(4,1);
for i = 1:1:N
    azimuth_angles(i,1) = music(array, Xs_DOA{i}, 1);
end

% Phi1 + Phi2 + Phi3 + Phi4 = 360 degree
Phi1 = anglePhiCal(r4, r1, r2);
Phi2 = anglePhiCal(r1, r2, r3);
Phi3 = anglePhiCal(r2, r3, r4);
Phi4 = anglePhiCal(r3, r4, r1);

theta1 = angleThetaCal(r2, r1, azimuth_angles(1));
theta2 = angleThetaCal(r3, r2, azimuth_angles(2));
theta3 = angleThetaCal(r4, r3, azimuth_angles(3));
theta4 = angleThetaCal(r1, r4, azimuth_angles(4));
thetas = [theta1, theta2, theta3, theta4];

rho12 = norm(r1 - r2);
rho23 = norm(r2 - r3);
rho34 = norm(r3 - r4);
rho41 = norm(r4 - r1);

matrix = [cos(theta1 * pi / 180), cos((Phi2-theta2) * pi / 180), 0, 0;
    0, cos(theta2 * pi / 180), cos((Phi3-theta3) * pi / 180), 0;
    0, 0, cos(theta3 * pi / 180), cos((Phi4-theta4) * pi / 180);
    cos((Phi1-theta1) * pi / 180), 0, 0, cos(theta4 * pi / 180)];

rhos = matrix \ [rho12;rho23;rho34;rho41];

rho1 = rhos(1);
rho2 = rhos(2);
rho3 = rhos(3);
rho4 = rhos(4);

H = kron(ones(numOfRx,1), eye(3));
b = zeros(3*numOfRx,1);
for i =1:1:numOfRx
    b(3*i-2:3*i,1) = [rs(1,i)+rhos(i)*cos(azimuth_angles(i) * pi / 180); rs(2,i)+rhos(i)*sin(azimuth_angles(i) * pi / 180); 0];
end

rm_doa = H \ b;

% Plot lines
figure();
subplot(1,2,1);
hold on;
for i = 1:1:N
    txt = sprintf(" r%d", i);
    linePlot(rs(:,i), azimuth_angles(i));   text(rs(1,i),rs(2,i),txt);
end

plot(rm_doa(1),rm_doa(2),'sk','MarkerFaceColor','r'); text(rm_doa(1),rm_doa(2),' Tx');
hold off;
xlim([min(rs(1,:))-10, max(rs(1,:))+10]);
ylim([min(rs(2,:))-10, max(rs(2,:))+10]);
daspect([1,1,1]);
xlabel("x");    ylabel("y");
title(sprintf("DOA results: r_m = [%.2f, %.2f, %.2f]^T", rm_doa(1), rm_doa(2), rm_doa(3)));

subplot(1,2,2);
hold on;
for i = 1:1:N
    txt = sprintf(" r%d", i);
    circlePlot(rs(:,i), rhos(i));   text(rs(1,i),rs(2,i),txt);
end

plot(rm_doa(1),rm_doa(2),'sk','MarkerFaceColor','r'); text(rm_doa(1),rm_doa(2),' Tx');
hold off;
daspect([1,1,1]);
xlabel("x");    ylabel("y");
title(sprintf("DOA results: r_m = [%.2f, %.2f, %.2f]^T", rm_doa(1), rm_doa(2), rm_doa(3)));

function linePlot(pos, azimuth)
    % pos (3x1 vector) = position of Rx
    % azimuth (float) = the azimuth of estimated DOA

    temp = -200:1:200;
    x = temp * cos(azimuth*pi/180) + pos(1);
    y = temp * sin(azimuth*pi/180) + pos(2);
    plot(x, y, pos(1), pos(2),'ok','MarkerFaceColor','y');
end

function theta_mid = angleThetaCal(right, mid, azimuth_angle)
    % Inputs:
    % right (3x1 vector) = the coordinate of Rx near mid
    % mid (3x1 vector) = the coordinate of target Rx
    % azimuth_angle = the azimuth angle of DOA estimated by mid
    % Outputs:
    % theta_mid = the angle in degree formed by mid and its adjacnt Rx and
    % the Tx.
    theta_comp = angle( (right(1)-mid(1)) + 1i*(right(2)-mid(2)) );
    theta_comp = theta_comp * 180 / pi;
    theta_mid = azimuth_angle - theta_comp;
    if theta_mid > 360
        theta_mid = theta_mid - 360;
    end
end

function Phi_mid = anglePhiCal(left, mid, right)
    % Inputs:
    % left, right (3x1 vector) = the coordinate of two Rxs near mid
    % mid (3x1 vector) = we are going to calculate the its angle phi formed
    % with two adjacnt Rxs
    % Outputs:
    % Phi_mid = the angle phi in degree formed by mid and its two adjacnt Rxs
    Phi_mid = angle( (left(1)-mid(1)) + 1i*(left(2)-mid(2)) ) - angle( (right(1)-mid(1)) + 1i*(right(2)-mid(2)) );
    Phi_mid = abs(Phi_mid);
    if Phi_mid > pi
        Phi_mid = 2*pi-Phi_mid;
    end
    Phi_mid = Phi_mid * 180 / pi;
end

function azimuth_angles=music(array, x_DOA, M)
% Inputs
% array (N*3 vector) = the array of antennas in receiver. N is the number
% of Rx antennas.
% x_DOA (N*L) = signals received by different antennas in the certain Rx
% M = number of source - 1
% Outputs
% azimuth_angles = the azimuth angle of DOA. The elevation is always 0.
L = size(x_DOA,2);             % number of snapshots
Rxx = (x_DOA * x_DOA') / L;   % covariance matrix

% Step 0: Assume M and array geometry are known
% Step 1: receive signal vector - N by L vector
% Step 2: find the covariance matrix  - Rxx
% Step 3: find the eigenvalue and eigenvector of Rxx
[eigenVec,eigenVal] = eig(Rxx);
[~, pos] = sort(abs(diag(eigenVal)));
% Step 4: form the matrix Es with columns the eigevectors which correspond to
% the M largest eigenvalues
signalEigVec = eigenVec(:, pos(end-M+1: end));
% Step 5: find the arg of the M minima of the function
N=10;   % adjust the resolution 
azimuths = 0:1/N:360;
elevation = 0;
cost_fun = zeros(size(azimuths));
Pn = eye(size(signalEigVec*signalEigVec')) - signalEigVec*signalEigVec';

for i = 1:1:length(azimuths)
    S = spv(array, [azimuths(i), elevation]);
    cost_fun(i) = abs(S' * Pn * S);
end

[~, azimuth_angles] = mink(cost_fun, M);
azimuth_angles = (azimuth_angles-1) / N;

end
