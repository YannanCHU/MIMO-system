clc;
clear;
close all;

load yc3021.mat

% Task A
phi = phi_mod;
phi = phi * pi / 180;
% phi = phi_mod;

% PN-codes
% two primitive polynomials that are used for gold-sequence production
% D^5 + D^2 + 1
% D^5 + D^3 + D^2 + D + 1
c1 = [1; 0; 0; 1; 0; 1];
c2 = [1; 0; 1; 1; 1; 1];

% two m-sequences
m1 = fMSeqGen(c1);
m2 = fMSeqGen(c2);
m = size(c1,1)-1;
Nc = 2 ^ m - 1;

% desired user's delay for gold-sequence:
d = phase_shift;

% gold-sequence generation
goldseq = fGoldSeq(m1, m2, d);

% number of path of each source.
paths = length(Beta_1);

% complex fading coefficient of each path
betas = Beta_1;

% Uniform Circular Array - hald-wavelength inter-antenna spacing
numOfAe = 5;
angleDiff = 360 / numOfAe;
radius = 0.5 / sin(deg2rad(angleDiff*2));
array = zeros(numOfAe,3);
for Ae_ith = 1:1:numOfAe
    array(Ae_ith, 1:1:2) = radius * [cos(deg2rad(30+angleDiff*(Ae_ith-1))) sin(deg2rad(30+(Ae_ith-1)*angleDiff))];
end

%% channel estimation
symbolsIn = Xmatrix.';
[delay_estimate, DOA_estimate] = fChannelEstimation(symbolsIn, goldseq, paths, array);

disp("the estimated delay is " + sprintf("%d, %d, and %d", delay_estimate));
disp("The estimated DOAs are " + sprintf("[%d, %d]", DOA_estimate(1,1), DOA_estimate(1,2))+ ... 
    sprintf("[%d, %d]", DOA_estimate(2,1), DOA_estimate(2,2)) + ... 
    sprintf("[%d, %d]", DOA_estimate(3,1), DOA_estimate(3,2)));

%% Spatio-Temporal Raker Beamformer
weights = STbeamformer(array, delay_estimate, DOA_estimate, betas, goldseq, paths);

%% Symbol demodulator
symbolsInMatrix = symbolExtension(symbolsIn, size(array, 1), length(goldseq), 1);
% Combine the signals received by different antennas
symbolsIn = (weights' * symbolsInMatrix).';
% QPSK Demodulation
bitsOut=fDSQPSKDemodulator(symbolsIn, phi);

%% Reconstruct and display the message
charBitNum = 8;
messageBits = uint8(reshape(bitsOut, charBitNum, length(bitsOut)/charBitNum));
messageChars = char(bi2de(messageBits', "left-msb")');
disp(messageChars);