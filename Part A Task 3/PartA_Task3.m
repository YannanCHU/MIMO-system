clc;
clear;
close all;

% Task A
% X = 3 (C), Y = 25 (Y)
X = 3; Y = 25;
% initisl angle: phi = X + 2Y
phi = X + 2 * Y;
phi = phi * pi / 180;
% QPSK Modulation
% x_div_demod = pskdemod(s_div_sc, 4, phi);

% PN-codes
% two primitive polynomials that are used for gold-sequence production
% D^4 + D + 1
% D^4 + D^3 + 1
c1 = [1; 0; 0; 1; 1];
c2 = [1; 1; 0; 0; 1];

% two m-sequences
m1 = fMSeqGen(c1);
m2 = fMSeqGen(c2);
m = size(c1,1)-1;
Nc = 2 ^ m - 1;

% desired user's delay for gold-sequence: d >= 1 + (X+Y) mod 12
% other two delaies: d2 = d+1, d3 = d+2
d = 1 + mod((X+Y), 12);


% gold-sequence generation
% the loop is used to find the balanced gold-sequence
for k = d:1:Nc
    a1 = fGoldSeq(m1, m2, k);
%     a1 = 1 - 2 * a1;
    a2 = fGoldSeq(m1, m2, k+1);
%     a2 = 1 - 2 * a2;
    a3 = fGoldSeq(m1, m2, k+2);
%     a3 = 1 - 2 * a3;
    if sum(1 - 2 * a1,1) == -1
%         disp("Delay is " + k);
        break;
    end
end
% as stores three gold-sequences
as = [a1, a2, a3];

% the 3-layer RGB image is used and its size is limited within 160 x 112
% each pixel is represented by 8 bits (intensity ranges from 0 to 255)
bitLengthMax = 160 * 112 * 3 * 8;
sumbolNumber = bitLengthMax / 2;
R = Nc * sumbolNumber;  % Number of channel symbol chips after DS-QPSK Modulation

bitsOuts = zeros(bitLengthMax, 3);
xs = zeros(3,1);
ys = zeros(3,1);
imageBitLength = zeros(3,1);
symbolsOut = zeros(R, 3);

% Display the all three images
figure(1);
for i = 1:1:3
    filename = ['photo-', num2str(i), '.jpg'];
    img_ith = imread(filename); subplot(1, 3, i); imshow(img_ith); 
    if i == 1
        title("Desired image (" + size(img_ith, 1) + " by " + size(img_ith, 2) + " pixels)");
    else
        title("Multiple-Access Interference (" + size(img_ith, 1) + " by " + size(img_ith, 2) + " pixels)");
    end
end

for i = 1:1:3
    filename = ['photo-', num2str(i), '.jpg'];
    [bitsOuts(:,i),xs(i,1),ys(i,1)] = fImageSource(filename, bitLengthMax);
    imageBitLength(i,1) = xs(i,1) * ys(i,1) * 3 * 8;
%     fImageSink(bitsOuts(:,i),imageBitLength(i,1),xs(i,1),ys(i,1));
    goldseq = as(:,i);
    symbolsOut(:,i) = fDSQPSKModulator(bitsOuts(:,i), goldseq, phi);
end

% channel parameters
% number of path of each source.
paths = [1; 1; 1];
% delay of each path
delays = {5; 7; 12};
% complex fading coefficient of each path
betas = {0.4; 0.7; 0.2};
% Due to SISO working mode, this DOA information is meaningless actually
DOAs = {[30, 0], [90 0], [150 0]};
% Uniform Circular Array - hald-wavelength inter-antenna spacing
numOfAe = 5;
angleDiff = 360 / numOfAe;
radius = 0.5 / sin(deg2rad(angleDiff*2));
array = zeros(numOfAe,3);
for Ae_ith = 1:1:numOfAe
    array(Ae_ith, 1:1:2) = radius * [cos(deg2rad(30+angleDiff*(Ae_ith-1))) sin(deg2rad(30+(Ae_ith-1)*angleDiff))];
end
SNR = [0 40];

%% Received signals with different SNR values
% symbolsIn_Rec = zeros(length(symbolsOut)+max(cell2mat(delays)), length(SNR));

figure();
for snr_ith = 1:1:length(SNR)
    % The output is symbolsIn is Fx1 complex vector
    symbolsIn=fChannel(paths,symbolsOut.',delays,betas,DOAs,SNR(snr_ith),array);
    
    % the gold-sequence of desired signal
    i = 1;
    goldseq = as(:,i);
    
    % channel parameter (i.e., delay) estimation
    [delay_estimate, DOA_estimate] = fChannelEstimation(symbolsIn,goldseq, paths(i), array);
    disp("When the SNR = " + SNR(snr_ith) + " dB, the estimated delay is " + sprintf("%d, %d, and %d", delay_estimate));
    disp("The estimated DOA is " + sprintf("[%d, %d]", DOA_estimate(1), DOA_estimate(2)));
    % obtain the DOA of two MAI signals
    [delay_estimate2, DOA_estimate2] = fChannelEstimation(symbolsIn,as(:,2), paths(2), array);
    [delay_estimate3, DOA_estimate3] = fChannelEstimation(symbolsIn,as(:,3), paths(3), array);
    
    % Superresolution subspace beamformer
    weights = superresolution(array, DOA_estimate, [DOA_estimate; DOA_estimate2; DOA_estimate3]);
    
    % Combine the signals received by different antennas
    symbolsIn = symbolsIn * conj(weights);
    
    % QPSK Demodulation
    bitsOut=fDSQPSKDemodulator(symbolsIn, goldseq, phi, delay_estimate, paths(i), betas.');
    
    % bit error rate
    [~,ber] = biterr(bitsOut,bitsOuts(:,1));
    disp("The Bit Error Rate is " + 100 * ber + "%.");
    
    % Reconstruction of the image with desired image's properties
    subplot(1,2,snr_ith);
    fImageSink(bitsOut,imageBitLength(i,1),xs(i,1),ys(i,1));
    if snr_ith == 1
        title(["The received image when SNR = 0 dB", "(Bit Error Rate is " + sprintf("%.2f", 100*ber) + " %)"], 'fontsize', 16);
    else
        title(["The received image when SNR = 40 dB", "(Bit Error Rate is " + sprintf("%.2f", 100*ber) + " %)"], 'fontsize', 16);
    end
end
