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
% *** the selection of d is very important, d = 5 can also provide balanced
% gold sequence but cannot perfectly estimate all delaies and recover the image. 
% While d = 6, 8 ,9 can achieve perfect image reconstruction. ***
% d = 9;

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
        disp("Delay is " + k);
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
paths = [3; 1; 1];
% delay of each path
delays = {[mod(X+Y, 4); 4+mod(X+Y, 5); 9+mod(X+Y, 6)]; 8; 13};
% complex fading coefficient of each path
betas = {[0.8; 0.4*exp(-1j*40*pi/180); 0.8*exp(1j*80*pi/180)]; 0.5; 0.2};
% Due to SISO working mode, this DOA information is meaningless actually
DOAs = {[30, 0; 45, 0; 20, 0], [80 0], [150 0]};
% Since there is no antenna array but just a single antenna
array = [0,0,0];
SNR = [0 40];

% Received signals with different SNR values
symbolsIn_Rec = zeros(length(symbolsOut)+max(cell2mat(delays)), length(SNR));

for snr_ith = 1:1:length(SNR)
    % The output is symbolsIn is Fx1 complex vector
    symbolsIn=fChannel(paths,symbolsOut.',delays,betas,DOAs,SNR(snr_ith),array);
    symbolsIn_Rec(:, snr_ith) = symbolsIn;
end

figure();
for snr_ith = 1:1:length(SNR)
    symbolsIn = symbolsIn_Rec(:, snr_ith);
    
    % the gold-sequence of desired signal
    i = 1;
    goldseq = as(:,i);
    
    % channel parameter (i.e., delay) estimation
    [delay_estimate] = fChannelEstimation(symbolsIn,goldseq, paths(i));
%     delay_estimate = [0 7 13]
    disp("When the SNR = " + SNR(snr_ith) + " dB, the estimated delays are " + sprintf("%d, %d, and %d", delay_estimate));
    
    % map the estimated delay to the correct fading coefficient
    beta_estimate = zeros(paths(i), 1);
    for j = 1:1:paths(i)
        delay_index = find(delays{i} == delay_estimate(j));
        if isempty(delay_index)
            beta_estimate(j) = 0;
        else
            beta_estimate(j) = betas{i}(delay_index);
        end
    end

%     beta_estimate = zeros(paths(i), 1);
%     for j = 1:1:paths(i)
%         delay_diff = abs(delay_estimate(j) - delays{i});
%         [~, delay_index] = min(delay_diff);
%         beta_estimate(j) = betas{i}(delay_index);
%     end
    
    % QPSK Demodulation
    bitsOut=fDSQPSKDemodulator(symbolsIn, goldseq, phi, delay_estimate, paths(i), beta_estimate);
    
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
