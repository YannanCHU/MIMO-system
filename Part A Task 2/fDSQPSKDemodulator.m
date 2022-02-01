% Yannan, GROUP (EE4/MSc), 2010, Imperial College.
% 2020, Jan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform demodulation of the received data using <INSERT TYPE OF RECEIVER>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% symbolsIn (Fx1 Integers) = R channel symbol chips received
% goldseq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence of the desired signal to be used in the demodulation process
% phi (Integer) = Angle index in radians of the QPSK constellation points
% delay_estimate (Integer) = the estimated delay which is a channel
% parameter
% path (Integer) = the number of path of the desired signal
% beta_estimate (path x 1 complex vector) = the corresponding fading
% coefficients of each path of the desired signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% bitsOut (Px1 Integers) = P demodulated bits of 1's and 0's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bitsOut]=fDSQPSKDemodulator(symbolsIn,GoldSeq,phi,delay_estimate, pathNum, beta_estimate)
% Convert the gold-sequence from 0, 1 bits to -1, 1 bit stream
goldseq = 1 - 2 * GoldSeq;
% Length of gold-sequence - Nc = W
W = length(goldseq);

% The length of received symbols
F = length(symbolsIn);
% The number of channel symbol chips after DS-QPSK Modulation
R = floor(F / W) * W;
SymbolNumber = R/W;

% pathNum = 1;
symbolMulti = zeros(SymbolNumber,pathNum);
for j = 1:1:pathNum
    symbolsIn_ = symbolsIn(1+delay_estimate(j):1:R+delay_estimate(j));
    % after de-spreading the input symbols
    symbolMulti(:, j) = goldseq' * reshape(symbolsIn_.', W, SymbolNumber);
end

% if there are multiple paths, we need a combination rule (MRC) to combine
% signals from different paths.
if pathNum > 1
    weights = conj(beta_estimate);
else
    weights = 1;
end

symbols = symbolMulti * weights;

% adjust the first angle of QPSK constellation to pi/4
symbols = symbols*exp(1j*(pi/4-phi));

% Each angle corresponds to a pair of bits
phi1 = pi/4;
angleRef = [phi1; phi1+pi/2; phi1-pi; phi1-pi/2];

% measure the angle of symbols (these angle ranges from -pi to pi)
angles = angle(symbols);

% QPSK demodulation
P = 2 * SymbolNumber;

bitsOut = zeros(P, 1);
for i = 1:1:SymbolNumber
    % measure the difference between detetced symbol angle and four
    % reference angles
    angleDif = abs(angles(i) - angleRef);
    [~, pos] = min(angleDif);
    if pos == 1
        bitsOut(2*i-1:2*i,1) = [0; 0];
    elseif pos == 2
        bitsOut(2*i-1:2*i,1) = [0; 1];
    elseif pos == 3
        bitsOut(2*i-1:2*i,1) = [1; 1];
    elseif pos == 4
        bitsOut(2*i-1:2*i,1) = [1; 0];
    end 
end
end
