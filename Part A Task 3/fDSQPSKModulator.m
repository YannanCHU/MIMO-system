% Yannan CHU, GROUP (EE4/MSc), 2010, Imperial College.
% 2022, Jan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform DS-QPSK Modulation on a vector of bits using a gold sequence
% with channel symbols set by a phase phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% bitsIn (Px1 Integers) = P bits of 1's and 0's to be modulated
% goldseq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence to be used in the modulation process
% phi (Integer) = Angle index in radians of the QPSK constellation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% symbolsOut (Rx1 Complex) = R channel symbol chips after DS-QPSK Modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [symbolsOut]=fDSQPSKModulator(bitsIn,goldseq,phi)
% because QPSK is used, the number of total bits should be a even number
if mod(length(bitsIn),2) == 1
    bitsIn = [bitsIn; 0];
end

% Number of symbols (One symbol consists of two bits)
symbolNumber = length(bitsIn) / 2;
symbols = zeros(symbolNumber, 1);

% Convert the gold-sequence from 0, 1 bits to -1, 1 bit stream
goldseq = 1 - 2 * goldseq;

% QPSK modulation
for i = 1:symbolNumber
    twoBits = bitsIn(2*i-1: 2*i,1);
    if isequal(twoBits, [0; 0])
        symbols(i,1) = sqrt(2)*(cos(phi) + 1i * sin(phi));
    elseif isequal(twoBits, [0; 1])
        symbols(i,1) = sqrt(2)*(cos(phi + pi/2) + 1i * sin(phi + pi/2));
    elseif isequal(twoBits, [1; 1])
        symbols(i,1) = sqrt(2)*(cos(phi + pi) + 1i * sin(phi + pi));
    elseif isequal(twoBits, [1; 0])
        symbols(i,1) = sqrt(2)*(cos(phi + 3 * pi/2) + 1i * sin(phi + 3 * pi/2));
    end
end

symbolsOut = symbols * goldseq';
R = symbolNumber * length(goldseq);
symbolsOut = reshape(symbolsOut.', R, 1);
end