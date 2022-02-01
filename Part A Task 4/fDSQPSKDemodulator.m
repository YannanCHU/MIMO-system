% Yannan, GROUP (EE4/MSc), 2010, Imperial College.
% 2022, Jan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform demodulation of the received data using <STAR subspace receiver>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% symbolsIn (P/2x1 Complex) = R channel symbol chips received.
% phi (Integer) = Angle index in radians of the QPSK constellation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% bitsOut (Px1 Integers) = P demodulated bits of 1's and 0's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bitsOut]=fDSQPSKDemodulator(symbolsIn,phi)

% The length of symbols
SymbolNumber = length(symbolsIn);

% adjust the first angle of QPSK constellation to pi/4
symbolsIn = symbolsIn*exp(1j*(pi/4-phi));

% Each angle corresponds to a pair of bits
phi1 = pi / 4;
angleRef = [phi1; phi1+pi/2; phi1-pi; phi1-pi/2];

% measure the angle of symbols (these angle ranges from -pi to pi)
angles = angle(symbolsIn);

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
