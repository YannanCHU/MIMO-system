% Yannan CHU, GROUP (EE4/MSc), 2010, Imperial College.
% 2022, Jan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs channel estimation for the desired source using the received signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% symbolsIn (Fx1 Complex) = R channel symbol chips received
% goldseq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence of the desired source used in the modulation process
% paths (Integer) = number of paths of the desired signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% delay_estimate = Vector of estimates of the delays of each path of the
% desired signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delay_estimate]=fChannelEstimation(symbolsIn,goldseq,paths)
F = length(symbolsIn);
% W = Nc = length of gold-sequence
goldseq = 1 - 2 * goldseq;
W = length(goldseq);
% The number of channel symbol chips after DS-QPSK Modulation
R = floor(F / W) * W;

maxPossibleDelay = F - R;

% The correlation between symbolsIn and goldseq
Rsg = zeros(maxPossibleDelay+1, 1);

for delay_ith = 0:1:maxPossibleDelay
    symbolsIn_ith = symbolsIn(1+delay_ith: R+delay_ith, 1);
    Rsg(delay_ith+1) = sum(abs(goldseq' * reshape(symbolsIn_ith.', W, R/W)));
end

[~, index] = maxk(Rsg, paths);

delay_estimate = mod(index-1, W);

% the fading coefficients are assumed to be known in advance
% the DOA are useless and their values have no meaning.


