% Yannan CHU, GROUP (EE4/MSc), 2010, Imperial College.
% 2022, Jan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs channel estimation for the desired source using the received signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% symbolsIn (FxAeNum Complex) = R channel symbol chips received
% AeNum is the number of antenna in receiver
% goldseq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence of the desired source used in the modulation process
% paths (Integer) = number of paths of the desired signal
% array (AeNumx3) = the position of each antenna in the array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% delay_estimate = Vector of estimates of the delays of each path of the
% desired signal
% DOA_estimate = Estimates of the azimuth and elevation of each path of the
% desired signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delay_estimate, DOA_estimate]=fChannelEstimation(symbolsIn,goldseq,paths,array)
F = length(symbolsIn);

% AeNum is the number of antenna in receiver - 5
AeNum = size(symbolsIn,2);
% W = Nc = length of gold-sequence
goldseq = 1 - 2 * goldseq;
W = length(goldseq);

% The number of channel symbol chips after DS-QPSK Modulation
R = floor(F / W) * W;

% If the array received signal vector is discretised by a chip rate
% sampler, i.e., Ts = Tc, number of chips (Nchip)
Nchip = W;
q = 1;

% After extension. Each column vector is transformed to a discretised
% signal and then vectorised. symbolsInMatrix is (N Next) by L
Next = 2  * q * Nchip;  % a scaler - 30
symbolsInMatrix = symbolExtension(symbolsIn, AeNum, Nchip, q);

% number of snapshots - L
L = size(symbolsInMatrix,2);

% The corvariance matrix
Rxx = (symbolsInMatrix * symbolsInMatrix') / L;

% Source Detection
% Number of sources - M
M = sourceDetection(Rxx, L);

% MUSIC Algorithm
[eigenVec,eigenVal] = eig(Rxx);
[~, pos] = maxk(abs(diag(eigenVal)), M);

signalEigVec = eigenVec(:, pos);

% Projection on signal and noise sub-space
Ps = fpo(signalEigVec);
Pn = eye(size(Ps)) - Ps;

% possible DOA
azimuths = 0:1:180; 
elevation = 0;

maxPossibleDelay = F - R;

% STAR subspace cost function
Zp = zeros(length(azimuths), maxPossibleDelay);

% The shift matrix J which models relative delay
J = [zeros(1, Next-1), 0; eye(Next-1, Next-1), zeros(Next-1,1)];

% double gold-sequence by padding zeros
goldseqDouble = [goldseq; zeros(W,1)];

for azimuth_ith = azimuths
    for delay_jth = 0:1:maxPossibleDelay
        % The array manifold response
        Sij = spv(array, [azimuth_ith elevation]);
        % The Basic STAR manifold response
        hij =  kron(Sij, (J^delay_jth) * goldseqDouble);
        Zp(azimuth_ith+1, delay_jth+1) = 1 / (hij' * Pn * hij);
    end
end

% find a few maximum gains
[delay_candidates, azimuth_index] = max(abs(Zp));
[~, delay_estimate] = maxk(delay_candidates, paths);

delay_estimate = sort(delay_estimate);

azimuth_estimate = zeros(1, paths);

for path_ith = 1:1:paths
    [~, azimuth_estimate(path_ith)] = max(Zp(:, delay_estimate(path_ith)));
end

azimuth_estimate = azimuth_estimate.' - 1;
delay_estimate = mod(delay_estimate.' - 1, W);

elevation_estimate = zeros(size(azimuth_estimate));

DOA_estimate = [azimuth_estimate, elevation_estimate];

% the fading coefficients are assumed to be known in advance


