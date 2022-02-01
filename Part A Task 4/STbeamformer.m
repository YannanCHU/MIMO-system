% Yannan, GROUP (EE4/MSc), 2010, Imperial College.
% 2022, Jan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement the Spatiotemporal Beaformer to generate the weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% array (Nx3) the antenna position in uniform circular array. N is the 
% number of Rx antennas in the array.
% delay_estimate (Lpx1 Integer) Delays of multiple paths.
% Lp (Integer) is the number of paths.
% DOA_estimate (Lpx2 Integer) DOAs of multiple paths.
% beta (Lpx1 complex) Fading coefficients of each channel path.
% goldSeq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence of the desired signal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% weights (2*N*Ncx1 complex) the weights for signals received by different antennas.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function weights = STbeamformer(array, delay_estimate, DOA_estimate, beta, goldSeq, Lp)
goldSeq = 1-2*goldSeq;

Nc = length(goldSeq);
N = size(array,1);

% ST manifold matrix
Hi = zeros(2*N*Nc, Lp);

Next = 2 * Nc;

% The shift matrix J which models relative delay
J = [zeros(1, Next-1), 0; eye(Next-1, Next-1), zeros(Next-1,1)];

% double gold-sequence by padding zeros
goldseqDouble = [goldSeq; zeros(Nc,1)];

for path_ith = 1:1:Lp
    Si = spv(array, DOA_estimate(path_ith, :));
    Hi(:, path_ith) = kron(Si, (J^delay_estimate(path_ith)) * goldseqDouble);
end

weights = Hi * beta;
