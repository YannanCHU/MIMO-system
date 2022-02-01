% Yannan, GROUP (EE4/MSc), 2010, Imperial College.
% CHU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models the channel effects in the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% paths (Mx1 Integers) = Number of paths for each source in the system.
% For example, if 3 sources with 1, 3 and 2 paths respectively then
% paths=[1;3;2]
% symbolsIn (MxR Complex) = Signals being transmitted in the channel
% delay (Cx1 cell) = Delay for each path in the system starting with
% source 1
% beta (Cx1 cell) = Fading Coefficient for each path in the system
% starting with source 1
% DOA (Cx1 cell) = Direction of Arrival for each source in the system in the form
% [Azimuth, Elevation]
% SNR (Integer) = Signal to Noise Ratio in dB
% array (Nx3) = Array locations in half unit wavelength. If no array then should
% be [0,0,0]. N is the number of Rx Antennas.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% symbolsOut (FxN Complex) = F channel symbol chips received from each antenna
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [symbolsOut]=fChannel(paths,symbolsIn,delay,beta,DOA,SNR,array)
seed = 3;       rng(seed);

% Number of sources
M = length(paths);

snr = 10^(SNR/10);

% only 1 receiver
N = 1;
% Due to the delay, F is slightly greater than R.
R = size(symbolsIn, 2);
F = R + max(cell2mat(delay));

symbolsOut = zeros(F, N);

for i = 1:1:M
    % number of paths per source
    pathNumber = paths(i);
    symbols_ith_multi = zeros(pathNumber, F);
    delays_ith_source = delay{i};
    betas_ith_source = beta{i};
    
    for j = 1:1:pathNumber
        % add the delay by shifting the row vector
        % add the fading coefficient (propagation gain) - beta
%         disp("The delay is " + delays_ith_source(j) + ";  beta is " + betas_ith_source(j));
        symbols_ith_multi(j, 1+delays_ith_source(j): R+delays_ith_source(j)) = symbolsIn(i, :);
        symbols_ith_multi(j,:) = betas_ith_source(j) * symbols_ith_multi(j,:);
    end
    
    % Array Manifold Vector - Si = 1 all the time due to the SISO system
    Si = spv(array, DOA{i});
    
    symbols_ith = Si * symbols_ith_multi;
    
    % compute the signal power
    symbolsPower_ith = mean(abs(symbols_ith(:)).^2);
    
    % compute the noise power
    noisePower_ith = symbolsPower_ith / snr;
    
    % A random variable's power equals its mean-squared value - E(X^2)
    % Since the noise is zero-mean AWGN, its power equals to its variance
    % The noise is random complex number and the nosie power is equally 
    % distributed in real and imaginary part.
    noise = sqrt(noisePower_ith/2) * (randn(size(symbols_ith)) + 1i * randn(size(symbols_ith)));
    
    symbols_ith = symbols_ith + noise;
    symbolsOut = symbolsOut + symbols_ith.';
end


