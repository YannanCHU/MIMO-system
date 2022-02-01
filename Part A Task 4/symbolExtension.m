% Yannan, GROUP (EE4/MSc), 2010, Imperial College.
% 2022, Jan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discritization and Extension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% symbolsIn (FxAeNum Complex) = R channel symbol chips received
% AeNum (Integer) is the number of antenna in receiver
% Nchip (Integer) is code period
% q (Integer) is the oversampling factor
% Outputs
% symbolsInMatrix ((N Next)xL Complex), where N = AeNum, Next = 2*q*Nchip


function symbolsInMatrix = symbolExtension(symbolsIn, AeNum, Nchip, q)
    if mod(length(symbolsIn), Nchip) ~= 0
        symbolsIn = [symbolsIn; zeros(Nchip-mod(length(symbolsIn), Nchip), size(symbolsIn,2))];
    end
    
    % Number of snapshots
    L = length(symbolsIn) / Nchip - 1;
    
    symbolsInMatrix = zeros(2 * Nchip * q * AeNum, L);
    
    for rx_ae_ith = 1:1:AeNum
        for symbol_jth = 1:1:L
            symbolsInMatrix(1+2*Nchip*(rx_ae_ith-1): 2*Nchip*rx_ae_ith, symbol_jth) = ...
                symbolsIn(1+Nchip*(symbol_jth-1):Nchip*(symbol_jth+1), rx_ae_ith);
        end
    end