% Yannan, GROUP (EE4/MSc), 2010, Imperial College.
% 2022, Jan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes MDL to detect the number of sources
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% Rxx (2*N*Nc x 2*N*Nc complex matrix) = The covariance of received signals
% after extension.
% L (Integer) = The number of snapshots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% M (Integer) = The number of detected sources. If the array received
% signal vector is transformed to a discretised signal, the M may be larger
% than the real number of sources.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = sourceDetection(Rxx, L)
[eigenVec,eigenVal] = eig(Rxx);
eigenVal = sort(abs(diag(eigenVal)), 'descend');
N = length(eigenVal);

aic = zeros(1, N);
mdl = zeros(1, N);

for k = 0:N-1
    geometric_mean = 1;
    arithmetic_mean = 0;
    
    for ll = k+1:N
        d_l = eigenVal(ll);
        geometric_mean = geometric_mean * d_l ^ (1/(N-k));
        arithmetic_mean = arithmetic_mean + (1/(N-k)) * d_l;
    end
    aic(1, k+1) = -2 * ((N-k)*L) * log((geometric_mean / arithmetic_mean))...
        + 2*k*(2*N-k);
    mdl(1, k+1) = -1 * ((N-k)*L) * log((geometric_mean / arithmetic_mean))...
        + 0.5*k*(2*N-k)*log(L);
end

% Number of detected sources:
[~, aic_M] = min(aic);  aic_M = aic_M - 1;
[~, mdl_M] = min(mdl);  mdl_M = mdl_M - 1;
% disp("The number of sources detected by AIC is " + aic_M);
% disp("The number of sources detected by MDL is " + mdl_M);

% Similar results can be achieved by either AIC or MDL. But MDL is used.
M = mdl_M;