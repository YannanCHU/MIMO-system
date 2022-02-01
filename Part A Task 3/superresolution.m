% Yannan, GROUP (EE4/MSc), 2010, Imperial College.
% 2022, Jan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Superresolution subspace beamformer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% array (Nx3) the antenna position in uniform circular array. N is the 
% number of Rx antennas in the array.
% desiredDirection (1x2) the DOA of the desired signal.
% detectedDirections (Mx2) DOAs of all signals including the desired signal
% and interference signals. M is the number of sources.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% ws (Nx1 complex) the weights for signals received by different antennas.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ws = superresolution(array, desiredDirection, detectedDirections)
    Sd = spv(array, desiredDirection);
    Sj = spv(array, setdiff(detectedDirections, desiredDirection, 'rows'));
    %finds the Projection operator
    ws = (eye(size(fpo(Sj))) - fpo(Sj)) * Sd;
end