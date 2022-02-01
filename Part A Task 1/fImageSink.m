% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the received image by converting bits back into R, B and G
% matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% bitsIn (Px1 Integers) = P demodulated bits of 1's and 0's
% Q (Integer) = Number of bits in the image
% x (Integer) = Number of pixels in image in x dimension
% y (Integer) = Number of pixels in image in y dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function fImageSink(bitsIn,Q,x,y)
% number of bit per pixel in each layer
bitPerPixel = 8;

% extract the first meaningful bits from the input bit-stream
bits_vector = bitsIn(1:Q, 1);

% reshape the meaning bits from vector form to matrix form
bits_matrix = reshape(bits_vector, length(bits_vector)/bitPerPixel, bitPerPixel);
bits_matrix = uint8(bits_matrix);

img = reshape((bi2de(bits_matrix)), y, x, 3);
% figure();
imshow(img);
