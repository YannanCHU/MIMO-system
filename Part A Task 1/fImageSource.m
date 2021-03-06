% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads an image file with AxB pixels and produces a column vector of bits
% of length Q=AxBx3x8 where 3 represents the R, G and B matrices used to
% represent the image and 8 represents an 8 bit integer. If P>Q then
% the vector is padded at the bottom with zeros.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% filename (String) = The file name of the image
% P (Integer) = Number of bits to produce at the output - Should be greater
% than or equal to Q=AxBx3x8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% bitsOut (Px1 Integers) = P bits (1's and 0's) representing the image
% x (Integer) = Number of pixels in image in x dimension
% y (Integer) = Number of pixels in image in y dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [bitsOut,x,y]=fImageSource(filename,P)
% import the image
img = imread(filename);
% figure(); imshow(img);

% number of bit per pixel in each layer
bitPerPixel = 8;

% length of column vector of bits
% Q = numel(img) * bitPerPixel;
[A, B, ~] = size(img);
Q = A * B * 3 * bitPerPixel;

% convert the pixel intensity of RGB image from decimal to binary version
img_binary = de2bi(img);

% generate the bit stream
img_bits = reshape(img_binary, Q, 1);
bitsOut = zeros(P, 1);

% if P > Q, implement the zero padding
bitsOut(1:Q, 1) = img_bits;

y = A;
x = B;
end

