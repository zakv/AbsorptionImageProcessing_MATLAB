function [ OD ] = get_OD_cheng( image_in, basis, back_region )
%Returns the OD of image_in using get_residual_cheng()
%   image_in should be an array containing the pixel counts
%
%   basis should be a basis array made by make_basis()
%
%   back_region should be a 2d array with 1's in the pixels that should be
%   considered as background and used, and 0's in the pixels that should be
%   ignored (e.g. if there are atoms there).  The easiest way to make this
%   matrix is to use make_back_region().  This is an optional argument.  If
%   it is not provided, the entire image will be used.

%Let get_residual_cheng() determine default back_region
if nargin<3
    [~,projection] = get_residual_cheng(image_in,basis);
else
    [~,projection] = get_residual_cheng(image_in,basis,back_region);
end
OD=-log(image_in./projection);
end