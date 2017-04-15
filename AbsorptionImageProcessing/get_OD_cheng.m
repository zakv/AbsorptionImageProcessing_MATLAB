function [ OD ] = get_OD_cheng( image_in, basis, back_region )
%Returns the OD of image_in using get_residual_cheng()
%   image_in should be an array containing the pixel counts
%
%   basis should be a basis array made by make_basis()
%
%   back_region should the same back_region that was given to
%   make_basis_cheng() in order to generate the basis.  It should be a 2D
%   array with 1's in the pixels that should be considered as background
%   and used, and 0's in the pixels that should be ignored (e.g. if there
%   are atoms there).  The easiest way to make this matrix is to use
%   make_back_region().

%Get the re-construced background and take the log of the image ratios to
%get the optical depth
[~,projection] = get_residual_cheng(image_in,basis,back_region);
OD=-log(image_in./projection);
end