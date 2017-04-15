function [ OD ] = get_OD_eig( image_in, basis, mean_back, back_region )
%Returns the OD of image_in using get_residual_eig()
%   image_in should be an array containing the pixel counts
%
%   basis should be a basis array made by make_basis_eig()
%
%   mean_back should be the mean of the background images (with zeros wherever
%   back_region is zero).  This is analogous to <\Gamma> mentioned in the
%   paper.
%
%   back_region should the same back_region that was given to
%   make_basis_eig() in order to generate the basis.  It should be a 2D
%   array with 1's in the pixels that should be considered as background
%   and used, and 0's in the pixels that should be ignored (e.g. if there
%   are atoms there).  The easiest way to make this matrix is to use
%   make_back_region().
%
%   Note that this function properly accounts for the effects of ignoring the
%   regions with atoms as long as the basis used was constructed with make_basis_eig()
%   and the same back_region that was used to make the basis is used here
%
%   Note that this function uses the eigenfaces algorithm, as demonstrated
%   in "Reduction of interference fringes in absorption imaging of cold atom
%   cloud using eigenface method" by Li et. al.

%Get the re-construced background and take the log of the image ratios to
%get the optical depth
[~,projection] = get_residual_eig(image_in,basis,mean_back,back_region);
OD=-log(image_in./projection);
end