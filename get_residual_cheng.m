function [ residual, projection ] = get_residual_cheng( image_in, basis, back_region )
%Removes the projection of image_in onto basis
%   image_in should be an array containing the pixel counts
%
%   basis should be a basis array made by make_basis()
%
%   back_region should be a 2d array with 1's in the pixels that should be
%   considered as background and used, and 0's in the pixels that should be
%   ignored (e.g. if there are atoms there).  The easiest way to make this
%   matrix is to use make_back_region().  This is an optional argument.  If
%   it is not provided, the entire image will be used.
%
%   residual is the image with the background removed
%
%   projection is the constructed background image that was used
%
%
%   Note that this function properly accounts for the effects of ignoring the
%   regions with atoms as long as the basis used was constructed with make_basis_cheng()
%   and the same back_region that was used to make the basis is used here

%If back_region is not provided, default to using the entire image
%This implies that back_region should be all 1's
if nargin<3
    [n_rows, n_cols] = size(image_in);
    back_region=ones(n_rows,n_cols);
end

%Do projection and subtract the pojection from the image
original_shape=size(image_in);
used_region=back_region.*image_in;

used_region=used_region(:);
%This next line works for orthogonal basis vectors. Parentheses make evaulation 
%~1000 times faster.  Without them matlab calculate basis*transpose(basis) first
%which creates a massive array.  With the parentheses matlab calculates
%traspose(basis)*used_region first which outputs a vector, easing the memory
%requirements and number of multiplications/additions in the operation.
projection=basis*(transpose(basis)*used_region);
residual=image_in(:)-projection;
projection=reshape(projection,original_shape);
residual=reshape(residual,original_shape);
end