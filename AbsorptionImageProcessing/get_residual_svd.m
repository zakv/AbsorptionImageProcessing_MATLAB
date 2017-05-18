function [ residual, projection ] = get_residual_svd( image_in, basis, back_region )
%Removes the projection of image_in onto basis
%   === Inputs ===
%   image_in should be the 2D array containing the pixel counts of an
%   image.
%
%   basis should be a basis array made by make_basis_svd().
%
%   back_region should the same back_region that was given to
%   make_basis_svd() in order to generate the basis.  It should be a 2D
%   array with 1's in the pixels that should be considered as background
%   and used, and 0's in the pixels that should be ignored (e.g. if there
%   are atoms there).  The easiest way to make this matrix is to use
%   make_back_region().
%
%   === Outputs ===
%   residual is the image with the background removed.
%
%   projection is the constructed background image that was subtracted from
%   image_in to give residual.
%
%   === Notes ===
%   This function properly accounts for the effects of ignoring the regions
%   with atoms as long as the basis used was constructed with
%   make_basis_svd() and the same back_region that was used to make the
%   basis is used here.  See the doc string for make_basis_svd() for more
%   details.
%
%   === Example Usage ===
%   >> %Get the image we'd like to analyze
%   >> filename = fullfile('20170405','Cool100d100d80PGCZ4.4_1_raw.ascii');
%   >> image_in = load_image(filename);
%   >> 
%   >> %Select a background region
%   >> row_min=40; row_max=60; col_min=50; col_max=80;
%   >> back_region = make_back_region(image_in,row_min,row_max,col_min,col_max);
%   >> 
%   >> %Make a basis
%   >> max_vectors = 20; %20 is typically a good number for this
%   >> ls_pattern = fullfile('20170405','*_back.ascii');
%   >> file_list = get_file_list(ls_pattern);
%   >> basis_svd = make_basis_svd(file_list,back_region,max_vectors);
%   >> 
%   >> %Use the basis to get the residual and projection of image_in
%   >> [residual,projection] = get_residual_svd(image_in,basis_svd,back_region);
%   >> 
%   >> %Plot the results
%   >> plot_image(image_in,['Original ',filename]);
%   >> plot_image(residual,'svd Residual');
%   >> plot_image(projection,'svd Projection');

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