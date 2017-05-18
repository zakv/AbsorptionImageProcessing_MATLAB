function [ OD ] = get_OD_eig( image_in, basis, mean_back, back_region )
%Returns the optical depth OD of image_in using get_residual_eig()
%   === Inputs ===
%   image_in should be the 2D array containing the pixel counts of an
%   image.
%
%   basis should be a basis array made by make_basis_eig().
%
%   mean_back should be the mean of the background images.  This is
%   analogous to $\Psi$ mentioned in the paper mentioned below.
%
%   back_region should the same back_region that was given to
%   make_basis_eig() in order to generate the basis.  It should be a 2D
%   array with 1's in the pixels that should be considered as background
%   and used, and 0's in the pixels that should be ignored (e.g. if there
%   are atoms there).  The easiest way to make this matrix is to use
%   make_back_region().
%
%   === Outputs ===
%   OD is a 2D array containing the optical depth of the cloud for each
%   pixel.
%
%   === Notes ===
%   This function uses the eigenfaces algorithm, as demonstrated in
%   "Reduction of interference fringes in absorption imaging of cold atom
%   cloud using eigenface method" by Li et. al., but with slight
%   modifications. The main difference between their algorithm and the
%   implementation here is that this code properly accounts for the atom
%   region.  See the discussion in the doc string of make_basis_eig() for
%   more details.
%
%   This function properly accounts for the effects of ignoring
%   the regions with atoms as long as the basis used was constructed with
%   make_basis_eig() and the same back_region that was used to make the
%   basis is used here.
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
%   >> [basis_eig, mean_back] = make_basis_eig(file_list,back_region,max_vectors);
%   >> 
%   >> %Use the basis to get the atomic cloud's optical depth
%   >> OD_eig = get_OD_eig(image_in,basis_eig,mean_back,back_region);
%   >> 
%   >> %Plot the results
%   >> plot_image(image_in,['Original ',filename]);
%   >> plot_image(OD_eig,'Eigenfaces Optical Depth');

%Get the re-construced background and take the log of the image ratios to
%get the optical depth
[~,projection] = get_residual_eig(image_in,basis,mean_back,back_region);
OD=-log(image_in./projection);
end