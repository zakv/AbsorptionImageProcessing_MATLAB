function [ residual, projection ] = get_residual_eig( image_in, basis, mean_back, back_region )
%Removes the projection of image_in onto basis using eigenfaces
%   image_in should be an array containing the pixel counts
%
%   basis should be a basis array made by make_basis_eig()
%
%   mean_back should be the mean of the background images (with zeros wherever
%   back_region is zero).  This is analogous to <\Gamma> mentioned in the
%   paper below.
%
%   back_region should the same back_region that was given to
%   make_basis_eig() in order to generate the basis.  It should be a 2D
%   array with 1's in the pixels that should be considered as background
%   and used, and 0's in the pixels that should be ignored (e.g. if there
%   are atoms there).  The easiest way to make this matrix is to use
%   make_back_region().
%
%   residual is the image with the background removed
%
%   projection is the constructed background image that was used
%
%
%   Note that this function properly accounts for the effects of ignoring the
%   regions with atoms as long as the basis used was constructed with make_basis_eig()
%   and the same back_region that was used to make the basis is used here
%
%   Note that this function uses the eigenfaces algorithm, as demonstrated
%   in "Reduction of interference fringes in absorption imaging of cold atom
%   cloud using eigenface method" by Li et. al.

%In accordance with the eigenface algorithm we first subtract off the mean face.
%We'll keep image_in 2D so that we can use get_residual_cheng() below
size_in=size(image_in);
mean_back=reshape(mean_back,size_in);
image_in=image_in-mean_back;

%get_residual_cheng() works almost exactly the same way so we'll re-use the code
% if nargin<=3
%     [residual,projection] = get_residual_cheng(image_in,basis);
% elseif nargin==4
%     [residual,projection] = get_residual_cheng(image_in,basis,back_region);
% end
% 
%Projection method mentioned in paper (replaced by get_residual_cheng() above)
% if nargin<3
%     [n_rows, n_cols] = size(image_in);
%     back_region=ones(n_rows,n_cols);
% end
% 
% %Do projection and subtract the pojection from the image
% used_region=back_region.*image_in;
% 
% used_region=used_region(:);
% 
% n_vectors=size(basis,2);
% projection=zeros(size(basis,1),1); %preallocate
% for k=1:n_vectors
%     u_i=basis(:,k);
%     projection=projection + (u_i'*used_region)*u_i;
% end
% residual=image_in(:)-projection;

%Same projection as in get_residual_cheng()
used_region=back_region.*image_in;
used_region=used_region(:);
projection=basis*(transpose(basis)*used_region);
residual=image_in(:)-projection;

%Return images to original shapes
projection=reshape(projection,size_in);
residual=reshape(residual,size_in);

%Now we add back in the mean face
projection=projection+mean_back;

end