function [ basis, mean_back, eigenvalues ] = make_basis_eig( file_list, back_region, max_vectors, show_progress )
%Given a list of filenames opens the files and forms an eigenfaces basis from the data
%   file_list should be a cell array of file names
%
%   back_region should be a 2d array with 1's in the pixels that should be
%   considered as background and used, and 0's in the pixels that should be
%   ignored (e.g. if there are atoms there).  The easiest way to make this
%   matrix is to use make_back_region().  Alternatively one can use
%   back_region=ones( size(image_in) ) to use the full image, where
%   image_in is an image of the appropriate size.
%
%   max_vectors is the maximum number of basis vectors to use and is an
%   optional argument.  Less may be used if file_list has less files.
%
%   show_progress is optional and should be true or false.  If true, a
%   window will pop up displaying the progress of the basis creation. This
%   defaults to false
%
%   file_list should be a linear cell array with a filename in each cell.
%   It is best to get file_list from get_file_list()
%
%   basis is an orthonormal basis (normalized over the region specified by,
%   but orthogonal over the full image back_region)
%   that is constructed using the eigenfaces algorithm. This is analagous to
%   the matrix B from the paper mentioned below. Each column is a basis image. An
%   image P is flattened into a column vector by writing P(:).  Any column
%   can be returned to it's original matrix form using Matlab's reshape()
%   function.  The columns are eigenvectors listed in order of decreasing
%   eigenvalues.  If more there are more than max_vectors in the basis, the ones
%   with the largest eigenvalues are kept and the rest are discarded.
%
%   mean_back is the mean of the background images (with zeros wherever
%   back_region is zero).  This is analogous to <\Gamma> mentioned in the
%   paper below.
%
%   eigenvalues is a list of the eigenvalues corresponding to the elements of
%   the basis (the first eigenvalue is for the first eigenvector etc.).  These
%   are not necessary for further calculation but may be helpful for deciding
%   how many vectors to keep
%
%   Note that this function uses the eigenfaces algorithm, as demonstrated
%   in "Reduction of interference fringes in absorption imaging of cold atom
%   cloud using eigenface method" by Li et. al.  In that paper, the atom region
%   is ignored by setting all its pixels to zero before constructing the
%   background.  Simply doing this messes up the algorithm.  This is clear when
%   you consider how the algorithm would treat an image that has no atoms.  For
%   example, imagine plugging in one of the basis vectors e_1 into the algorithm.
%   We know that the projection onto e_1 should be 1, however if part of the
%   input is set to zero before taking the inner product, the result will be
%   less than one, and the reconstructed background will not be as accurate as it
%   should be.
%   To remedy this, we set the atom region to 0 in the background images before
%   constructing an orthonormal basis with them.  Next we figure out how to write
%   that orthonormal basis as a superposition of those masked background images.
%   We then take that superposition of input images, but this time without the
%   mask, and that creates our basis.  This basis is not orthonormal in the sense
%   that B'*B=1, rather it's orthonormal in the sense it forms an orthonormal basis
%   for the subspace spanned by the background images in the region where
%   back_region is 1.
%   One can then use this basis to reconstruct the back_region of image_in, and this
%   this immediately gives us the background of the full image.  This is then
%   subtracted from image_in to give the picture of the atom cloud. This projection
%   and subtraction is performed in a separate function

%Call get_images_array() with the proper number of arguments
%This gives us the data in the files as a matrix with each column containing
%the data from one image
if nargin<3
    images_array=get_images_array(file_list);
    max_vectors=Inf;
elseif nargin==3
    images_array=get_images_array(file_list);
elseif nargin==4
    images_array=get_images_array(file_list,show_progress);
end

% n_vectors=size(images_array,2);
% %Set atom region to zero by multiplying each column with back_region
% back_region=back_region(:); %make into column vector
% images_array_masked=zeros(size(images_array)); %preallocate
% for k=1:n_vectors
%     images_array_masked(:,k)=back_region.*images_array(:,k);
% end
% 
% %Figure out mean face (analog of the mean <\Gamma> in the paper)
% A=zeros(size(images_array_masked)); %preallocate
% mean_back=mean(images_array_masked,2);
% for k=1:n_vectors
%     A(:,k)=images_array_masked(:,k)-mean_back;
% end

%Build the A matrix from the paper (but with the atom region masked)
n_vectors=size(images_array,2);
back_region=back_region(:); %make into column vector
A=zeros(size(images_array)); %preallocate
mean_back=mean(images_array,2);
% mean_back_masked=back_region.*mean_back;
for k=1:n_vectors
    A(:,k)=back_region.*(images_array(:,k)-mean_back);
end

%This is the clever, numerically-efficient step of eigenfaces

% [V,D]=eig(A'*A); %V is matrix with eigenvectors as columns, D is diagonal with eigenvalues
% %These are not necessarily sorted from largest to smallest eigenvalue, let's do that
% [eigenvalues,ind]=sort(abs(diag(D)),'descend'); %ind gives us the indices of the eigenvalues large->small
% V=V(:,ind); %Now V is sorted with the eigenvectors of decreasing eigenvalues.
% U=A*V; %The columns of U are the eigenfaces, except with the masked sections set to zero.
% if max_vectors<n_vectors %We may as well trim down our basis now
%     U=U(:,1:max_vectors);
% end

%This gives the eigenvectors with the largest values
%However they are ordered with the largest eigenvalues to the right, so we'll flip things around
%for consistency with the previous work
n_eigs=min(n_vectors,max_vectors); %eigs() complains if n_eigs>n_vectors
[V,D]=eigs(A'*A,n_eigs,'lm');
% eigenvalues=fliplr(diag(D)); %D is diagonal with eigenvalues. fliplr reverses order
% V=fliplr(V);
%These are not necessarily sorted from largest to smallest eigenvalue, let's do that
eigenvalues=diag(D);
[~,ind]=sort(abs(eigenvalues),'descend'); %ind gives us the indices of the eigenvalues large->small
eigenvalues=eigenvalues(ind); %Now eigenvalues are sorted by magnitude
V=V(:,ind); %Now V is sorted with the eigenvectors of decreasing eigenvalues.
% U=A*V; %The columns of U are the eigenfaces, except with the masked sections set to zero.
% if max_vectors<n_vectors %We may as well trim down our basis now
%     U=U(:,1:max_vectors);
% end

%The eigenfaces U_i are superpositions of Phi_i as can be seen by
% the relation U_i = A*V_i = [Phi_1,...Phi_M]*V_i = Phi_1*V_i1 +...+ Phi_M*V_iM
%Where V_ij is element number j of vector V_i.
%Now that we know how to write U_i in terms of Phi_i, let's do it using the unmasked Phi_i
A=zeros(size(images_array)); %preallocate
for k=1:n_vectors
    A(:,k)=images_array(:,k)-mean_back;
end
basis=A*V; %This columns of this basis are analogous to the U_i in the paper
for k=1:n_eigs
    vec=basis(:,k);
    %normalize it over the background region
    %i.e. so that (back_region.*vec) has magnitude 1
    basis(:,k)=vec/norm(back_region.*vec);
end

end