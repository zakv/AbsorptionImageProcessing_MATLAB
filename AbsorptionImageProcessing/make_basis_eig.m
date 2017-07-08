function [ basis, mean_back, eigenvalues ] = make_basis_eig( file_list, back_region, max_vectors, image_size, show_progress )
%Given a list of filenames opens the files and forms an eigenfaces basis from the data
%   === Inputs ===
%   file_list should be a linear cell array with a filename in each cell.
%   It is best to get file_list from get_file_list()
%
%   back_region should be a 2d array with 1's in the pixels that should be
%   considered as background and used, and 0's in the pixels that should be
%   ignored (e.g. if there are atoms there).  The easiest way to make this
%   matrix is to use make_back_region().  Alternatively one can use
%   back_region=ones( size(image_in) ) to use the full image, where
%   image_in is an image of the appropriate size.
%
%   max_vectors (optional) is the maximum number of basis vectors to make.
%   Less may be used if file_list has less files.
%
%   image_size (optional) should be a 1D array with two elements (i.e. what
%   the size() function would return for a 2D array) specifying [n_rows,
%   n_columns] for the images used to make the basis.  These dimensions
%   should of course be the same as for the images which will be analyzed
%   using this basis.  When image_size is provided, images of different
%   dimensions are ignored when contstructing the basis, which prevents
%   errors from dimension mismatch.  Therefore, this argument is necessary
%   if some of the files in file_list may have different dimensions than
%   others.  When image_size is not provided, no checks are done on the
%   dimensions of the background images which will lead to errors if they
%   are not all the same
%
%   show_progress (optional) should be true or false.  If true, a window
%   will pop up displaying the progress of loading the images from the hard
%   drive. This defaults to false
%
%   === Outputs ===
%   basis is a 2D arary giving an orthonormal basis that is constructed
%   using the eigenfaces algorithm with slight modifications.  An image P
%   is flattened into a column vector by writing P(:), and such a column
%   vector can be returned to its original matrix form using Matlab's
%   reshape() function.  Each column of basis is an eigenvector of the
%   covariance matrix of the background images.  This implies that each
%   column of basis is a basis vector for the subspace of variations from
%   the mean background image mean_back. The columns of basis are listed in
%   order of decreasing eigenvalue. If there are more than max_vectors in
%   the basis, the ones with the largest eigenvalues are kept and the rest
%   are discarded.  Each column of basis is a basis image analagous to the
%   eigenvectors U_i from the paper mentioned below.  The difference is how
%   the normalization and orthogonalization is defined.  The U_i are
%   orthonormal in the usual sense.  However, if vec1 and vec2 are columns
%   of basis, then they are orthonormal in the sense that
%   (back_region.*vec1) has magnitude 1 and
%   (back_region.*vec1)*(back_region.*vec2)=0.  In this sense the columns
%   of basis are orthonormal over the region specified by back_region. This
%   normalization is used to properly account for the effects of ignoring
%   the atom region of the image.
%
%   mean_back is the mean of the background images.  This is analogous to
%   $\Psi$ mentioned in the paper mentioned below.
%
%   eigenvalues is a list of the eigenvalues corresponding to the elements
%   of the basis (the first eigenvalue is for the first eigenvector etc.).
%   These are not necessary for further calculation but may be helpful for
%   deciding how many vectors to keep
%
%   === Notes ===
%   This function uses the eigenfaces algorithm, as demonstrated in
%   "Reduction of interference fringes in absorption imaging of cold atom
%   cloud using eigenface method" by Li et. al.  In that paper, the atom
%   region is ignored by setting all its pixels to zero before constructing
%   the background.  Simply doing this messes up the algorithm.  This is
%   clear when one considers how the algorithm would treat an image that
%   has no atoms.  For example, imagine plugging in one of the basis
%   vectors e_1 into the algorithm. We know that the projection onto e_1
%   should be 1, however if part of the input is set to zero before taking
%   the inner product, the result will be less than one, and the
%   reconstructed background will not be as accurate as it should be.
%   
%   To remedy this, we set the atom region to 0 in the background images
%   before constructing an orthonormal basis with them.  Next we figure out
%   how to write that orthonormal basis as a superposition of those masked
%   background images. We then take that superposition of input images, but
%   this time without masking them, and that creates our basis.  This basis
%   is not orthonormal in the sense that B'*B=1, rather it's orthonormal in
%   the sense it forms an orthonormal basis for the subspace spanned by the
%   background images in the region where back_region is 1.  Specifically,
%   if vec1 and vec2 are columns of basis, then they are orthonormal in the
%   sense that (back_region.*vec1) has magnitude 1 and
%   (back_region.*vec1)*(back_region.*vec2)=0.  One can then use this basis
%   to reconstruct the back_region of image_in, and that gives the
%   background of the full image.  This can then be subtracted from
%   image_in to give the picture of the atom cloud. That projection and
%   subtraction is performed in separate functions; this function simply
%   contructs the basis that those other functions require.
%
%   === Example Usage ===
%   >> %Get an image we'd like to analyze in the future
%   >> %(we need to know what dimensions we need to make back_region)
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
%   >> [basis_eig,mean_back] = make_basis_eig(file_list,back_region,max_vectors);
%   >>
%   >> %Let's plot a few to see what we've made
%   >> image_size = size(image_in);
%   >> mean_back_reshaped = reshape(mean_back,image_size);
%   >> plot_image(mean_back_reshaped,'Mean background image');
%   >> plot_basis_vectors(basis_eig,1:3,image_size);

%Call get_images_array() with the proper number of arguments
%This gives us the data in the files as a matrix with each column containing
%the data from one image
if nargin<3
    images_array=get_images_array(file_list);
    max_vectors=Inf;
elseif nargin==3
    images_array=get_images_array(file_list);
elseif nargin==4
    images_array=get_images_array(file_list,image_size);
elseif nargin==5
    images_array=get_images_array(file_list,image_size,show_progress);
end

%Build the A matrix from the paper (but with the atom region masked)
n_vectors=size(images_array,2);
back_region=back_region(:); %make into column vector
A=zeros(size(images_array)); %preallocate
mean_back=mean(images_array,2);
for k=1:n_vectors
    A(:,k)=back_region.*(images_array(:,k)-mean_back);
end

%This is the clever, numerically-efficient step of eigenfaces

%This gives the eigenvectors with the largest values
n_eigs=min(n_vectors,max_vectors); %eigs() complains if n_eigs>n_vectors
[V,D]=eigs(A'*A,n_eigs,'lm');
%These are not necessarily sorted from largest to smallest eigenvalue, let's do that
eigenvalues=diag(D);
[~,ind]=sort(abs(eigenvalues),'descend'); %ind gives us the indices of the eigenvalues large->small
eigenvalues=eigenvalues(ind); %Now eigenvalues are sorted by magnitude
V=V(:,ind); %Now V is sorted with the eigenvectors of decreasing eigenvalues.

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