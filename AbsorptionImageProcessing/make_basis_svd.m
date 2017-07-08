function [ basis ] = make_basis_svd( file_list, back_region, max_vectors, image_size, show_progress )
%Given a list of filenames opens the files and forms a basis from the data
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
%   show_progress (optional) and should be true or false.  If true, a
%   window will pop up displaying the progress of loading the images from
%   the hard drive. This defaults to false
%
%   === Outputs ===
%   basis is a 2D arary giving an orthonormal basis that is constructed
%   using  an svd decomposition.  An image P is flattened into a column
%   vector by writing P(:), and such a column vector can be returned to its
%   original matrix form using Matlab's reshape() function.  Each column is
%   a basis vector for the subspace of image backgrounds. The columns of
%   basis are listed in order of decreasing singular value. If there are
%   more than max_vectors in the basis, the ones with the largest singular
%   values are kept and the rest are discarded.  The columns of basis are
%   not orthonormal in the usual way.  If vec1 and vec2 are columns of
%   basis, then they are orthonormal in the sense that (back_region.*vec1)
%   has magnitude 1 and (back_region.*vec1)*(back_region.*vec2)=0.  In this
%   sense the columns of basis are orthonormal over the region specified by
%   back_region. This normalization is used to properly account for the
%   effects of ignoring the atom region of the image.
%
%   === Notes === 
%   It is tempting to believe that simply setting the region of an image
%   with atoms to zero before projecting it onto the subspace of background
%   images will properly account for the atom region.  However, simply
%   doing this messes up the algorithm.  This is clear when one considers
%   how the algorithm would treat an image that has no atoms.  For example,
%   imagine plugging in one of the basis vectors e_1 into the algorithm. We
%   know that the projection onto e_1 should be 1, however if part of the
%   input is set to zero before taking the inner product, the result will
%   be less than one, and the reconstructed background will not be as
%   accurate as it should be.
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
%   >> basis_svd = make_basis_svd(file_list,back_region,max_vectors);
%   >>
%   >> %Let's plot a few to see what we've made
%   >> image_size = size(image_in);
%   >> plot_basis_vectors(basis_svd,1:3,image_size);

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

%Take non-orthogonal basis and make it orthogonal over the background region
n_vectors=size(images_array,2);
%Set atom region to zero by multiplying each column with back_region
images_array_masked=zeros(size(images_array)); %preallocate
back_region=back_region(:); %make into column vector
for k=1:n_vectors
    images_array_masked(:,k)=back_region.*images_array(:,k);
end
%Take svd decomposition of matrix.  The svd decomposition for a matrix
%A is given by A=U*S*V' where S is diagonal and U and V are orthonormal,
%and V' is the transpose of V.  The columns of U will be an orthonormal
%basis for the columnspace of A.  We'll use images_array_masked for A.
%Note that U and V being orthonormal implies U*U'=1 and V*V'=1
%
%Once we have the svd decomposition for the masked basis, we figure out
%How to transform A into U (i.e. how to take superpositions of the columns
%of A to make U).  We then apply this to the unmasked images_array to
%get a basis which is orthonormal over the background region but still
%contains info about the atom region.  To subtract the background (which
%will be done in a different function) we mask the image_in with
%back_region and then project it onto the supspace spanned by the basis
%created here and subtract that from the unmasked image_in
%
%Here's how we figure out how to transform A into U.  We right multiply
%the svd decomposition by V, which gives A*V=U*S.  Now S is diagonal, but
%it may have zeros on the diagonal so it might not be possible to invert.
%The work-around is to note that since S is diagonal, (U*S) is just the
%same as U except that each column is scaled by a scalar.  Since we know
%that the columns of U should be normalized, we just normalize the columns
%of A*V to get the basis.  Note that when I say normalize here I mean
%normalize over the background subspace (i.e. don't include the atom
%region when calculating the norm).

[~,~,V] = svd(images_array_masked,'econ');
basis=images_array*V; %this is our basis, but with non-normalized columns
%Let's trim the basis now that we can, then we'll get back to normalizing
%Due to svd's output ordering, these will be the vectors with the largest
%singular values.
if max_vectors<n_vectors
    basis=basis(:,1:max_vectors);
end
k_max=size(basis,2);
for k=1:k_max
    vec=basis(:,k);
    %normalize it over the background region
    %i.e. so that (back_region.*vec) has magnitude 1
    basis(:,k)=vec/norm(back_region.*vec);
end
end