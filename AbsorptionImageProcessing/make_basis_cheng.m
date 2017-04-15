function [ basis ] = make_basis_cheng( file_list, back_region, max_vectors, show_progress )
%Given a list of filenames opens the files and forms a basis from the data
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
%   This function makes an orthonormal basis for the images ignoring the
%   regions for which back_region is 0
%   The output basis is a matrix where each column is a basis image. An
%   image P is flattened into a column vector by writing P(:).  Any column
%   can be returned to it's original matrix form using Matlab's reshape()
%   function.
%   Note that the singular value decomposition used orders the entries in
%   its S matrix from largest to smallest, which means that if you're going
%   to only take N of these vectors for your basis, you're best off taking
%   the first N columns of the output basis

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