classdef AbsorptionImageProcessor
%Library for functions related to processing data from absorption imaging
%   Detailed explanation goes here
%   TODO: add basic instructions here
    
methods(Static)
    function import_functions()
        %Imports the functions defined in this class into the caller's
        %namespace
        %   This allows the functions to be called directly rather than via
        %   the dot notation.  In other words, entering
        %   'AbsorptionImageProcessor.import_functions() at the command
        %   line allows users to write 'plot_image(...)' instead of the
        %   lengthier 'AbsorptionImageProcessor.plot_image(...)'
        %
        %   A convenient alternative to using this function would be to
        %   enter 'aip=AbsorptionImageProcessor' at the command line, then
        %   functions can be called like so 'aip.plot_image(...)'
        
        %First get the names of all of the methods in this class file
        methods_list=methods(AbsorptionImageProcessor);
        for k=1:length(methods_list)
            method_name=methods_list{k};
            %Now that we have the names let's make function handles
            func=str2func(['AbsorptionImageProcessor.',method_name]);
            %Now let's put those function handles in the caller's namespace
            assignin('base',method_name,func);
        end
    end


function [ file_list ] = get_file_list( ls_pattern )
%Takes a file name pattern with wildcards and returns a list of file names
%   This function used to have separate code for windows and unix systems that made the
%   unix version try to imitate the behavior of "ls" on windows.  Now both use dir instead
%
%   file_list should be a linear cell array with a filename in each cell.
%   It is best to get file_list from get_file_list()

% The "ls" function works differently on linux vs windows so we have separate code for each
% if ispc()
%     %Windows code
%     file_list=ls(ls_pattern);
% end
if isunix()  || ispc()
    temp=dir(ls_pattern); %cells of objects, one object for each file
    k_max=size(temp,1); %k_max= number of files
    file_list=cell(k_max,1);
    for k=1:k_max
        file_list{k}=fullfile(temp(k).folder,temp(k).name);
    end
    %     %"ls" section for linux compatability
    %     temp=dir(ls_pattern); %cells of objects, one object for each file
    %     k_max=size(temp,1); %k_max= number of files
    %     longest_length=0; %Need to figure out longest filename's length
    %     for k=1:k_max
    %         longest_length=max(longest_length,size(temp(k).name,2));
    %     end
    %     file_list=zeros(k_max,longest_length);
    %     for k=1:k_max
    %         name_length=size(temp(k).name,2);
    %         file_list(k,1:name_length)=temp(k).name;
    %     end
    %     file_list=char(file_list); %This should be the same thing as "ls" returns on windows
end
end



function [ file_list ] = indices_to_file_list( filename_pattern, k_list )
%Returns a cell array of filenames with the k_list values substituted
%   filename_pattern is a string to be passed to sprintf that contains
%   a formatting operator (e.g. '%d' for integers)
%
%   k_list is the list of values to substitute into filename_pattern
%
%   file_list is a cell array of filenames created by substituting the
%   values in k_list into filename_pattern using sprintf
file_list=cell(length(k_list),1);
l=1;
for k=k_list
    filename=sprintf(filename_pattern,k);
    file_list{l}=filename;
    l=l+1;
end
end



function [ fig ] = plot_image( image_in, plot_name, clim )
%Plots image_in with a colorbar
%   image_in should be a 2D array
%   plot_name (optional) should be a string for the plot title
%   clim (optional) is a 2-element [cmin, cmax] giving the upper and lower
%       limits of the color scale
fig=figure();
image(image_in,'CDataMapping','scaled');
if nargin>1
    title(plot_name,'Interpreter','none');
end
colorbar();
if nargin>2
    set(gca,'clim',clim);
end
end



function plot_basis_vectors( basis, indices )
%Plots image_in with a colorbar
%   basis should be a 2D array of basis vectors from get_basis()
%   indices should be a 1D array of indices
for index=indices
    figure();
    vector=basis(:,index);
    %There's probably a smarter way to do the next line, but I'm not sure
    %exactly what that way is right now
    image_in=reshape(vector,121,121);
    image(image_in,'CDataMapping','scaled');
    plot_name=sprintf('Basis Vector #%d',index);
    title(plot_name,'Interpreter','none');
    colorbar();
end
end



function [ fig ] = plot_cross_sections( image_in, DisplayName, limits, fig )
%Plots the integrated horizontal and vertical cross sections of image_in
%   image_in should be a 2D array
%
%   DisplayName should be a string giving the name for the legend
%
%   limits should be [row_min,row_max,col_min,col_max] where each refers
%   to the limis of the horizontal or vertical integration directions.
%   For example, the vertical cross section will be the sum of the rows
%   where only indices row_min:row_max are included in the integration
%
%   fig is an optional argument.  If desired the data can be plotted on
%   a figure previously output by this function.  It will be assumed
%   that the figure has the same subplot structure as this function
%   outputs.  If fig is not provided, a new figure will be created and
%   returned

%open a new figure if necessary
if nargin<4
    fig=figure();
end

%Unpack limits
row_min=limits(1);
row_max=limits(2);
col_min=limits(3);
col_max=limits(4);

%Start with horizontal cross section (i.e. Integrate columns)
subplot(2,1,1);
hold on
title('Horizontal Cross Section');
plot(sum(image_in(row_min:row_max,:),1),'DisplayName',DisplayName);
legend('show');
hold off

%Vertical cross section (i.e. Integrate rows)
subplot(2,1,2);
hold on
title('Vertical Cross Section');
plot(sum(image_in(:,col_min:col_max),2),'DisplayName',DisplayName);
hold off
end



function [ OD_simple ] = get_OD_simple( image_name )
%Uses our previous method to get the optical depth
%   image_name should be a string with the name of the file without '_raw.ascii'
filename_raw=strcat(image_name,'_raw.ascii');%
data_raw=importdata(filename_raw);
filename_back=strcat(image_name,'_back.ascii');%
data_back=importdata(filename_back);
filename_noise=strcat(image_name,'_noise.ascii');%
data_noise=importdata(filename_noise);

dataI= data_raw-data_noise;
dataII= data_back-data_noise;

OD_simple = -1*log(abs(dataI)./abs(dataII));
end



function [ images_array, image_size ] = get_images_array( file_list, show_progress )
%Given a list of filenames opens the files and turns each into a column of images_array
%   file_list should be a linear cell array with a filename in each cell.
%   It is best to get file_list from get_file_list()
%
%   show_progress is optional and should be true or false.  If true, a
%   window will pop up displaying the progress of the basis creation. This
%   defaults to false
%
%   images_array is a matrix where each column is a basis image. It is
%   NOT orthogonal or even normalized.  An image P is flattened into a column
%   vector by writing P(:).  Any column can be returned to it's original matrix
%   form using Matlab's reshape() function.
%
%   image_size is a vector giving the dimensions of the 2D images.  It is equivalent
%   to the output of size(image_in) where image_in is one of the images in its inital 2D
%   form.

if nargin<2
    show_progress=false; %default show_progress to be false
end

%Arrange image data into an array (columns form a non-orthogonal basis)
n_files=length(file_list);
k_max=n_files;
if show_progress
    waitbar_handle = waitbar(0,sprintf('Starting...\n')); %creates a progress bar window
    set(findall(waitbar_handle,'type','text'),'Interpreter','none'); %Don't interpret underscores as subscripts
end

%Figure out Matrix sizes and pre-allocate basis array
data=load_image(file_list{1});
image_size=size(data);
images_array=zeros(numel(data),k_max); %Each column will be a basis vector, which is an image reshaped into a column vector
for k=1:k_max
    current_filename=file_list{k};
    if show_progress %Update waitbar
        waitbar_text=strcat(sprintf('Loading Image #%d of %d\n %s',k,k_max,current_filename));
        waitbar(k/k_max,waitbar_handle,waitbar_text);
    end
    %Get data for the current image
    data=load_image(current_filename);
    %Re-arrange data into column vector and append it to basis matrix
    images_array(:,k)=data(:);
end
if show_progress
    close(waitbar_handle)
end
end



function [ mean_image ] = get_mean_image( file_list )
%Averages the data in the images from file_list and returns it as a 2D array
%   file_list should be a linear cell array with a filename in each cell.
%   It is best to get file_list from get_file_list()
%
%   mean_image is the average of all of the images in file_list.  It is
%   returned as a 2D array with the same dimensions as the input images.

%get the images as series of column vectors and get their initial dimensions
[images_array,image_size]=get_images_array(file_list);
mean_image=mean(images_array,2); %average the different images
mean_image=reshape(mean_image,image_size); %return the image to its 2D form
end



function [ back_region ] = make_back_region( image_in, row_min, row_max, col_min ,col_max)
%Creates an array that marks a background region for an image for a reactangular atom region
%   image_in should be an image for which you'd like to mark the background
%   region.  It is only used to determine how many rows and columns back_region
%   should have
%
%   row_min, etc. should be indices that correspond to the region of the image
%   which should be ignored due to the possible presence of atoms
%
%   This function only works for an image in which you'd like to ignore a single
%   rectangular region.  The background region of the array contains 1's while
%   the region to be ignored contains 0's
%
%   Example: >> back_region=make_back_region(image_in,61,70,1,121);
[n_rows, n_cols] = size(image_in);
back_region=ones(n_rows,n_cols); %Start with an array of ones
back_region(row_min:row_max,col_min:col_max)=0; %Set the appropriate region to 0
end



function [ basis ] = make_basis_cheng( file_list, back_region, max_vectors, show_progress )
%Given a list of filenames opens the files and forms a basis from the data
%   file_list should be a cell array of file names
%
%   back_region should be a 2d array with 1's in the pixels that should be
%   considered as background and used, and 0's in the pixels that should be
%   ignored (e.g. if there are atoms there).  The easiest way to make this
%   matrix is to use make_back_region().  This is an optional argument.  If
%   it is not provided, the entire image will be used.
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
if nargin==1
    images_array=get_images_array(file_list);
    %make default back_region.  Note that we're just going to make it
    %a column vector, so we don't need to worry about making it a 2D array
    %that's the same size as the images
    back_region=ones( size(images_array,1),1 );
    max_vectors=Inf;
elseif nargin==2
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


function [ OD ] = get_OD_cheng( image_in, basis, back_region )
%Returns the OD of image_in using get_residual_cheng()
%   image_in should be an array containing the pixel counts
%
%   basis should be a basis array made by make_basis()
%
%   back_region should be a 2d array with 1's in the pixels that should be
%   considered as background and used, and 0's in the pixels that should be
%   ignored (e.g. if there are atoms there).  The easiest way to make this
%   matrix is to use make_back_region().  This is an optional argument.  If
%   it is not provided, the entire image will be used.

%Let get_residual_cheng() determine default back_region
if nargin<3
    [~,projection] = get_residual_cheng(image_in,basis);
else
    [~,projection] = get_residual_cheng(image_in,basis,back_region);
end
OD=-log(image_in./projection);
end



function [ basis, mean_back, eigenvalues ] = make_basis_eig( file_list, back_region, max_vectors, show_progress )
%Given a list of filenames opens the files and forms an eigenfaces basis from the data
%   file_list should be a cell array of file names
%
%   back_region should be a 2d array with 1's in the pixels that should be
%   considered as background and used, and 0's in the pixels that should be
%   ignored (e.g. if there are atoms there).  The easiest way to make this
%   matrix is to use make_back_region().  This is an optional argument.  If
%   it is not provided, the entire image will be used.
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
if nargin==1
    images_array=get_images_array(file_list);
    %make default back_region.  Note that we're just going to make it
    %a column vector, so we don't need to worry about making it a 2D array
    %that's the same size as the images
    back_region=ones( size(images_array,1),1 );
    max_vectors=Inf;
elseif nargin==2
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
if nargin<3
    [n_rows, n_cols] = size(image_in);
    back_region=ones(n_rows,n_cols);
end
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
%   back_region should be a 2d array with 1's in the pixels that should be
%   considered as background and used, and 0's in the pixels that should be
%   ignored (e.g. if there are atoms there).  The easiest way to make this
%   matrix is to use make_back_region().  This is an optional argument.  If
%   it is not provided, the entire image will be used.
%
%
%   Note that this function properly accounts for the effects of ignoring the
%   regions with atoms as long as the basis used was constructed with make_basis_eig()
%   and the same back_region that was used to make the basis is used here
%
%   Note that this function uses the eigenfaces algorithm, as demonstrated
%   in "Reduction of interference fringes in absorption imaging of cold atom
%   cloud using eigenface method" by Li et. al.

%Let get_residual_eig() determine default back_region
if nargin<4
    [~,projection] = get_residual_eig(image_in,basis,mean_back);
else
    [~,projection] = get_residual_eig(image_in,basis,mean_back,back_region);
end
OD=-log(image_in./projection);
end
end %End of Static Methods
    
end %End of Class