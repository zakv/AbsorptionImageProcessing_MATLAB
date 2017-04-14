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