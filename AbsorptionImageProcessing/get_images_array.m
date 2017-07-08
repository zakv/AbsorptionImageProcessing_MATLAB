function [ images_array, image_size ] = get_images_array( file_list, image_size, show_progress )
%Given a list of filenames opens the files and turns each into a column of
%images_array
%   === Inputs ===
%   file_list should be a linear cell array with a filename in each cell.
%   It is best to get file_list from get_file_list()
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
%   will pop up displaying the progress of the basis creation. This
%   defaults to false.
%
%   === Outputs ===
%   images_array is a matrix where each column is a basis image. It is NOT
%   orthogonal or even normalized.  An image P is flattened into a column
%   vector by writing P(:).  Any column can be returned to it's original
%   matrix form using Matlab's reshape() function.
%
%   image_size is a vector giving the dimensions of the 2D images.  It is
%   equivalent to the output of size(image_in) where image_in is one of the
%   images in its inital 2D form.
%
%   === Example Usage ===
%   >> ls_pattern = fullfile('20170405','Savefile_4*_back.ascii');
%   >> file_list = get_file_list(ls_pattern);
%   >> images_array = get_images_array(file_list);

image_size_given=true;
if nargin<2
    image_size_given=false; %The input image_size was not provided
end
if nargin<3
    show_progress=false; %default show_progress to be false
end

%Arrange image data into an array (columns form a non-orthogonal basis)
n_files=length(file_list);
k_max=n_files;
if show_progress
    waitbar_handle = waitbar(0,sprintf('Starting...\n')); %creates a progress bar window
    set(findall(waitbar_handle,'type','text'),'Interpreter','none'); %Don't interpret underscores as subscripts
end

%Figure out how many pixels there are in each image
if image_size_given
    %Can calculate the number of pixels in an image from its size
    n_pixels=image_size(1)*image_size(2);
else
    %Figure out Matrix sizes if they weren't given
    %Assume they're all the size of the first one
    data=load_image(file_list{1});
    image_size=size(data);
    n_pixels=numel(data);
end

%Pre-allocate basis array
images_array=zeros(n_pixels,k_max); %Each column will be a basis vector, which is an image reshaped into a column vector
n_used_images=0; %Keep track of how many images were used so we can trim the array at the end if necessary

%Load all of the images
for k=1:k_max
    current_filename=file_list{k};
    if show_progress %Update waitbar
        waitbar_text=strcat(sprintf('Loading Image #%d of %d\n %s',k,k_max,current_filename));
        waitbar(k/k_max,waitbar_handle,waitbar_text);
    end
    %Get data for the current image
    data=load_image(current_filename);
    if image_size_given
        %Make sure the image is the right size before trying to append it
        if isequal( image_size, size(data) )
            %In this case, the image is the correct size
            %Re-arrange data into column vector and append it to basis matrix
            n_used_images=n_used_images+1;
            images_array(:,n_used_images)=data(:);
        end
    else
        %If image_size was not given, just blindly append the data
        %Re-arrange data into column vector and append it to basis matrix
        images_array(:,k)=data(:);
    end
end

if image_size_given
    %Trim the images_array to the correct size given the number of images
    %that were used
    images_array=images_array(:,1:n_used_images);
end


if show_progress
    close(waitbar_handle)
end

end