function save_image(image_in,filename)
%Saves the data in image_in to a .ascii csv file
%   image_in should be a 2D array of data
%
%   filename should be the full filename to be written including the path
%   and the file extension '.ascii'

dlmwrite(filename,image_in,'delimiter','\t')
end

