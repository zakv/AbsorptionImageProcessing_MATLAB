function save_image(image_in,filename)
%Saves the data in image_in to a tab-delimited .ascii csv file
%   === Inputs ===
%   image_in should be a 2D array of data
%
%   filename should be the full filename to be written including the path
%   and the file extension '.ascii'
%
%   === Example Usage ===
%   >> image = load_image( fullfile('20170412','Savefile_45_back.ascii') ); %load an image
%   >> save_image(image,'save_image_test.ascii'); %Writes data to a file

dlmwrite(filename,image_in,'delimiter','\t')
end

