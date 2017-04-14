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