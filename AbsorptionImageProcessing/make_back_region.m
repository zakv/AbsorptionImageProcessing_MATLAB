function [ back_region ] = make_back_region( image_in, row_min, row_max, col_min ,col_max)
%Creates an array that marks a background region around a rectangular atom
%region
%   === Inputs ===
%   image_in should be an image for which you'd like to mark the background
%   region.  It is only used to determine how many rows and columns
%   back_region should have
%
%   row_min, etc. should be indices that correspond to the region of the
%   image which should be ignored due to the possible presence of atoms.
%   The rows/columns specified by row_min, row_max, etc. are included in
%   the atom region.
%
%   === Outputs ===
%
%   === Notes ===
%   This function is only designed for an image in which you'd like to
%   ignore a single rectangular region.  The background region of the array
%   contains 1's while the region to be ignored contains 0's.  More
%   complicated background regions can be created by the user.  Simply
%   start with the array ones(size(image_in)) and then mark all the pixels
%   you'd like to ignore (due the presence of atoms) with a zero.  The
%   algorithms here do not have any requirements on the size or shape of
%   the atom region, as long as there is some portion of the image is
%   background to analyze.
%
%   === Example Usage ===
%   >> %Get an image we'd like to analyze in the future
%   >> %(we need to know what dimensions we need to make back_region)
%   >> filename = fullfile('20170405','Cool100d100d80PGCZ4.4_1_raw.ascii');
%   >> image_in = load_image(filename);
%   >> 
%   >> %Suppose the atoms are in rows 60:70 and columns 50:80
%   >> back_region=make_back_region(image_in,61,70,50,80);

[n_rows, n_cols] = size(image_in);
back_region=ones(n_rows,n_cols); %Start with an array of ones
back_region(row_min:row_max,col_min:col_max)=0; %Set the appropriate region to 0
end