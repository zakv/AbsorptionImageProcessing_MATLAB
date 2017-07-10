function [ file_list ] = get_recent_file_list( ls_pattern, max_files )
%Returns a list of at most max_files files whose names match ls_pattern
%   === Inputs ===
%   ls_pattern should be a pattern that could be passed to the ls function.
%   It can (and usually should) include the "*" wildcard to match multiple
%   files.
%
%   max_files is the maximum number of file names to return.  Fewer than
%   max_files may be returned if not enough file names match ls_pattern.
%   If you don't want to limit the number of returned files, use the
%   function get_file_list() instead.
%
%   === Outputs ===
%   file_list is a linear cell array (a column vector) with a filename in
%   each cell.  The list will include the most recently modified files
%   based on the files' modification date returned by dir().  The length of
%   the list will be at most max_files.
%
%   === Notes ===
%   The modification times of the files correspond to when the files were
%   written to the hard drive.  This may not correspond to when the time
%   the images were taken with the camera.  For example, if the files are
%   copied from one computer to another, the modification times will
%   correspond to when each of the files were copied over.  If they are
%   copied in one big batch, they may not be copied written to the hard
%   drive in the order in which they were taken.
%
%   === Example Usage ===
%   >> ls_pattern = fullfile('20170405','Savefile_4*_back.ascii');
%   >> max_files = 100;
%   >> file_list = get_recent_file_list(ls_pattern, max_files);

%First get the list of all files that match ls_pattern
file_list=get_file_list(ls_pattern);

%Now let's get all of the modification times
%There's probably a way to vectorize this, but I coudn't figure it out
n_files_in=length(file_list);
mod_times=zeros(1,n_files_in);
for k=1:n_files_in
    dir_output=dir( file_list{k} );
    mod_times(k)=dir_output.datenum;
end

%Now let's figure out the ordering of the modification times
[~, sorting_indices]=sort(mod_times,'descend');

%Now rearrange the files in file_list so that they are sorted
file_list=file_list(sorting_indices);

%Now let's trim file_list to be at most max_files long
file_list=file_list( 1:min(end,max_files) ); %Take up to 50 most recent files

end