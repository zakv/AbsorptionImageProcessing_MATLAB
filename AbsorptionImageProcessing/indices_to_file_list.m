function [ file_list ] = indices_to_file_list( filename_pattern, k_list )
%Returns a cell array of filenames with the k_list values substituted
%   === Inputs ===
%   filename_pattern is a string to be passed to sprintf that contains
%   a formatting operator (e.g. '%d' for integers)
%
%   k_list is the list of values to substitute into filename_pattern
%
%   === Outputs ===
%   file_list is a cell array of filenames created by substituting the
%   values in k_list into filename_pattern using sprintf
%
%   === Example Usage ===
%   >> filename_pattern=fullfile('20170405','Savefile_%d_back.ascii');
%   >> k_list=1:20;
%   >> file_list=indices_to_file_list(filename_pattern,k_list);

file_list=cell(length(k_list),1);
l=1;
for k=k_list
    filename=sprintf(filename_pattern,k);
    file_list{l}=filename;
    l=l+1;
end
end