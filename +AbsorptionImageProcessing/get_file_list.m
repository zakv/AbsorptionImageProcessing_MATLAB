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