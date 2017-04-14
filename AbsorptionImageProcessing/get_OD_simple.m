function [ OD_simple ] = get_OD_simple( image_name )
%Uses our previous method to get the optical depth
%   image_name should be a string with the name of the file without
%   '_raw.ascii' or any ending like that
%
%   This function simply subtracts the _noise from the _raw and _back files
%   and then takes the log of their ratios to get the optical depth.  It
%   does not do anything clever to reconstruct the background so it does
%   not handle vibrations or lab temperature drifts well.  In particular
%   some fringes in the image do not go away even after extensive
%   averaging.
%
%   Example Usage:
%   >> OD_simple = get_OD_simple( fullfile('20170406','Savefile_1') );
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