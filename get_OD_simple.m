function [ OD_simple ] = get_OD_simple( image_name )
%Uses our previous method to get the optical depth
%   image_name should be a string with the name of the file without '_raw.ascii'
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